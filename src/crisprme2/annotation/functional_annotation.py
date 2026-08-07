"""
annotation_transformer.py
-------------------------
Functional-annotation transform: the pipeline hands it a batch of resolved
alignments and it fills the ``u32`` annotation slots in place, exactly the way
the scoring transform fills the score slots.

Per row it does three things:

1. **interval** — derive the target's genomic footprint ``[start, stop)`` from
   the occurrence position, the miner offset, the strand, and the run's PAM
   placement (see :func:`_target_intervals`). This mirrors the sink's
   ``PamContext::target_start`` so annotations line up with the reported
   ``start`` column exactly.
2. **fetch** — for each input BED, ``TabixFile.fetch`` the records overlapping
   that interval and read their 4th-column terms.
3. **encode** — OR each term's bit (via the cached :class:`FeatureRegistry`
   mapping) into the alignment's annotation slot for that BED.

Slot order is BED input order, which is also the order of ``annotation_names``
and therefore of the annotation columns in the report.

Ownership & lifecycle
~~~~~~~~~~~~~~~~~~~~~~~
The transformer owns its annotation resources. Given only the BED paths, its
``__init__`` builds one :class:`FeatureRegistry` per file (scanning column 4)
and opens one :class:`AnnotationBed` (tabix handle; the index is computed on
construction if absent). Those handles are held in an ``ExitStack`` and released
when the transformer is closed, so it is used as a context manager::

    with AnnotationTransformer(bed_paths, pam_len, pam_upstream,
                               contig_ids, loggers) as annotate:
        with Pipeline(..., transforms=[scorer, annotate]) as pipeline:
            pipeline.submit(batcher)   # drains, calling annotate(batch) per chunk

``AnnotationBed`` also closes itself on ``__del__``, so a forgotten ``close()``
is a leak until GC, not a correctness problem.

Concurrency
~~~~~~~~~~~
``pysam.TabixFile`` is not thread-safe. This transform assumes it runs as a
single GIL-serialized ``PyTransform`` stage (one instance, one set of open
BEDs). If the pipeline ever fans annotation across worker threads, each worker
needs its own transformer with its own opened BEDs.
"""

from __future__ import annotations

from collections.abc import Sequence
from contextlib import ExitStack
from types import TracebackType
from typing import Any, Dict

import numpy as np

import os


from ..crisprme_core_api import AlignmentBatch, FeatureRegistry
from ..logger import CrisprmeLoggers
from ..protocol import Transformer

from .bedfile import AnnotationBed
from .crisprme2_annotation_error import Crisprme2FunctionalAnnotationError


# ASCII bytes used while scanning the aligned target / strand columns.
_GAP = 0x2D  # '-' : RNA-bulge gap in rseq (consumes no reference base)
_MINUS = 0x2D  # '-' : reverse-strand marker in the strand column
_NUL = 0x00  # row terminator; bytes past the first NUL are stale

# Maximum annotation columns (== u32 annotation slots on the alignment frame).
_MAX_ANNOTATIONS = 10


def _target_intervals(
    pos: np.ndarray,
    offset: np.ndarray,
    strand: np.ndarray,
    rseq2d: np.ndarray,
    pam_upstream: bool,
    pam_len: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized ``[start, stop)`` for every row of a batch (0-based, half-open).

    ``start`` follows the same rule as the sink's ``target_start``::

        scanned_on_revcomp = is_reverse XOR pam_upstream   # == Strand::scanned_on_revcomp
        start = pos            if scanned_on_revcomp        # offset dropped
              = pos + offset    otherwise

    ``stop = start + gapless_nucs + pam_len``, where ``gapless_nucs`` is the
    protospacer's reference footprint: the count of non-gap bases in the
    NUL-terminated ``rseq`` row (``'-'`` marks RNA bulges, which consume no
    reference base; DNA-bulge bases are present and counted).
    """
    pos = pos.astype(np.int64, copy=False)
    offset = offset.astype(np.int64, copy=False)
    # gapless_nucs: live = positions before the first NUL; of those, count the
    # non-gap bases
    live = np.logical_and.accumulate(rseq2d != _NUL, axis=1)
    gapless = np.count_nonzero((rseq2d != _GAP) & live, axis=1).astype(np.int64)
    scanned_on_revcomp = (strand == _MINUS) ^ pam_upstream
    start = np.where(scanned_on_revcomp, pos, pos + offset)
    stop = start + gapless + pam_len
    return start, stop


class FunctionalAnnotator(Transformer):
    """Fills each alignment's ``u32`` annotation slots from a set of BEDs.

    The transformer builds and owns its own registries and tabix handles; the
    caller supplies only file paths.

    Attributes
    ----------
    annotations : Sequence[str]
        Annotation BED paths (plain or bgzip-compressed), in the order the
        annotation columns should appear in the report (== ``annotation_names``
        order). At most 10 (the number of annotation slots).
    pam_len : int
        Reference length of the run's PAM.
    pam_upstream : bool
        PAM placement for the run (``True`` for Cas12a-style upstream PAMs).
    contig_ids : Sequence[str]
        Dense contig-id -> name table (the one ``search._compute_contig_ids``
        feeds the sink), indexed by a row's ``contig_id``.
    loggers : CrisprmeLoggers
        Shared logger bundle for error routing.
    """

    __slots__ = (
        "_anns",
        "_pam_len",
        "_pam_upstream",
        "_contigs",
        "_loggers",
        "_stack",
        "_rows_seen",
        "_fetch_groups",
    )

    def __init__(
        self,
        annotations: Sequence[str],
        pam_len: int,
        pam_upstream: bool,
        contig_ids: Dict[str, int],
        loggers: CrisprmeLoggers,
    ) -> None:
        self._loggers = loggers
        self._pam_len = pam_len
        self._pam_upstream = pam_upstream
        self._contigs = {i: c for c, i in contig_ids.items()}
        # Build a registry and open a tabix handle per BED. The ExitStack owns
        # every handle so a failure part-way through unwinds the ones already
        # opened instead of leaking them.
        self._stack = ExitStack()
        anns: list[tuple[AnnotationBed, FeatureRegistry]] = []
        try:
            for path in annotations:
                registry = FeatureRegistry(path, loggers)  # scans column 4
                bed = self._stack.enter_context(AnnotationBed(path, loggers))  # opens
                anns.append((bed, registry))
        except Exception as e:
            self._stack.close()
            self._loggers.errorlog.log_raise_exception(
                f"Failed opening annotation files: {e}",
                os.EX_IOERR,
                Crisprme2FunctionalAnnotationError,
            )
        self._anns = tuple(anns)
        # Instrumentation for the interval cache: how many rows asked for an
        # annotation vs how many distinct intervals actually reached tabix.
        self._rows_seen = 0
        self._fetch_groups = 0
        loggers.verboselog.info(
            f"AnnotationTransformer ready with {len(self._anns)} annotation BED(s)"
        )

    # -- resource lifecycle ----------------------------------------------------

    def close(self) -> None:
        if self._rows_seen:
            self._loggers.basiclog.info(
                f"Annotation interval cache: {self._rows_seen} rows -> "
                f"{self._fetch_groups} distinct intervals "
                f"({100.0 * self._fetch_groups / self._rows_seen:.1f}% fetch rate, "
                f"{self._rows_seen / max(self._fetch_groups, 1):.2f}x reuse)"
            )
        self._stack.close()

    def __enter__(self) -> "FunctionalAnnotator":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> bool:
        self.close()
        return False  # never suppress

    # -- transform -------------------------------------------------------------

    def __call__(self, raw_batch: Any) -> None:
        if not self._anns:
            return  # nothing to annotate
        batch = AlignmentBatch(raw_batch, self._loggers)
        n_rows = batch.n_rows
        if n_rows == 0:
            return  # nothing to annotate
        contig_id = batch.contig_id  # (N,) uint16
        pos = batch.pos  # (N,) uint32  (occurrence fwd-left)
        offset = batch.offset  # (N,) uint8   (miner offset)
        strand = batch.strand  # (N,) uint8   (ASCII '+' / '-')
        rseq = batch.rseq  # (N*L,) uint8, flat
        n = int(contig_id.shape[0])
        if n == 0:
            return
        rseq2d = rseq.reshape(n_rows, -1)
        start, stop = _target_intervals(
            pos, offset, strand, rseq2d, self._pam_upstream, self._pam_len
        )
        # Writeable u32 slot per annotation column; assigned (not OR'd) so the
        # result never depends on prior slot contents
        slots = [batch.feature(s) for s in range(len(self._anns))]
        out = [np.zeros(n, dtype=np.uint32) for _ in self._anns]
        contigs = self._contigs
        n_contigs = len(contigs)
        n_ann = len(self._anns)

        # Visit rows in genomic order so rows sharing an interval land next to
        # each other. Ordering is internal only: every result is written back to
        # its original row index `i`, so `out` is bit-identical.
        order = np.lexsort((stop, start, contig_id)).tolist()

        # One-entry cache. `fetch_features` is a pure function of
        # (contig, lo, hi) for a given open BED, so reuse across identical keys
        # is exact; the sort guarantees duplicates are adjacent, so a single
        # slot catches all of them.
        prev_key: tuple[int, int, int] | None = None
        prev_masks: list[int] = []
        for i in order:
            cid = int(contig_id[i])
            if cid >= n_contigs:
                continue  # unmapped contig -> leave this row unannotated
            lo, hi = int(start[i]), int(stop[i])
            self._rows_seen += 1
            key = (cid, lo, hi)
            if key != prev_key:
                contig = contigs[cid]
                prev_masks = [
                    registry.accumulate(bed.fetch_features(contig, lo, hi))
                    for bed, registry in self._anns
                ]
                prev_key = key
                self._fetch_groups += 1
            for s in range(n_ann):
                out[s][i] = prev_masks[s]
        for s in range(len(self._anns)):
            slots[s][:] = out[s]  # bulk copy into the zero-copy Rust buffer
