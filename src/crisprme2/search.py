"""
scanner.py
----------
Genome scanning pipeline: FASTA -> target extraction -> alignment.

This module is the top-level orchestrator that connects:

1. :func:`extract_targets`   — iterate over FASTA contigs, feed chunks to a
                               :class:`~crisprme2.crisprme_core_api.TargetBatcher`,
                               and submit full batches to the alignment
                               :class:`~crisprme2.crisprme_core_api.Pipeline`.
2. :func:`scan_fasta_reference_genome` — public entry-point that wires
                               everything together from the CLI/API layer.

Data flow
~~~~~~~~~
::

    FASTA files
        │
        ▼
    read_fasta_files()          one Fasta handle per contig
        │
        ▼  (contig loop)
    fasta.fetch(contig)         raw nucleotide string
        │
        ▼  (chunk loop, CHUNKSIZE = 10 Mbp, CHUNKOVERLAP = size-1)
    TargetBatcher.feed_chunk()  IUPAC encode + PAM scan -> accumulate windows
        │
        ├─ flushed? ──► pipeline.submit(batcher)   transfer batch to GPU pipeline
        │
        ▼  (after all contigs)
    TargetBatcher.finalize()    flush tail + clear internal map
"""

from .crisprme_core_api import TargetBatcher, Pipeline, Thresholds
from .crisprme2_error import Crisprme2SearchError
from .fasta_utils import read_fasta_files
from .logger import CrisprmeLoggers
from .fasta import Fasta
from .guide import Guide
from .pam import PAM
from .protocol import Transformer

from typing import List, Dict, Tuple
from time import time

import os


# ==============================================================================
# Chunk geometry constants
# ==============================================================================

#: Number of base-pairs in each FASTA sub-chunk fed to the batcher.
CHUNKSIZE: int = 100_000

#: Number of overlapping base-pairs kept between consecutive chunks.
#: Must satisfy: CHUNKOVERLAP >= window_size - 1.
#: Window size is at most guide(20) + PAM(3) + max_bulge(2) = 25,
#: so 29 is a safe conservative default.
CHUNKOVERLAP: int = 29  # updated at runtime to max(size - 1, 29)

#: Default pipeline memory-pool chunk count.
_PIPELINE_CHUNKS: int = 10_000


# ==============================================================================
# Internal search helpers
# ==============================================================================


def _safe_fasta_contig(fasta: Fasta, contig: str, loggers: CrisprmeLoggers) -> str:
    """
    Return the contig name as it appears in an open *fasta* handle,
    normalising "chr"-prefix mismatches between the dict key and pyfaidx.

    Tries *contig* first; if absent, falls back to the normalised
    single-contig name exposed by ``fasta.contig``.  Closes the handle
    and raises before returning if neither name is found.
    """
    c = contig
    if c not in fasta:
        contig_alt = fasta.contig  # normalized single-contig name from file
        if contig_alt in fasta:
            c = contig_alt
        else:
            fasta.close()  # ensure file is closed before raising exception
            loggers.errorlog.log_raise_exception(
                f"Contig {contig} not found in FASTA {fasta._filepath} (available: {fasta.contig})",
                os.EX_DATAERR,
                Crisprme2SearchError,
            )
    return c


def _compute_contig_ids(contigs: List[str]) -> Dict[str, int]:
    """Assign a stable integer id to each contig name"""
    return {c: i for i, c in enumerate(contigs)}


def _compute_overlap(size: int) -> int:
    """
    Return the chunk overlap that satisfies the Rust batcher constraint
    ``overlap >= size - 1``.  Always at least ``CHUNKOVERLAP`` so the
    constant is never silently violated.
    """
    return max(size - 1, CHUNKOVERLAP)


def _chunk_sequence(
    fasta: Fasta, contig: str, overlap: int, loggers: CrisprmeLoggers
) -> Tuple[List[str], int]:
    """
    Fetch a contig sequence from an already-open *fasta* handle and split
    it into overlapping sub-chunks.

    .. note::
        This function must be called **inside** an open ``with fasta``
        block (i.e. from :func:`_process_contig`).  The *fasta* parameter
        is the live handle ``fa``, not the outer wrapper.
    """
    c = _safe_fasta_contig(fasta, contig, loggers)
    sequence = fasta.fetch(c)
    return sequence.chunk(CHUNKSIZE, overlap), len(sequence)


def _submit_and_log(
    pipeline: Pipeline, batcher: TargetBatcher, label: str, loggers: CrisprmeLoggers
) -> None:
    """
    Submit *batcher* to *pipeline* and log the action.

    Separating this into a helper keeps the chunk loop readable and gives
    a single point to add metrics / tracing in the future.
    """
    stats = batcher.stats()
    loggers.verboselog.debug(
        f"{label}: submitting batch - "
        f"{stats.hits_in_batch} hits, {stats.unique_windows} unique windows"
    )
    pipeline.submit(batcher)


def _process_contig(
    fasta: Fasta,
    batcher: TargetBatcher,
    pipeline: Pipeline,
    contig: str,
    contig_id: int,
    overlap: int,
    size: int,
    loggers: CrisprmeLoggers,
):
    """
    Open a FASTA handle, chunk its sequence, and feed each chunk to
    *batcher*, submitting to *pipeline* whenever the batch is full.

    This function owns the ``with fasta`` context for one contig.
    It delegates chunking to :func:`_chunk_sequence` and submission
    to :func:`_submit_and_log`.
    """
    with fasta as fa:
        chunk_seqs, seqlen = _chunk_sequence(fa, contig, overlap, loggers)
        for i, chunk_seq in enumerate(chunk_seqs):
            # absolute genomic start of the full chunk (including left overlap for i > 0)
            core_start: int = i * CHUNKSIZE
            core_len: int = min(CHUNKSIZE, seqlen - core_start)
            chunk_start: int = 0 if i == 0 else core_start - overlap
            if len(chunk_seq) < size:
                # chunk too short to contain even one window; skip rather than
                # sending empty work to rust
                loggers.verboselog.debug(
                    f"Contig {contig!r}, chunk {i}: sequence ({len(chunk_seq)} bp) "
                    f"shorter than window size ({size}), skipping"
                )
                continue
            result = batcher.feed_chunk(contig_id, chunk_start, chunk_seq, core_len)
            if result.flushed:
                _submit_and_log(
                    pipeline, batcher, f"contig={contig!r} chunk={i}", loggers
                )


def _scan_reference_genome(
    fastas: Dict[str, Fasta],
    contig_ids: Dict[str, int],
    guide: Guide,
    pam: PAM,
    size: int,
    upstream: bool,
    outdir: str,
    threads: int,
    thresholds: Thresholds,
    transforms: List[Transformer],
    loggers: CrisprmeLoggers,
) -> None:
    """
    Scan every contig in *fastas* for off-target candidates and route
    full batches through the alignment pipeline.

    Manages three nested levels of state:

    - **Pipeline context** — one :class:`Pipeline` for the entire genome
      run, opened once and closed (workers joined) after all contigs.
    - **Contig loop** — each contig is opened, chunked, and fully
      processed before the next contig begins.
    - **Chunk loop** — :data:`CHUNKSIZE`-bp sub-sequences (with overlap)
      are fed to :class:`TargetBatcher` one at a time.  When
      ``feed_chunk`` returns ``flushed=True``, the batch is submitted to
      the pipeline before the next chunk is processed.

    After all contigs are exhausted, a final tail flush submits any
    windows that did not trigger an automatic flush, and
    :meth:`~crisprme2.crisprme_core_api.TargetBatcher.finalize` clears
    the internal Rust map.

    Parameters
    ----------
    fastas : dict[str, Fasta]
        Mapping from normalised contig name to an unopened
        :class:`~crisprme2.fasta.Fasta` handle.
    contig_ids : dict[str, int]
        Mapping from contig name to its integer index.
    guide : Guide
        Guide RNA object; ``.sequence`` forwarded to the batcher.
    pam : PAM
        PAM object; ``.pam`` forwarded to the batcher.
    size : int
        Window extraction width (guide + PAM + bulge offset).
    upstream : bool
        ``True`` if the PAM is 3' of the protospacer (e.g. SpCas9 NGG).
    outdir : str
        Path of the CSV report. Truncated on open.
    threads : int
        Number of parallel scanner threads inside the batcher.
    thresholds : Thresholds
        Alignment thresholds (max mismatches, DNA bulges, RNA bulges)
        forwarded to the pipeline and used at flush time.
    transforms : list[callable]
        Ordered transform callables forming the pipeline's scoring and
        annotation stage chain.
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Raises
    ------
    Crisprme2SearchError
        If any contig scan fails (FASTA I/O, position overflow, etc.).
        The error message includes the contig name and the underlying cause.
    """
    overlap = _compute_overlap(size)
    # build batcher - one per genome run; reset between flushes by Rust
    batcher = TargetBatcher.create(
        pam, guide, size, upstream, overlap, threads, loggers
    )
    loggers.verboselog.debug(
        f"TargetBatcher ready (id={batcher.id}, size={size}, overlap={overlap})"
    )
    # pipeline: one context for the entire genome run
    with Pipeline.create(_PIPELINE_CHUNKS, thresholds, transforms, pam, upstream, outdir, loggers) as pipeline:
        for contig, fasta in fastas.items():
            contig_id = contig_ids[contig]
            loggers.verboselog.debug(
                f"Processing contig {contig!r} "
                f"(id={contig_id}, threads={threads}, upstream={upstream}, size={size})"
            )
            contig_start = time()  # trace contig processing running time
            try:
                _process_contig(
                    fasta, batcher, pipeline, contig, contig_id, overlap, size, loggers
                )
            except Crisprme2SearchError:
                raise  # already formatted; propagate as-is
            except Exception as e:
                loggers.errorlog.log_raise_exception(
                    f"Processing contig {contig!r} failed: {e}",
                    os.EX_DATAERR,
                    Crisprme2SearchError,
                )
            finally:
                loggers.verboselog.debug(
                    f"Contig {contig!r} processed in {time() - contig_start:.2f}s"
                )
        # tail flush: submit whatever remains after the last auto-flush
        tail_stats = batcher.stats()
        if tail_stats.hits_in_batch > 0 or tail_stats.unique_windows > 0:
            _submit_and_log(pipeline, batcher, "tail flush", loggers)
        # finalize clears internal rust states; log what was flushed in the tail
        final_stats = batcher.finalize()
        loggers.basiclog.info(
            f"Processing complete - batcher id = {batcher.id}, "
            f"total chunks={batcher.total_chunks_fed}, "
            f"total flushes={batcher.total_flushes}, "
            f"tail residual: hits={final_stats.hits_in_batch}, "
            f"unique windows={tail_stats.unique_windows}"
        )
    # pipeline.__exit__ signals EOF and joins all worker threads here


# ==============================================================================
# Public API
# ==============================================================================


def search_offtargets_reference_genome(
    fasta_files: List[str],
    pam: PAM,
    guide: Guide,
    upstream: bool,
    outdir: str,
    threads: int,
    thresholds: Thresholds,
    transforms: List[Transformer],
    loggers: CrisprmeLoggers,
) -> None:
    """
    Full reference-genome off-target scanning pipeline.

    Reads FASTA files, computes window size, assigns contig ids, then
    delegates to :func:`extract_targets` which manages the batcher and
    pipeline lifecycle.

    Parameters
    ----------
    fasta_files : list[str]
        Paths to one or more FASTA files (one per chromosome or all-in-one).
    pam : PAM
        Parsed PAM object.
    guide : Guide
        Guide RNA object.
    upstream : bool
        ``True`` if the PAM is 3' of the protospacer (e.g. SpCas9 NGG).
    outdir : str
        Path of the CSV report. Truncated on open.
    threads : int
        Number of parallel scanner threads.
    thresholds : Thresholds
        Alignment thresholds forwarded to the pipeline.
    transforms : list[callable]
        Transform callables forming the pipeline's scoring/annotation chain.
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Raises
    ------
    Crisprme2ScannerError
        On FASTA I/O errors or scanning failures.
    """
    loggers.verboselog.debug(
        "Starting reference-genome/assembly off-target extraction pipeline"
    )
    fastas = read_fasta_files(fasta_files, loggers)
    contig_ids = _compute_contig_ids(list(fastas.keys()))
    size = len(guide) + len(pam) + max(thresholds.bdna, thresholds.brna)
    loggers.verboselog.debug(
        f"Contigs: {list(fastas.keys())}"
        f" | window size: {size}"
        f" | thresholds: {thresholds}"
    )
    # extract targets from reference genome fasta files
    _scan_reference_genome(
        fastas,
        contig_ids,
        guide,
        pam,
        size,
        upstream,
        outdir,
        threads,
        thresholds,
        transforms,
        loggers,
    )
