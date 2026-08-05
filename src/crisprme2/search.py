""" """

from pathlib import Path
from time import time
from typing import Dict, List, Optional, Tuple

import os
import shutil
import tempfile


from .annotation import FunctionalAnnotator
from .crisprme_core_api import (
    init_native_logging,
    partition_offtargets,
    TargetBatcher,
    Pipeline,
    Thresholds,
)
from .crisprme2_inputargs import Crisprme2SearchInputArgs
from .crisprme2_error import Crisprme2SearchError
from .dna_alphabet import reverse_complement
from .fasta import Fasta, read_fasta_files
from .guide import Guide, GuidesList, read_guides
from .logger import CrisprmeLoggers
from .pam import PAM, read_pam, SPCAS9, XCAS9
from .protocol import Transformer
from .scores import CfdScorer
from .utils import TOOLNAME


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

#: Prefix for the hidden, per-run scratch directory created inside outdir to
#: hold the transient intermediate report. Leading dot keeps it out of the way.
_TMP_DIR_PREFIX: str = ".crisprme2_tmp_"

# ==============================================================================
# Internal search helpers
# ==============================================================================


def _compute_report_name(guide: Guide, pam: PAM, outdir: str) -> str:
    return os.path.join(outdir, f"{guide.sequence}_{pam.pam}.tsv")


def _safe_fasta_contig(fasta: Fasta, contig: str, loggers: CrisprmeLoggers) -> str:
    """
    Return the contig name as it appears in an open *fasta* handle,
    normalising "chr"-prefix mismatches between the dict key and pyfaidx.

    Tries *contig* first; if absent, falls back to the normalised
    single-contig name exposed by ``fasta.contig``.  Closes the handle
    and raises before returning if neither name is found.
    """
    c = contig
    if c not in fasta.contigs:
        contig_alt = f"chr{contig}"  # normalized single-contig name from file
        if contig_alt in fasta:
            c = contig_alt
        else:
            fasta.close()  # ensure file is closed before raising exception
            loggers.errorlog.log_raise_exception(
                f"Contig {contig} not found in FASTA {fasta._filepath}",
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


def _compute_chunk_seq_stran_pairs(
    chunk_seq: str, upstream: bool, loggers: CrisprmeLoggers
) -> List[Tuple[str, int]]:
    chunk_seq_rc = reverse_complement(chunk_seq, loggers)
    if upstream:  # force pam downstream the target sequence
        return [(chunk_seq_rc, 1), (chunk_seq, 0)]
    return [(chunk_seq, 1), (chunk_seq_rc, 0)]


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
    upstream: bool,
    loggers: CrisprmeLoggers,
) -> None:
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
            # assign chunk sequence-strand pairs based on PAM position
            for seq, strand in _compute_chunk_seq_stran_pairs(
                chunk_seq, upstream, loggers
            ):
                result = batcher.feed_chunk(
                    contig_id, chunk_start, seq, strand, core_len
                )
                if result.flushed:
                    _submit_and_log(
                        pipeline, batcher, f"contig={contig!r} chunk={i}", loggers
                    )


def _partition_report_names(prefix: str, outdir: str) -> Tuple[str, str]:
    base = os.path.join(outdir, prefix)
    return (f"{base}_primary.tsv", f"{base}_alternative.tsv")


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
    annotations: List[str],
    annotation_names: List[str],
    cluster_dist: int,
    criteria: List[str],
    output_prefix: Optional[str],
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
    # hidden, unique scratch dir inside outdir holds the transient intermediate.
    # The finally guarantees the dir and its contents are removed on success or
    # on any failure during pipeline setup, mining, or partitioning
    tmpdir = tempfile.mkdtemp(prefix=_TMP_DIR_PREFIX, dir=outdir)
    loggers.verboselog.debug(f"staging intermediate report in {tmpdir}")
    try:
        outpath = _compute_report_name(guide, pam, tmpdir)  # intermediate -> hidden dir
        # pipeline: one context for the entire genome run
        with Pipeline.create(
            _PIPELINE_CHUNKS,
            thresholds,
            transforms,
            pam,
            upstream,
            outpath,
            contig_ids,
            annotations,
            annotation_names,
            loggers,
        ) as pipeline:
            for contig, fasta in fastas.items():
                contig_id = contig_ids[contig]
                loggers.verboselog.debug(
                    f"Processing contig {contig!r} "
                    f"(id={contig_id}, threads={threads}, upstream={upstream}, size={size})"
                )
                contig_start = time()  # trace contig processing running time
                try:
                    _process_contig(
                        fasta,
                        batcher,
                        pipeline,
                        contig,
                        contig_id,
                        overlap,
                        size,
                        upstream,
                        loggers,
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
        # Pipeline closed -> the intermediate report is fully flushed,
        # split it into primary and alternative reports
        prefix = output_prefix or f"{guide.sequence}_{pam.pam}"
        primary_path, alternative_path = _partition_report_names(prefix, outdir)
        partition_offtargets(
            outpath, primary_path, alternative_path, criteria, cluster_dist, loggers
        )
    finally:
        # remove the hidden dir and everything in it, whatever happened
        shutil.rmtree(tmpdir, ignore_errors=True)


def _build_pam_and_guides(
    args: Crisprme2SearchInputArgs, loggers: CrisprmeLoggers
) -> Tuple[GuidesList, PAM]:
    """
    Initialise PAM and guide data structures from validated CLI arguments.

    Parameters
    ----------
    args : Crisprme2SearchInputArgs
        Validated argument namespace.
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Returns
    -------
    tuple[GuidesList, PAM]
        ``(guides, pam)`` ready for use in the search pipeline.
    """
    loggers.basiclog.info("Initialising PAM and guide data structures")
    pam = read_pam(args.pam, loggers)
    guides = read_guides(args, loggers)
    loggers.verboselog.debug(f"PAM: {pam} | guides: {len(guides)}")
    return guides, pam


def _build_thresholds(
    args: Crisprme2SearchInputArgs, loggers: CrisprmeLoggers
) -> Thresholds:
    """
    Construct a :class:`~crisprme2.crisprme_core_api.Thresholds` instance
    from validated CLI arguments.

    Parameters
    ----------
    args : Crisprme2SearchInputArgs
        Validated argument namespace.
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Returns
    -------
    Thresholds
        Alignment thresholds for this run.
    """
    loggers.verboselog.debug(
        f"Building Thresholds(max_mm={args.mm}, bdna={args.bdna}, "
        f"brna={args.brna}, max_edit_dist={args.max_edit_dist})"
    )
    return Thresholds(
        max_mm=args.mm,
        max_bdna=args.bdna,
        max_brna=args.brna,
        max_edit_dist=args.max_edit_dist,
        loggers=loggers,
    )


def _build_transforms(
    pam: PAM,
    annotations: List[str],
    contig_ids: Dict[str, int],
    loggers: CrisprmeLoggers,
) -> List[Transformer]:
    transforms: List[Transformer] = []
    # ---- scoring transform
    if pam.cas_system in [SPCAS9, XCAS9]:
        # CFD score + slot 0
        # CFD pam is the last two bases of the PAM sequence
        # For NGG the key is "GG"; for NGA it is "GA", etc.
        transforms.append(CfdScorer(pam=pam.pam, loggers=loggers))

    # ---> add future scorers here <---

    # ---- annotation transform
    if annotations:
        transforms.append(
            FunctionalAnnotator(
                annotations, len(pam), pam.upstream, contig_ids, loggers
            )
        )

    # ---> add gene annotation here <---

    loggers.verboselog.debug(
        "Transform chain assembled: " f"{[type(t).__name__ for t in transforms]}"
    )
    return transforms


def _search_offtargets_reference_genome(
    fasta_files: List[str],
    pam: PAM,
    guide: Guide,
    upstream: bool,
    outdir: str,
    threads: int,
    thresholds: Thresholds,
    annotations: List[str],
    annotation_names: List[str],
    cluster_dist: int,
    criteria: List[str],
    output_prefix: Optional[str],
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
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Raises
    ------
    Crisprme2ScannerError
        On FASTA I/O errors or scanning failures.
    """

    fastas = read_fasta_files(fasta_files, loggers)
    contig_ids = _compute_contig_ids(list(fastas.keys()))
    # initialize transforms
    transforms = _build_transforms(pam, annotations, contig_ids, loggers)
    size = len(guide) + len(pam) + max(thresholds.bdna, thresholds.brna)
    loggers.verboselog.debug(
        "Starting reference-genome/assembly off-target extraction pipeline\n"
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
        annotations,
        annotation_names,
        cluster_dist,
        criteria,
        output_prefix,
        loggers,
    )


# ==============================================================================
# Public API
# ==============================================================================


def execute_offtargets_search(args: Crisprme2SearchInputArgs) -> None:
    """
    Run the full CRISPRme2 complete-search pipeline.

    This is the composition root: it wires CLI arguments to pipeline
    components and delegates execution to specialised modules.  The call
    graph is::

        execute_offtargets_search(args)
            ├── CrisprmeLoggers(args.outdir)
            ├── _build_pam_and_guides(args)     -> GuidesList, PAM
            ├── _build_thresholds(args)         -> Thresholds
            ├── _build_transforms(pam)          -> list[Transformer]
            └── (per guide)
                └── search_offtargets_reference_genome(...)

    Parameters
    ----------
    args : Crisprme2SearchInputArgs
        Fully validated CLI argument namespace produced by
        :func:`~crisprme2.__main__.create_parser_crisprme2`.

    Raises
    ------
    Crisprme2SearchError
        If any component of the search pipeline fails.
    """
    loggers = CrisprmeLoggers(args.outdir)  # initialize loggers
    init_native_logging(loggers)  # initialize rust-level logging
    loggers.basiclog.info(f"Start {TOOLNAME} search")
    # initialize pam and guide objects
    guides, pam = _build_pam_and_guides(args, loggers)
    # initialize thresholds object
    thresholds = _build_thresholds(args, loggers)
    for guide in guides:
        # retrieve candidate off-targets for current guide
        loggers.verboselog.debug(
            f"Starting off-target search for guide {guide.sequence}"
        )
        if args.vcfs:
            # variant and haplotype aware search path (not yet implemented)
            loggers.verboselog.debug(
                "VCF files provided - variant-aware search path "
                "not yet implemented (skipping)"
            )
            continue
        else:
            # reference-only search path
            _search_offtargets_reference_genome(
                args.fastas,
                pam,
                guide,
                args.upstream,
                args.outdir,
                args.threads,
                thresholds,
                args.annotations,
                args.annotation_names,
                args.cluster_dist,
                args.prioritization_criteria,
                args.output_prefix,
                loggers,
            )
