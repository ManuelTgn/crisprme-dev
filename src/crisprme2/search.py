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
from .utils import flatten_list, OFFTARGETLEN
from .logger import CrisprmeLoggers
from .fasta import Fasta
from .guide import Guide
from .pam import PAM

from typing import List, Dict, Tuple
from time import time

import sys
import os


# ---------------------------------------------------------------------------
# Chunk geometry constants
# ---------------------------------------------------------------------------

# Number of base-pairs in each FASTA sub-chunk fed to the batcher.
CHUNKSIZE: int = 100_000

# Number of overlapping base-pairs kept between consecutive chunks.
# Must satisfy: CHUNKOVERLAP >= window_size - 1.
# Window size is at most guide(20) + PAM(3) + max_bulge(2) = 25,
# so 29 is a safe conservative default.
CHUNKOVERLAP: int = 29  # updated at runtime to max(size - 1, 29)

# Default pipeline memory-pool chunk count.
_PIPELINE_CHUNKS: int = 10_000


def _safe_fasta_contig(fasta: Fasta, contig: str, loggers: CrisprmeLoggers) -> str:
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
                Crisprme2ScannerError,
            )
    return c


def receive_data(pipeline):
    complete = False
    while not complete:
        complete, res = pipeline.receive()


def _extract_and_align(fasta: Fasta, contig: str, loggers: CrisprmeLoggers):
    with fasta as f:
        # ensure we fec=tch using a reference that exists in the opened handle
        c = _safe_fasta_contig(fasta, contig, loggers)
        sequence = f.fetch(c)  # fecth contig sequence
        seqlen = len(sequence)  # avoid lazy compute
        chunkedseq = sequence.chunk(CHUNKSIZE, CHUNKOVERLAP)


def extract_targets(
    fastas: Dict[str, Fasta],
    contig_ids: Dict[str, int],
    guide: Guide,
    pam: PAM,
    size: int,
    right: bool,
    threads: int,
    loggers: CrisprmeLoggers,
) -> None:
    for contig, fasta in fastas.items():  # iterate over single fasta
        loggers.verboselog.debug(
            f"Scanning contig {contig} for targets (threads = {threads}, right = {right}, size = {size})"
        )
        start = time()  # trace scanner run time on current contig
        contig_id = contig_ids[contig]  # retrieve contig id
        try:  # Fasta.contig normalizes "chr" prefix; dict key are already be normalized
            with fasta as f:
                # ensure we fetch using a reference that exists in the opened handle
                c = _safe_fasta_contig(fasta, contig, loggers)
                sequence = f.fetch(c)  # fetch contig sequence
                seqlen = len(sequence)
                chunkedseq = sequence.chunk(CHUNKSIZE, CHUNKOVERLAP)
                for i, chunkseq in enumerate(chunkedseq):  # iterate over subchunks
                    core_start = (
                        i * CHUNKSIZE
                    )  # compute chunk start (e.g. 0, 10M, etc.)
                    core_len = min(CHUNKSIZE, seqlen - core_start)
                    # initialize batcher data for subsequence
                    chunk_start = 0 if i == 0 else core_start - CHUNKOVERLAP
                    if len(chunkseq) < size:
                        continue
                    # pass subchunk to rust API batcher
                    status = batcher.feed_chunk(
                        contig_id, chunk_start, chunkseq, core_len
                    )
                    # st = batcher.stats()
                    # print(
                    #     "hits_in_batch", st.hits_in_batch, "unique", st.unique_windows
                    # )

                    # TODO: remove after debugging
                    # if status.flushed:
                    # batcher.flush_and_align_placeholder_tsv(guide.sequence, 4, otfname)
                    # batcher.flush_and_align()

                    # 1. ------> pipeline = pipeline()

                    # 2. ------> pipeline.submit(batcher)

                    # 3. ------> trash batcher -> finalize()

                    # complete = False
                    # while not complete:
                    #     complete, result = pipeline.receive()

                    # windows_written, rows_written = batcher.flush_and_align_placeholder_tsv(
                    #     guide.sequence,
                    #     4,                 # max_mm
                    #     str(otfname),       # PathBuf/str ok
                    # )
                    # loggers.verboselog.debug(
                    #     f"Flushed: windows={windows_written}, rows={rows_written}, "
                    #     f"at contig={contig}, chunk={i}, chunk_start={chunk_start}"
                    # )
        except Exception as e:
            # raise to stop the pipeline
            loggers.errorlog.log_raise_exception(
                f"Scanning contig {contig} failed: {e}",
                os.EX_DATAERR,
                Crisprme2ScannerError,
            )
        finally:

            loggers.verboselog.debug(
                f"Contig {contig} scanned in {time() - start:.2f}s"
            )

    # pipeline.close()

    # Final tail flush: you must explicitly flush remaining windows too
    # windows_written, rows_written = batcher.flush_and_align_placeholder_tsv(
    #     guide.sequence,
    #     4,
    #     str(otfname),
    # )
    # batcher.flush_and_align_placeholder_tsv(guide.sequence, 4, otfname)
    batcher.flush_and_align()
    tail = batcher.finalize()  # clears internal state

    # loggers.verboselog.debug(
    #     f"Final flush: windows={windows_written}, rows={rows_written}; "
    #     f"tail stats: hits={tail.hits_in_batch}, unique={tail.unique_windows}"
    # )


def _compute_target_size(guide: Guide, pam: PAM, offset: int) -> int:
    return len(guide) + len(pam) + offset


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


def scan_reference_genome(
    fastas: Dict[str, Fasta],
    contig_ids: Dict[str, int],
    guide: Guide,
    pam: PAM,
    size: int,
    upstream: bool,
    threads: int,
    thresholds: Thresholds,
    transforms: List,
    loggers: CrisprmeLoggers,
) -> None:
    overlap = _compute_overlap(size)
    # build batcher - one per genome run; reset between flushes by Rust
    batcher = TargetBatcher.create(
        pam, guide, size, upstream, overlap, threads, loggers
    )
    loggers.verboselog.debug(
        f"TargetBatcher ready (id={batcher.id}, size={size}, overlap={overlap})"
    )
    # pipeline: one context for the entire genome run
    with Pipeline.create(_PIPELINE_CHUNKS, thresholds, transforms, loggers) as pipeline:
        pass


def search_offtargets_reference_genome(
    fasta_files: List[str],
    pam: PAM,
    guide: Guide,
    upstream: bool,
    threads: int,
    thresholds: Thresholds,
    transforms: List,
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
    size = 30  # TODO: define as constant?
    loggers.verboselog.debug(
        f"Contigs: {list(fastas.keys())}"
        f" | window size: {size}"
        f" | thresholds: {thresholds}"
    )
    # extract targets from reference genome fasta files
    scan_reference_genome(
        fastas,
        contig_ids,
        guide,
        pam,
        size,
        upstream,
        threads,
        thresholds,
        transforms,
        loggers,
    )
