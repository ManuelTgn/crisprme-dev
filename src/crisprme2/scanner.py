""" """

from .crisprme2_error import Crisprme2ScannerError
from .fasta_utils import read_fasta_files
from .utils import flatten_list, OFFTARGETLEN
from .logger import CrisprmeLoggers
from .fasta import Fasta
from .guide import Guide
from .pam import PAM

from ._crisprme2_native import TargetBatcher

from typing import List, Dict, Tuple
from time import time

import sys
import os


# define sequence chunk size
CHUNKSIZE = 10_000_000

# define overlap size between chunks
CHUNKOVERLAP = 29  # 30 - 1 (we extract 30-mers)

# define batch hits
BATCHITS = 1_000_000


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


def extract_targets(
    fastas: Dict[str, Fasta],
    contig_ids: Dict[str, int],
    guide: Guide,
    pam: PAM,
    size: int,
    right: bool,
    threads: int,
    loggers: CrisprmeLoggers,
):

    batcher = TargetBatcher(
        pam.pam, size, right, threads, BATCHITS, 250_000, CHUNKOVERLAP
    )

    otfname = "results/test-run/test.txt"

    print(contig_ids)

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
                    st = batcher.stats()
                    print(
                        "hits_in_batch", st.hits_in_batch, "unique", st.unique_windows
                    )

                    # TODO: remove after debugging
                    if status.flushed:
                        # batcher.flush_and_align_placeholder_tsv(guide.sequence, 4, otfname)
                        batcher.flush_and_align()

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
    return {c: i for i, c in enumerate(contigs)}


def scan_fasta_reference_genome(
    fasta_files: List[str],
    pam: PAM,
    guide: Guide,
    offset: int,
    right: bool,
    threads: int,
    loggers: CrisprmeLoggers,
):
    loggers.verboselog.debug(
        "Follow reference genome-based off-targets extraction pipeline"
    )
    fastas = read_fasta_files(fasta_files, loggers)  # read input fasta files
    contig_ids = _compute_contig_ids(list(fastas.keys()))  # compute contig ids
    # compute off-target size for extraction
    size = _compute_target_size(guide, pam, offset)  # offset is max(bdna, brna)
    loggers.verboselog.debug(f"Off-targets extraction size: {size}")
    # extract targets from reference genome fasta files
    extract_targets(fastas, contig_ids, guide, pam, size, right, threads, loggers)
