""" """

from .crisprme2_error import Crisprme2ScannerError
from .fasta_utils import read_fasta_files
from .utils import flatten_list, OFFTARGETLEN
from .logger import CrisprmeLoggers
from .fasta import Fasta
from .guide import Guide
from .pam import PAM

from .target_candidates_scanner_rs import extract_targets_rs

from typing import List, Dict, Tuple
from time import time

import sys
import os


# define sequence chunk size
CHUNKSIZE = 10_000_000

# define overlap size between chunks
CHUNKOVERLAP = 29  # 30 - 1 (we extract 30-mers)


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
    pam: PAM,
    size: int,
    right: bool,
    threads: int,
    loggers: CrisprmeLoggers,
):
    for contig, fasta in fastas.items():  # iterate over single fasta
        loggers.verboselog.debug(
            f"Scanning contig {contig} for targets (threads = {threads}, right = {right}, size = {size})"
        )
        start = time()  # trace scanner run time on current contig
        try:  # Fasta.contig normalizes "chr" prefix; dict key are already be normalized
            with fasta as f:
                # ensure we fetch using a reference that exists in the opened handle
                c = _safe_fasta_contig(fasta, contig, loggers)
                sequence = f.fetch(c)  # fetch contig sequence
                chunkedseq = sequence.chunk(CHUNKSIZE, CHUNKOVERLAP)
                # preallocate target candidates lists
                candidates_chunk: List[Tuple[List[int], List[int]]] = [None] * len(chunkedseq)  # type: ignore
                for i, chunkseq in enumerate(chunkedseq):
                    # extract targets in spwaning threads on each sequence chunk
                    # go down to rust to optimize threads spawning
                    candidates_chunk[i] = extract_targets_rs(
                        chunkseq, pam.pam, size, right, threads
                    )
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
        all_pos = []
        all_strand = []
        for pos, strand in candidates_chunk:
            all_pos.extend(pos)
            all_strand.extend(strand)

        # candidates = flatten_list(candidates_chunk)
        print(f"Number of candidates: {len(all_pos)}")
        print(f"Pos: {(sys.getsizeof(all_pos) + sum(sys.getsizeof(e) for e in all_pos)) / (1024 ** 3)}")
        print(f"strand: {(sys.getsizeof(all_strand) + sum(sys.getsizeof(e) for e in all_strand)) / (1024 ** 3)}")


def _compute_target_size(guide: Guide, pam: PAM, offset: int) -> int:
    return len(guide) + len(pam) + offset


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
    # compute off-target size for extraction
    size = _compute_target_size(guide, pam, offset)  # offset is max(bdna, brna)
    loggers.verboselog.debug(f"Off-targets extraction size: {size}")
    # extract targets from reference genome fasta files
    extract_targets(fastas, pam, size, right, threads, loggers)
