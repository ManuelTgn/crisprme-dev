""" """

from .crisprme2_error import Crisprme2FastaError, Crisprme2BitsetError
from .sequence import PADDING, GenomeFasta, Sequence
from .utils import DNA, IUPACTABLE
from .logger import CrisprmeLoggers
from .encoder import encode
from .bitset import Bitset
from .pam import PAM

from typing import List, Tuple, Set
from time import time

import sys
import os


def read_genome(fasta_fnames: List[str], loggers: CrisprmeLoggers) -> List[GenomeFasta]:
    genome = []  # construct genome (list) using input fasta files
    for fasta_fname in fasta_fnames:
        loggers.verboselog.debug(f"Loading FASTA file: {fasta_fname}")
        try:
            genome.append(GenomeFasta(fasta_fname, loggers))
        except (Crisprme2FastaError, Exception):
            loggers.errorlog.log_exception(
                f"Failed genome construction while loading {fasta_fname}", os.EX_DATAERR
            )
            sys.exit(os.EX_DATAERR)
    assert genome  # must not be empty
    return genome


def _match(
    bitset1: List[Bitset],
    bitset2: List[Bitset],
    position: int,
    loggers: CrisprmeLoggers,
) -> bool:
    # bitwise matching operation for input bitsets
    # assumes the two bitsets have the same length
    try:
        return all((ntbit & bitset2[i]).to_bool() for i, ntbit in enumerate(bitset1))
    except (ValueError, IndexError, Exception):
        loggers.errorlog.log_raise_exception(
            f"PAM bitwise matching with offtarget candidatefailed at position {position}",
            os.EX_DATAERR,
            Crisprme2BitsetError,
        )
        sys.exit(os.EX_DATAERR)


def _match2(pam, offtargetpam):
    for i, nt in enumerate(pam):
        if offtargetpam[i] not in IUPACTABLE[nt]:
            return False
    return True


def filter_offtarget(offtarget_pam: List[str], pam_patterns: List[Set[str]]) -> bool:
    assert len(offtarget_pam) == len(pam_patterns)  # they must match
    for i, nt in enumerate(offtarget_pam):
        if nt not in pam_patterns[i]:  # not valid pam, skip target
            return False
    return True  # valid pam


def _compute_pam_patterns(
    pam: PAM, loggers: CrisprmeLoggers
) -> Tuple[List[Set[str]], List[Set[str]]]:
    try:  # patterns used to filter targets based on input pam
        pam_patterns_fw = [set(IUPACTABLE[nt]) for nt in pam.pam]  # forward patterns
        pam_patterns_rc = [set(IUPACTABLE[nt]) for nt in pam.rc]  # reverse patterns
    except (KeyError, Exception):
        loggers.errorlog.log_exception(
            f"Failed computing matching patterns for PAM: {pam}", os.EX_DATAERR
        )
        sys.exit(os.EX_DATAERR)
    return pam_patterns_fw, pam_patterns_rc


def _retrieve_pam(offtarget: List[str], right: bool, pamlen: int) -> List[str]:
    return offtarget[:pamlen] if right else offtarget[-pamlen:]


def fetch_offtargets(
    sequence: Sequence,
    pam: PAM,
    guidepamlen: int,
    right: bool,
    loggers: CrisprmeLoggers,
) -> Tuple[List[str], List[str]]:
    offtargets_fw, offtargets_rc = (
        [],
        [],
    )  # iterate over sequence to fetch offtargets (with padding)
    total = sequence.stop_index - guidepamlen + 1 - sequence.start_index  # TODO: remove
    progress_interval = max(1, total // 100)
    # compute matching patterns for pam
    pam_patterns_fw, pam_patterns_rc = _compute_pam_patterns(pam, loggers)
    for i in range(sequence.start_index, sequence.stop_index - guidepamlen + 1):
        if i % progress_interval == 0:
            print(f"Progress: {((i + 1) / total) * 100:.2f}%", end="\r")
        candidate = sequence[i : i + guidepamlen]
        # recover pam sequence from offtarget on forward and reverse strands
        candidate_pam_fw = _retrieve_pam(candidate, right, len(pam))  # type: ignore
        candidate_pam_rc = _retrieve_pam(candidate, (not right), len(pam))  # type: ignore
        if filter_offtarget(candidate_pam_fw, pam_patterns_fw):  # check on fw
            offtargets_fw.append(sequence.fetch(i, i + guidepamlen))
        if filter_offtarget(candidate_pam_rc, pam_patterns_rc):  # check on rev
            offtargets_rc.append(sequence.fetch(i, i + guidepamlen))
    print()
    return offtargets_fw, offtargets_rc


def compute_offtargets(
    genome: List[GenomeFasta],
    pam: PAM,
    guidelen: int,
    right: bool,
    outdir: str,
    loggers: CrisprmeLoggers,
):
    guidepamlen = len(pam) + guidelen  # compute guide + pam length
    for contig in genome:  # iterate over each genome contig
        loggers.verboselog.debug(
            f"Fetching off-target candidates from contig: {contig.contig}"
        )
        start = time()
        contig.read()  # read contig sequence
        offtargets = fetch_offtargets(contig.sequence, pam, guidepamlen, right, loggers)
        loggers.verboselog.debug(
            f"Fetched {len(offtargets[0])} on 5'-3' and {len(offtargets[1])} on 3'-5' on contig {contig.contig}"
        )
        loggers.verboselog.debug(
            f"Off-target candidates fetched from contig {contig.contig} in {time() - start:.2f}s"
        )
        # # TODO: after check remove
        # with open(os.path.join(outdir, f"offtargets_fw_{contig.contig}.txt"), mode="w") as outfile:
        #     outfile.write("\n".join([ot[PADDING:PADDING+guidepamlen]for ot in offtargets[0]]))
        # with open(os.path.join(outdir, f"offtargets_rc_{contig.contig}.txt"), mode="w") as outfile:
        #     outfile.write("\n".join(ot[PADDING:PADDING+guidelen]for ot in offtargets[1]))


def process_genome(
    fasta_fnames: List[str],
    pam: PAM,
    guidelen: int,
    right: bool,
    outdir: str,
    loggers: CrisprmeLoggers,
):
    loggers.basiclog.info(
        f"Reconstructing alternative genomes and retrieving off-targets"
    )
    genome = read_genome(fasta_fnames, loggers)  # load input genome data
    # assumes input guides share the same length
    compute_offtargets(genome, pam, guidelen, right, outdir, loggers)
