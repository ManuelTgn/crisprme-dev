""" """

from .crisprme2_error import Crisprme2EnrichmentError
from .crisprme2_argparse import Crisprme2SearchInputArgs
from .logger import CrisprmeLoggers
from .variant import VariantRecord
from .fasta import Fasta
from .sample import Sample
from .vcf import VCF

from typing import List, Dict, Tuple, Optional
from concurrent.futures import Future
from pysam import TabixFile

import concurrent.futures
import numpy as np




from .utils import DNA, IUPACTABLE
from .pam import PAM

from glob import glob
from time import time

import sys
import os


# def read_genome(fasta_fnames: List[str], loggers: CrisprmeLoggers) -> List[GenomeFasta]:
#     genome = []  # construct genome (list) using input fasta files
#     for fasta_fname in fasta_fnames:
#         loggers.verboselog.debug(f"Loading FASTA file: {fasta_fname}")
#         try:
#             genome.append(GenomeFasta(fasta_fname, loggers))
#         except (Crisprme2FastaError, Exception):
#             loggers.errorlog.log_exception(
#                 f"Failed genome construction while loading {fasta_fname}", os.EX_DATAERR
#             )
#             sys.exit(os.EX_DATAERR)
#     assert genome  # must not be empty
#     return genome


# def _match(
#     bitset1: List[Bitset],
#     bitset2: List[Bitset],
#     position: int,
#     loggers: CrisprmeLoggers,
# ) -> bool:
#     # bitwise matching operation for input bitsets
#     # assumes the two bitsets have the same length
#     try:
#         return all((ntbit & bitset2[i]).to_bool() for i, ntbit in enumerate(bitset1))
#     except (ValueError, IndexError, Exception):
#         loggers.errorlog.log_raise_exception(
#             f"PAM bitwise matching with offtarget candidatefailed at position {position}",
#             os.EX_DATAERR,
#             Crisprme2BitsetError,
#         )
#         sys.exit(os.EX_DATAERR)


# def _match2(pam, offtargetpam):
#     for i, nt in enumerate(pam):
#         if offtargetpam[i] not in IUPACTABLE[nt]:
#             return False
#     return True


# def filter_offtarget(offtarget_pam: List[str], pam_patterns: List[Set[str]]) -> bool:
#     assert len(offtarget_pam) == len(pam_patterns)  # they must match
#     for i, nt in enumerate(offtarget_pam):
#         if nt not in pam_patterns[i]:  # not valid pam, skip target
#             return False
#     return True  # valid pam


# def _compute_pam_patterns(
#     pam: PAM, loggers: CrisprmeLoggers
# ) -> Tuple[List[Set[str]], List[Set[str]]]:
#     try:  # patterns used to filter targets based on input pam
#         pam_patterns_fw = [set(IUPACTABLE[nt]) for nt in pam.pam]  # forward patterns
#         pam_patterns_rc = [set(IUPACTABLE[nt]) for nt in pam.rc]  # reverse patterns
#     except (KeyError, Exception):
#         loggers.errorlog.log_exception(
#             f"Failed computing matching patterns for PAM: {pam}", os.EX_DATAERR
#         )
#         sys.exit(os.EX_DATAERR)
#     return pam_patterns_fw, pam_patterns_rc


# def _retrieve_pam(offtarget: List[str], right: bool, pamlen: int) -> List[str]:
#     return offtarget[:pamlen] if right else offtarget[-pamlen:]


# def fetch_offtargets(
#     # sequence: Sequence,
#     sequence: List[str],
#     pam: PAM,
#     guidepamlen: int,
#     right: bool,
#     loggers: CrisprmeLoggers,
# ) -> Tuple[List[str], List[str]]:
#     offtargets_fw, offtargets_rc = (
#         [],
#         [],
#     )  # iterate over sequence to fetch offtargets (with padding)
#     start_index, stop_index = PADDING, len(sequence) - PADDING
#     # total = sequence.stop_index - guidepamlen + 1 - sequence.start_index  # TODO: remove
#     total = stop_index - guidepamlen + 1 - start_index  # TODO: remove
#     progress_interval = max(1, total // 100)
#     # compute matching patterns for pam
#     pam_patterns_fw, pam_patterns_rc = _compute_pam_patterns(pam, loggers)
#     # for i in range(sequence.start_index, sequence.stop_index - guidepamlen + 1):
#     for i in range(start_index, stop_index - guidepamlen + 1):
#         if i % progress_interval == 0:
#             print(f"Progress: {((i + 1) / total) * 100:.2f}%", end="\r")
#         candidate = sequence[i - PADDING : i + guidepamlen + PADDING]
#         # candidate = sequence[i: i + guidepamlen]
#         # recover pam sequence from offtarget on forward and reverse strands
#         candidate_pam_fw = _retrieve_pam(candidate, right, len(pam))  # type: ignore
#         candidate_pam_rc = _retrieve_pam(candidate, (not right), len(pam))  # type: ignore
#         if filter_offtarget(candidate_pam_fw, pam_patterns_fw):  # check on fw
#             # offtargets_fw.append(sequence.fetch(i, i + guidepamlen))
#             offtargets_fw.append(candidate)
#         if filter_offtarget(candidate_pam_rc, pam_patterns_rc):  # check on rev
#             offtargets_fw.append(candidate)
#             # offtargets_rc.append(sequence.fetch(i, i + guidepamlen))
#     print()
#     return offtargets_fw, offtargets_rc


# def compute_offtargets(
#     genome: List[GenomeFasta],
#     pam: PAM,
#     guidelen: int,
#     right: bool,
#     outdir: str,
#     loggers: CrisprmeLoggers,
# ):
#     guidepamlen = len(pam) + guidelen  # compute guide + pam length
#     for contig in genome:  # iterate over each genome contig
#         loggers.verboselog.debug(
#             f"Fetching off-target candidates from contig: {contig.contig}"
#         )
#         start = time()
#         contig_seq = contig.read()  # read contig sequence
#         offtargets = fetch_offtargets(contig_seq, pam, guidepamlen, right, loggers)
#         loggers.verboselog.debug(
#             f"Fetched {len(offtargets[0])} on 5'-3' and {len(offtargets[1])} on 3'-5' on contig {contig.contig}"
#         )
#         loggers.verboselog.debug(
#             f"Off-target candidates fetched from contig {contig.contig} in {time() - start:.2f}s"
#         )
#         # # TODO: after check remove
#         # with open(os.path.join(outdir, f"offtargets_fw_{contig.contig}.txt"), mode="w") as outfile:
#         #     outfile.write("\n".join([ot[PADDING:PADDING+guidepamlen]for ot in offtargets[0]]))
#         # with open(os.path.join(outdir, f"offtargets_rc_{contig.contig}.txt"), mode="w") as outfile:
#         #     outfile.write("\n".join(ot[PADDING:PADDING+guidelen] for ot in offtargets[1]))


# # def process_genome(
# #     fasta_fnames: List[str],
# #     pam: PAM,
# #     guidelen: int,
# #     right: bool,
# #     outdir: str,
# #     loggers: CrisprmeLoggers,
# # ):
# #     loggers.basiclog.info(
# #         f"Reconstructing alternative genomes and retrieving off-targets"
# #     )
# #     genome = read_genome(fasta_fnames, loggers)  # load input genome data
# #     # assumes input guides share the same length
# #     compute_offtargets(genome, pam, guidelen, right, outdir, loggers)



def read_fasta_files(fasta_fnames: List[str], loggers: CrisprmeLoggers) -> Dict[str, Fasta]:
    fasta_map = {}  # contig - fasta map
    for fasta_fname in fasta_fnames:
        loggers.verboselog.debug(f"Create Fasta object: {fasta_fname}")
        with Fasta(fasta_fname,loggers) as fasta:
            if fasta.nreferences != 1:  # fasta files must be chromosome-separated
                loggers.errorlog.log_raise_exception(f"Fasta file {fasta_fname} contains multiple contig data", os.EX_DATAERR, Crisprme2EnrichmentError)
            contig = fasta.contig # assume one single contig
            if contig in fasta_map:
                loggers.errorlog.log_raise_exception(f"Multiple Fasta files with contig {contig}", os.EX_DATAERR, Crisprme2EnrichmentError)
            fasta_map[contig] = fasta
        loggers.verboselog.debug(f"Successfully created Fasta object: {fasta_fname}")
    return fasta_map


def read_vcf_files(vcf_fnames: List[str], loggers: CrisprmeLoggers) -> Dict[str, VCF]:
    vcf_map = {}
    for vcf_fname in vcf_fnames:
        loggers.verboselog.debug(f"Create VCF object: {vcf_fname}")
        vcf = VCF(vcf_fname, loggers)  # create vcf object
        contig = vcf.contig  # assume one single contig
        if contig in vcf_map:
            loggers.errorlog.log_raise_exception(f"Multiple VCF files with contig {contig}", os.EX_DATAERR, Crisprme2EnrichmentError)
        vcf_map[contig] = vcf
        loggers.verboselog.debug(f"Successfully created VCF object: {vcf_fname}")
    return vcf_map



def create_fasta_vcf_map(fasta_fnames: List[str], vcf_fnames: List[str], loggers: CrisprmeLoggers) -> Dict[str, Tuple[Fasta, Optional[VCF]]]:
    fasta_map = read_fasta_files(fasta_fnames, loggers)
    vcf_map = read_vcf_files(vcf_fnames, loggers)
    return {contig: (fasta_map[contig], vcf_map[contig]) if contig in vcf_map else (fasta_map[contig], None) for contig in fasta_map}
    

def construct_samples_list(fasta_vcf_map: Dict[str, Tuple[Fasta, Optional[VCF]]], loggers: CrisprmeLoggers) -> List[Sample]:
    vcf = None
    for _, (_, vcf_) in fasta_vcf_map.items():
        if vcf_ is not None:  # take first available vcf in dataset
            vcf = vcf_
            break
    samples = [Sample("REF", loggers)]  # create reference "sample"
    if vcf is not None:  # input vcf data
        # assumption: all input VCFs share the same samples set
        samples += [Sample(sample, loggers) for sample in vcf.get_samples()]
    return samples


def _split_ranges(length: int, threads: int, loggers: CrisprmeLoggers, overlap: int = 100) -> List[Tuple[int, int]]:
    if length <= 0 or threads <= 0:  # should never happen
        loggers.errorlog.log_raise_exception(f"Empty string or 0 threads used", os.EX_DATAERR, Crisprme2EnrichmentError)
    # compute the base size of the non-overlapping portion of each chunk
    # np.ceil ensures that the sum of the base sizes is at least the total string 
    # length, guaranteeing full coverage
    base_chunk_size = int(np.ceil(length / threads))
    ranges: List[Tuple[int, int]] = []
    for i in range(threads):
        # calculate the non-overlapping start index for the current chunk
        non_overlapping_start = i * base_chunk_size
        # calculate the non-overlapping end index (exclusive)
        non_overlapping_end = min((i + 1) * base_chunk_size, length)
        # apply overlap to start and end indexes
        start_index = max(0, non_overlapping_start - overlap)
        end_index = min(length, non_overlapping_end + overlap)
        # if the end index calculation results in a slice that is smaller than 
        # the required overlap ensure the end index is at least the string length
        # if it's the last chunk
        if i == threads - 1:
            end_index = length
        if start_index < end_index:  # skip if start and end are the same
            ranges.append((start_index, end_index))
    return ranges

def _collect_tasks(vcf: VCF, contig_length: int, threads: int, loggers: CrisprmeLoggers):
    return [
        (vcf, vcf.contig, start, stop, loggers)
        for start, stop in _split_ranges(contig_length, threads, loggers, overlap=0)
    ]  # collect tasks to perform in parallel (no overlap required)

def _retrieve_variants_range(vcf: VCF, contig: Optional[str], start: Optional[int], stop: Optional[int], samples: List[Sample], loggers: CrisprmeLoggers) -> int:
    try:
        with TabixFile(vcf.filepath, index=vcf.index) as tbx:
            # variants = [1 for v in tbx.fetch(contig, start, stop)]
            variants = [VariantRecord(v, samples, vcf.phased, vcf.ploidy, loggers) for v in tbx.fetch(contig, start, stop)]
            loggers.verboselog.debug(f"Retrieved {len(variants)} variants in {contig}:{start}-{stop}")
            return len(variants)
    except Exception as e:
        loggers.errorlog.log_exception(f"Error retrieving variants in {contig}:{start}-{stop}: {e}", os.EX_DATAERR)

def _collect_results(futures: List[Future], loggers: CrisprmeLoggers) -> int:
    results_all = []
    for future in concurrent.futures.as_completed(futures):
        try:
            results_all.append(future.result())
        except Exception as e:
            loggers.errorlog.log_exception(f"Error in parallel variant retrieval task: {e}", os.EX_DATAERR)
    return sum(results_all)

def retrieve_variants_contig(vcf: VCF, contig_length: int, samples: List[Sample], threads: int, loggers: CrisprmeLoggers) -> int:
    tasks = _collect_tasks(vcf, contig_length, threads, loggers)  # collect tasks
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(_retrieve_variants_range, vcf, contig, start, stop, samples, loggers)
            for vcf, contig, start, stop, loggers in tasks
        ]
        return _collect_results(futures, loggers)


def reconstruct_targets(fasta_vcf_map: Dict[str, Tuple[Fasta, Optional[VCF]]], samples: List[Sample], threads: int, loggers: CrisprmeLoggers):
    tasks = []  # tasks collector item
    for contig, (fasta, vcf) in fasta_vcf_map.items():
        with fasta as f:
            if vcf is not None:
                print(contig)
                variants = vcf.read(threads=threads)
                print(len(variants))
                




def enrich_genome(args: Crisprme2SearchInputArgs, loggers: CrisprmeLoggers):
    fasta_vcf_map = create_fasta_vcf_map(args.fastas, args.vcfs, loggers)
    samples = construct_samples_list(fasta_vcf_map, loggers)
    reconstruct_targets(fasta_vcf_map, samples, args.threads, loggers)
    
    