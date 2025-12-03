""" """

from .crisprme2_error import Crisprme2EnrichmentError
from .crisprme2_argparse import Crisprme2SearchInputArgs
from .logger import CrisprmeLoggers
from .sequence import ContigSequence
from .fasta import Fasta
from .sample import Sample
from .vcf import VCF
from .pam import PAM

from .target_candidates_parser import (
    find_target_candidates
)  # defined in rust/src/lib.rs

from typing import List, Dict, Tuple, Optional
from time import time

import os

# contig chunk size and overlap
CHUNKSIZE = 10_000_000
CHUNKOVERLAP = 500


def read_fasta_files(
    fasta_fnames: List[str], loggers: CrisprmeLoggers
) -> Dict[str, Fasta]:
    fasta_map = {}  # contig - fasta map
    for fasta_fname in fasta_fnames:
        loggers.verboselog.debug(f"Create Fasta object: {fasta_fname}")
        with Fasta(fasta_fname, loggers) as fasta:
            if fasta.nreferences != 1:  # fasta files must be chromosome-separated
                loggers.errorlog.log_raise_exception(
                    f"Fasta file {fasta_fname} contains multiple contig data",
                    os.EX_DATAERR,
                    Crisprme2EnrichmentError,
                )
            contig = fasta.contig  # assume one single contig
            if contig in fasta_map:
                loggers.errorlog.log_raise_exception(
                    f"Multiple Fasta files with contig {contig}",
                    os.EX_DATAERR,
                    Crisprme2EnrichmentError,
                )
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
            loggers.errorlog.log_raise_exception(
                f"Multiple VCF files with contig {contig}",
                os.EX_DATAERR,
                Crisprme2EnrichmentError,
            )
        vcf_map[contig] = vcf
        loggers.verboselog.debug(f"Successfully created VCF object: {vcf_fname}")
    return vcf_map


def create_fasta_vcf_map(
    fasta_fnames: List[str], vcf_fnames: List[str], loggers: CrisprmeLoggers
) -> Dict[str, Tuple[Fasta, Optional[VCF]]]:
    fasta_map = read_fasta_files(fasta_fnames, loggers)
    vcf_map = read_vcf_files(vcf_fnames, loggers)
    return {
        contig: (
            (fasta_map[contig], vcf_map[contig])
            if contig in vcf_map
            else (fasta_map[contig], None)
        )
        for contig in fasta_map
    }


def construct_samples_list(
    fasta_vcf_map: Dict[str, Tuple[Fasta, Optional[VCF]]], loggers: CrisprmeLoggers
) -> List[Sample]:
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


# def _split_ranges(length: int, threads: int, loggers: CrisprmeLoggers, overlap: int = 100) -> List[Tuple[int, int]]:
#     if length <= 0 or threads <= 0:  # should never happen
#         loggers.errorlog.log_raise_exception(f"Empty string or 0 threads used", os.EX_DATAERR, Crisprme2EnrichmentError)
#     # compute the base size of the non-overlapping portion of each chunk
#     # np.ceil ensures that the sum of the base sizes is at least the total string
#     # length, guaranteeing full coverage
#     base_chunk_size = int(np.ceil(length / threads))
#     ranges: List[Tuple[int, int]] = []
#     for i in range(threads):
#         # calculate the non-overlapping start index for the current chunk
#         non_overlapping_start = i * base_chunk_size
#         # calculate the non-overlapping end index (exclusive)
#         non_overlapping_end = min((i + 1) * base_chunk_size, length)
#         # apply overlap to start and end indexes
#         start_index = max(0, non_overlapping_start - overlap)
#         end_index = min(length, non_overlapping_end + overlap)
#         # if the end index calculation results in a slice that is smaller than
#         # the required overlap ensure the end index is at least the string length
#         # if it's the last chunk
#         if i == threads - 1:
#             end_index = length
#         if start_index < end_index:  # skip if start and end are the same
#             ranges.append((start_index, end_index))
#     return ranges

# def _collect_tasks(vcf: VCF, contig_length: int, threads: int, loggers: CrisprmeLoggers):
#     return [
#         (vcf, vcf.contig, start, stop, loggers)
#         for start, stop in _split_ranges(contig_length, threads, loggers, overlap=0)
#     ]  # collect tasks to perform in parallel (no overlap required)

# def _retrieve_variants_range(vcf: VCF, contig: Optional[str], start: Optional[int], stop: Optional[int], samples: List[Sample], loggers: CrisprmeLoggers) -> int:
#     try:
#         with TabixFile(vcf.filepath, index=vcf.index) as tbx:
#             # variants = [1 for v in tbx.fetch(contig, start, stop)]
#             variants = [VariantRecord(v, samples, vcf.phased, vcf.ploidy, loggers) for v in tbx.fetch(contig, start, stop)]
#             loggers.verboselog.debug(f"Retrieved {len(variants)} variants in {contig}:{start}-{stop}")
#             return len(variants)
#     except Exception as e:
#         loggers.errorlog.log_exception(f"Error retrieving variants in {contig}:{start}-{stop}: {e}", os.EX_DATAERR)

# def _collect_results(futures: List[Future], loggers: CrisprmeLoggers) -> int:
#     results_all = []
#     for future in concurrent.futures.as_completed(futures):
#         try:
#             results_all.append(future.result())
#         except Exception as e:
#             loggers.errorlog.log_exception(f"Error in parallel variant retrieval task: {e}", os.EX_DATAERR)
#     return sum(results_all)

# def retrieve_variants_contig(vcf: VCF, contig_length: int, samples: List[Sample], threads: int, loggers: CrisprmeLoggers) -> int:
#     tasks = _collect_tasks(vcf, contig_length, threads, loggers)  # collect tasks
#     with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
#         futures = [
#             executor.submit(_retrieve_variants_range, vcf, contig, start, stop, samples, loggers)
#             for vcf, contig, start, stop, loggers in tasks
#         ]
#         return _collect_results(futures, loggers)


# def chunk_contig_sequence(contig_sequence: ContigSequence, threads: int):
#     contig_chunks = [c for c in contig_sequence.chunk(CHUNKSIZE, 0)]
#     targets = []
#     tot = len(contig_chunks)
#     start_time = time()
#     for i, contig_chunk in enumerate(contig_chunks):

#         targets.extend(extract_targets_parallel(contig_chunk.sequence, 23, threads))
#         print(f"Progress: {(((i + 1) / tot) * 100):.2f}%%", end="\r")
#     print()
#     print(f"targets: {len(targets)}")
#     print(f"elapsed time {time() - start_time:.2f}s")
#     print()


# def reconstruct_targets(
#     fasta_vcf_map: Dict[str, Tuple[Fasta, Optional[VCF]]],
#     samples: List[Sample],
#     threads: int,
#     loggers: CrisprmeLoggers,
# ):
#     tasks = []  # tasks collector item
#     for contig, (fasta, vcf) in fasta_vcf_map.items():
#         with fasta as f:
#             print(contig, f.length)

#             chunk_contig_sequence(f.fetch(contig), threads)
#             # if vcf is not None:
#             #     print(contig)
#             #     start_time1 = time()
#             #     # t1 = threads // 2
#             #     # t2 = threads - t1
#             #     ranges = _split_ranges(fasta.length, threads, loggers, 500)
#             #     for start, stop in ranges:
#             #         start_time2 = time()
#             #         variants = vcf.read(start=start, stop=stop, threads=1)
#             #         # reader =  cyvcf2.VCF(vcf.filepath, mode="r", threads=1, lazy=True)
#             #         # variants = [v for v in reader]
#             #         # del reader
#             #         print(len(variants), f"region: {contig}:{start}-{stop}\ttime: {time() - start_time2:.2f}s")
#             #     print(f"contig: {contig}\ttotal time: {time() - start_time1:.2f}s")


def _chunk_contig_sequence(contig_sequence: ContigSequence) -> List[ContigSequence]:
    # chunk each contig in 10Mb chunks
    return [c for c in contig_sequence.chunk(CHUNKSIZE, CHUNKOVERLAP)]

def _find_target_candidates(contig_sequence: str, contig: str, pam_seq: str, offset: int, right: bool, is_first_iteration: bool, index_path: str, threads: int):
    find_target_candidates(contig_sequence, contig, pam_seq, offset, right, is_first_iteration, index_path, threads)

def _scan_sequence(contig_seq: ContigSequence, contig: str, pam_seq: str, offset: int, outdir: str, right: bool, threads: int, loggers: CrisprmeLoggers) -> None:
    # split contig sequence in 10 Mb long chunks
    contig_chunks = _chunk_contig_sequence(contig_seq)
    index_path = os.path.join(outdir, f"{contig}")  # define targets index path
    is_first_iteration = True  # first iteration (create indexes)
    for chunk in contig_chunks:  # scan sequence to extract targets
        _find_target_candidates(chunk.sequence.upper(), contig, pam_seq, offset, right, is_first_iteration, index_path, threads)
        is_first_iteration = False  # append to index files


def retrieve_targets(fasta_vcf_map: Dict[str, Tuple[Fasta, Optional[VCF]]], pam: PAM, guidelen: int, offset: int, right: bool, threads: int, outdir: str, loggers: CrisprmeLoggers):
    # use offset to account for bulges in alignments
    guidelen_offset = guidelen + len(pam) + offset
    for contig, (fasta, vcf) in fasta_vcf_map.items():
        with fasta as f:
            start = time()
            loggers.verboselog.debug(f"Scanning contig {contig} for targets (use {threads} threads)")
            _scan_sequence(f.fetch(contig), contig, pam.pam, guidelen_offset, outdir, right, threads, loggers)
            loggers.verboselog.debug(f"Scanning contig {contig} completed in {time() - start:.2f}s")

def retrieve_target_candidates(args: Crisprme2SearchInputArgs, pam: PAM, guidelen: int, offset: int, loggers: CrisprmeLoggers):
    # map each contig fasta to its variant data
    fasta_vcf_map = create_fasta_vcf_map(args.fastas, args.vcfs, loggers)
    retrieve_targets(fasta_vcf_map, pam, guidelen, offset, args.right, args.threads, args.outdir, loggers)
