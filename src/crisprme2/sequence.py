""" """

from .crisprme2_error import Crisprme2SequenceError, Crisprme2ContigSequenceError
from .logger import CrisprmeLoggers
from .utils import RC

from typing import Optional, Union, List, Generator

import sys
import os


class SequenceStats:

    def __init__(self, length: int, n_count: int) -> None:
        self._length = length
        self._n_count = n_count

    @property
    def length(self) -> int:
        return self._length
    
    @property
    def n_count(self) -> int:
        return self._n_count


class Sequence:
    
    def __init__(self, sequence: str, loggers: CrisprmeLoggers):
        self._loggers = loggers  # store loggers
        self._sequence = list(sequence)  # store sequence as list of str
        self._length = len(self._sequence)
        self._stats: Optional[SequenceStats] = None  # basic sequence stats

    def __len__(self) -> int:
        return self._length
    
    def __str__(self) -> str:
        return "".join(self._sequence)
    
    def subsequence(self, start: int, end: int) -> str:
        if start < 0 or end > len(self) or start >= end:
            self._loggers.errorlog.log_raise_exception(f"Invalid coordinates: start={start}, end={end}, length={len(self)}", os.EX_DATAERR, Crisprme2SequenceError)
        return "".join(self._sequence[start:end])
    
    def reverse_complement(self) -> List[str]:
        try:
            return [RC[nt] for nt in self._sequence[::-1]]
        except (KeyError, Exception) as e:
            self._loggers.errorlog.log_raise_exception(f"Error computing reverse complement: {str(e)}", os.EX_DATAERR, Crisprme2SequenceError)
    
    def calculate_statistics(self) -> SequenceStats:
        if self._stats is not None:
            return self._stats
        n_count = sum(1 for nt in self._sequence if nt.upper() == "N")    
        self._stats = SequenceStats(length=len(self), n_count=n_count)
        return self._stats
        
    @property
    def sequence(self) -> str:
        return "".join(self._sequence)
    

class ContigSequence(Sequence):
    
    def __init__(self, sequence: str, contig: str, start: int, stop: int, loggers: CrisprmeLoggers) -> None:
        super().__init__(sequence, loggers)
        self._contig = contig  # store sequence contig name
        self._start = start  # start sequence start position
        self._stop = stop  # store sequence stop position

    def chunk(self, size: int, overlap: int):
        if overlap >= size:
            self._loggers.errorlog.log_raise_exception(f"Overlap size ({overlap}) must be less than the chunk size ({size})", os.EX_DATAERR, Crisprme2ContigSequenceError)
        if size >= self._length:  # handle short contigs
            size = self._length
        step = size - overlap  # e.g., size=100 and overlap=10, step is 90
        for i in range(0, self._length, step):
            start, stop = i, i + size
            if stop >= self._length:
                stop = self._length
            yield ContigSequence(self.subsequence(start, stop), self._contig, start, stop, self._loggers)


    @property
    def contig(self) -> str:
        return self._contig
    
    @property
    def start(self) -> int:
        return self._start
    
    @property
    def stop(self) -> int:
        return self._stop
    







# from .crisprme2_error import Crisprme2FastaError
# from .logger import CrisprmeLoggers

# from typing import Union, List
# from pysam.utils import SamtoolsError
# from time import time

# import pysam
# import sys
# import os

# FAI = "fai"  # fasta index extension format
# # TODO: consider making it an input param
# PADDING = 100  # up/dowstream padding length

# # fasta extensions
# FASTAEXTENSIONS = {"fa", "fasta"}


# class Sequence:

#     def __init__(self, sequence: str, loggers: CrisprmeLoggers) -> None:
#         self._loggers = loggers  # store loggers
#         self._sequence = list(sequence)  # sequence as list
#         self._hash = None  # object hash (pre-computed for efficiency)
#         # define sequence start and stop boundaries
#         self._start_idx = PADDING  # assumes sequences longer than PADDING
#         self._stop_idx = len(sequence) - PADDING

#     def __repr__(self) -> str:
#         return f"<{self.__class__.__name__} object; sequence={self.sequence}>"

#     def __str__(self) -> str:
#         return "".join(self._sequence)

#     def __eq__(self, value: object) -> bool:
#         if not isinstance(value, Sequence):
#             return NotImplemented
#         return "".join(self._sequence) == value.sequence

#     def __hash__(self) -> int:
#         if self._hash is None:
#             self._hash = hash("".join(self._sequence))
#         return self._hash

#     def __len__(self) -> int:
#         return len(self._sequence)

#     def __getitem__(self, idx: Union[int, slice]) -> Union[str, List[str]]:
#         if not hasattr(self, "_sequence"):
#             self._loggers.errorlog.log_raise_exception(
#                 f"Missing _sequence attribute on class {self.__class__.__name__}",
#                 os.EX_DATAERR,
#                 AttributeError,
#             )
#         try:
#             return self._sequence[idx]
#         except IndexError:
#             self._loggers.errorlog.log_exception(
#                 f"Index {idx} out of bounds", os.EX_DATAERR
#             )
#             sys.exit(os.EX_DATAERR)

#     def __iter__(self) -> "SequenceIterator":
#         return SequenceIterator(self, self._loggers)

#     def fetch(self, start: int, stop: int) -> List[str]:
#         if start < PADDING:
#             self._loggers.errorlog.log_raise_exception(
#                 f"start index ({start}) out of range", os.EX_DATAERR, ValueError
#             )
#         if stop > len(self) - PADDING:
#             self._loggers.errorlog.log_raise_exception(
#                 f"stop index ({stop}) out of range", os.EX_DATAERR, ValueError
#             )
#         return self._sequence[start - PADDING : stop + PADDING]

#     @property
#     def sequence(self) -> str:
#         return "".join(self._sequence)

#     @property
#     def start_index(self) -> int:
#         return self._start_idx

#     @property
#     def stop_index(self) -> int:
#         return self._stop_idx


# class SequenceIterator:

#     def __init__(self, sequence: Sequence, loggers: CrisprmeLoggers) -> None:
#         self._loggers = loggers
#         if not hasattr(sequence, "_sequence_"):
#             self._loggers.errorlog.log_raise_exception(
#                 f"Missing _sequence_ attribute on class {self.__class__.__name__}",
#                 os.EX_DATAERR,
#                 AttributeError,
#             )
#         self._sequence = sequence  # sequence object to iterate
#         self._index = 0  # iterator index

#     def __next__(self) -> str:
#         if self._index < len(self._sequence):
#             result = self._sequence[self._index]
#             assert isinstance(result, str)  # slices not allowed here
#             self._index += 1  # go to next position in sequence
#             return result
#         raise StopIteration  # stop iteration over sequence object


# class Fasta:

#     def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
#         self._loggers = loggers  # store loggers
#         self._fname = fname  # store input file
#         self._faidx = self._index_fasta()  # initialize fasta index

#     def _index_fasta(self) -> str:
#         # look for index file for the current fasta file, if not found compute it
#         if _find_fai(self._fname):  # index found, return it
#             return f"{self._fname}.{FAI}"
#         self._loggers.verboselog.debug(
#             f"FASTA index not found for {self._fname}. Generating FASTA index"
#         )
#         start = time()  # measure indexing  time
#         try:
#             pysam.faidx(self._fname)  # index fasta using samtools
#         except (SamtoolsError, Exception):
#             self._loggers.errorlog.log_exception(
#                 f"An error occurred while indexing {self._fname}", os.EX_SOFTWARE
#             )
#         assert _find_fai(self._fname)  # index must be available now
#         self._loggers.verboselog.debug(
#             f"FASTA index for {self._fname} computed in {time() - start: .2f}s"
#         )
#         return f"{self._fname}.{FAI}"


# def _find_fai(fastafile: str) -> bool:
#     """Check if a FASTA index file exists for the given FASTA file.

#     Checks if a FASTA index file (.fai) exists for the given FASTA file in the same
#     directory.

#     Args:
#         fastafile: The path to the FASTA file.

#     Returns:
#         True if the index file exists and is not empty, False otherwise.
#     """
#     # search index for the input fasta file, assumes that the index is located
#     # in the same folder as the indexed fasta
#     fastaindex = f"{os.path.abspath(fastafile)}.{FAI}"  # avoid unexpected crashes
#     if os.path.exists(fastaindex):
#         return os.path.isfile(fastaindex) and os.stat(fastaindex).st_size > 0
#     return False


# class GenomeFasta(Fasta):

#     def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
#         super().__init__(fname, loggers)
#         self._contig = self._retrieve_contig()  # initialize contig name
#         self._length = self._compute_contig_length()  # initialize fasta length

#     def __len__(self) -> int:
#         return self._length

#     def _retrieve_contig(self) -> str:
#         f = pysam.FastaFile(self._fname, filepath_index=self._faidx)
#         if len(f.references) != 1:  # fastas are expected to be chromosome-wise
#             self._loggers.errorlog.log_raise_exception(
#                 f"Unexpected number of contigs ({len(f.references)}) found in {self._fname}",
#                 os.EX_DATAERR,
#                 Crisprme2FastaError,
#             )
#         return f.references[0]  # contig name

#     def _compute_contig_length(self) -> int:
#         f = pysam.FastaFile(self._fname, filepath_index=self._faidx)
#         if len(f.lengths) != 1:  # fastas are expected to be chromosome-wise
#             self._loggers.errorlog.log_raise_exception(
#                 f"Unexpected number of contigs ({len(f.references)}) found in {self._fname}",
#                 os.EX_DATAERR,
#                 Crisprme2FastaError,
#             )
#         return f.lengths[0]  # contig name

#     def read(self) -> List[str]:
#         try:  # read fasta file sequence content
#             with open(self._fname, mode="r") as infile:
#                 infile.readline()  # skip fasta header
#                 # read fasta content
#                 # self._sequence = Sequence(
#                 #     "".join([line.strip() for line in infile.readlines()]),
#                 #     self._loggers,
#                 # )
#                 return list("".join([line.strip() for line in infile.readlines()]))
#         except (IOError, Exception):
#             self._loggers.errorlog.log_exception(
#                 f"An error occurred while reading {self._fname}", os.EX_IOERR
#             )
#             sys.exit(os.EX_IOERR)

#     # @property
#     # def sequence(self) -> Sequence:
#     #     return self._sequence

#     @property
#     def contig(self) -> str:
#         return self._contig


# class GuideFasta(Fasta):

#     def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
#         super().__init__(fname, loggers)
#         self._read_guides()  # read guides in fasta

#     def _read_guides(self) -> None:
#         f = pysam.FastaFile(self._fname, filepath_index=self._faidx)
#         gnames = f.references  # retrieve guides seqnames (fasta headers)
#         try:
#             self._guides = list(
#                 {Sequence(f.fetch(gname), self._loggers) for gname in gnames}
#             )  # extract guide sequence
#         except (SamtoolsError, Exception):
#             self._loggers.errorlog.log_exception(
#                 f"Failed parsing guides from {self._fname}", os.EX_DATAERR
#             )

#     @property
#     def guides(self) -> List[Sequence]:
#         return self._guides

