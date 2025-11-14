""" """

from .crisprme2_error import Crisprme2FastaError
from .logger import CrisprmeLoggers

from typing import Union, List
from pysam.utils import SamtoolsError
from time import time

import pysam
import sys
import os

FAI = "fai"  # fasta index extension format
# TODO: consider making it an input param
PADDING = 100  # up/dowstream padding length


class Sequence:

    def __init__(self, sequence: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._sequence = list(sequence)  # sequence as list
        self._hash = None  # object hash (pre-computed for efficiency)
        # define sequence start and stop boundaries
        self._start_idx = PADDING  # assumes sequences longer than PADDING
        self._stop_idx = len(sequence) - PADDING

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; sequence={self.sequence}>"

    def __str__(self) -> str:
        return "".join(self._sequence)

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, Sequence):
            return NotImplemented
        return "".join(self._sequence) == value.sequence

    def __hash__(self) -> int:
        if self._hash is None:
            self._hash = hash("".join(self._sequence))
        return self._hash

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, idx: Union[int, slice]) -> Union[str, List[str]]:
        if not hasattr(self, "_sequence"):
            self._loggers.errorlog.log_raise_exception(
                f"Missing _sequence attribute on class {self.__class__.__name__}",
                os.EX_DATAERR,
                AttributeError,
            )
        try:
            return self._sequence[idx]
        except IndexError:
            self._loggers.errorlog.log_exception(
                f"Index {idx} out of bounds", os.EX_DATAERR
            )
            sys.exit(os.EX_DATAERR)

    def __iter__(self) -> "SequenceIterator":
        return SequenceIterator(self, self._loggers)

    def fetch(self, start: int, stop: int) -> List[str]:
        if start < PADDING:
            self._loggers.errorlog.log_raise_exception(
                f"start index ({start}) out of range", os.EX_DATAERR, ValueError
            )
        if stop > len(self) - PADDING:
            self._loggers.errorlog.log_raise_exception(
                f"stop index ({stop}) out of range", os.EX_DATAERR, ValueError
            )
        return self._sequence[start - PADDING : stop + PADDING]

    @property
    def sequence(self) -> str:
        return "".join(self._sequence)

    @property
    def start_index(self) -> int:
        return self._start_idx

    @property
    def stop_index(self) -> int:
        return self._stop_idx


class SequenceIterator:

    def __init__(self, sequence: Sequence, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers
        if not hasattr(sequence, "_sequence_"):
            self._loggers.errorlog.log_raise_exception(
                f"Missing _sequence_ attribute on class {self.__class__.__name__}",
                os.EX_DATAERR,
                AttributeError,
            )
        self._sequence = sequence  # sequence object to iterate
        self._index = 0  # iterator index

    def __next__(self) -> str:
        if self._index < len(self._sequence):
            result = self._sequence[self._index]
            assert isinstance(result, str)  # slices not allowed here
            self._index += 1  # go to next position in sequence
            return result
        raise StopIteration  # stop iteration over sequence object


class Fasta:

    def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._fname = fname  # store input file
        self._faidx = self._index_fasta()  # initialize fasta index

    def _index_fasta(self) -> str:
        # look for index file for the current fasta file, if not found compute it
        if _find_fai(self._fname):  # index found, return it
            return f"{self._fname}.{FAI}"
        self._loggers.verboselog.debug(
            f"FASTA index not found for {self._fname}. Generating FASTA index"
        )
        start = time()  # measure indexing  time
        try:
            pysam.faidx(self._fname)  # index fasta using samtools
        except (SamtoolsError, Exception):
            self._loggers.errorlog.log_exception(
                f"An error occurred while indexing {self._fname}", os.EX_SOFTWARE
            )
        assert _find_fai(self._fname)  # index must be available now
        self._loggers.verboselog.debug(
            f"FASTA index for {self._fname} computed in {time() - start: .2f}s"
        )
        return f"{self._fname}.{FAI}"


def _find_fai(fastafile: str) -> bool:
    """Check if a FASTA index file exists for the given FASTA file.

    Checks if a FASTA index file (.fai) exists for the given FASTA file in the same
    directory.

    Args:
        fastafile: The path to the FASTA file.

    Returns:
        True if the index file exists and is not empty, False otherwise.
    """
    # search index for the input fasta file, assumes that the index is located
    # in the same folder as the indexed fasta
    fastaindex = f"{os.path.abspath(fastafile)}.{FAI}"  # avoid unexpected crashes
    if os.path.exists(fastaindex):
        return os.path.isfile(fastaindex) and os.stat(fastaindex).st_size > 0
    return False


class GenomeFasta(Fasta):

    def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
        super().__init__(fname, loggers)
        self._contig = self._retrieve_contig()  # initialize contig name
        self._length = self._compute_contig_length()  # initialize fasta length

    def __len__(self) -> int:
        return self._length

    def _retrieve_contig(self) -> str:
        f = pysam.FastaFile(self._fname, filepath_index=self._faidx)
        if len(f.references) != 1:  # fastas are expected to be chromosome-wise
            self._loggers.errorlog.log_raise_exception(
                f"Unexpected number of contigs ({len(f.references)}) found in {self._fname}",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        return f.references[0]  # contig name

    def _compute_contig_length(self) -> int:
        f = pysam.FastaFile(self._fname, filepath_index=self._faidx)
        if len(f.lengths) != 1:  # fastas are expected to be chromosome-wise
            self._loggers.errorlog.log_raise_exception(
                f"Unexpected number of contigs ({len(f.references)}) found in {self._fname}",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        return f.lengths[0]  # contig name

    def read(self) -> List[str]:
        try:  # read fasta file sequence content
            with open(self._fname, mode="r") as infile:
                infile.readline()  # skip fasta header
                # read fasta content
                # self._sequence = Sequence(
                #     "".join([line.strip() for line in infile.readlines()]),
                #     self._loggers,
                # )
                return list("".join([line.strip() for line in infile.readlines()]))
        except (IOError, Exception):
            self._loggers.errorlog.log_exception(
                f"An error occurred while reading {self._fname}", os.EX_IOERR
            )
            sys.exit(os.EX_IOERR)

    # @property
    # def sequence(self) -> Sequence:
    #     return self._sequence

    @property
    def contig(self) -> str:
        return self._contig


class GuideFasta(Fasta):

    def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
        super().__init__(fname, loggers)
        self._read_guides()  # read guides in fasta

    def _read_guides(self) -> None:
        f = pysam.FastaFile(self._fname, filepath_index=self._faidx)
        gnames = f.references  # retrieve guides seqnames (fasta headers)
        try:
            self._guides = list(
                {Sequence(f.fetch(gname), self._loggers) for gname in gnames}
            )  # extract guide sequence
        except (SamtoolsError, Exception):
            self._loggers.errorlog.log_exception(
                f"Failed parsing guides from {self._fname}", os.EX_DATAERR
            )

    @property
    def guides(self) -> List[Sequence]:
        return self._guides
