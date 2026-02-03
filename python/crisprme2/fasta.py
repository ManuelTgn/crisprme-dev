""" """

from .crisprme2_error import (
    Crisprme2FastaError,
    Crisprme2FastaFileNotFoundError,
    Crisprme2SequenceError,
)
from .utils import FAI, warning, find_fai_index
from .sequence import ContigSequence
from .logger import CrisprmeLoggers

from typing import Optional, List
from pysam.utils import SamtoolsError
from pysam import FastaFile, faidx
from pathlib import Path

import os


# fasta file extensions
FASTAEXTENSIONS = valid_extensions = {"fasta", "fa", "fna", "ffn", "faa", "frn", "fas"}


class Fasta:

    def __init__(self, filepath: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._filepath = filepath  # fasta filename
        self._validate_file()  # validate fasta file structure
        self._index = self._search_index()  # fai index
        self._fasta_handle: Optional[FastaFile] = None
        self._is_open = False
        self._init_contig_length()  # initialize contig name and length

    def _validate_file(self) -> None:
        # check file extension
        if (
            os.path.splitext(os.path.basename(self._filepath))[1].replace(".", "")
            not in FASTAEXTENSIONS
        ):
            self._loggers.errorlog.log_raise_exception(
                f"File {self._filepath} does not have a standard FASTA extension",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )

    def _index_fasta(self, pytest: bool = False) -> str:
        if hasattr(self, "_index"):
            if self._index and not pytest:  # launch warning
                warning("FASTA index already present, forcing update", 1)
        try:  # create index in the same folder as the input fasta
            self._loggers.verboselog.debug(
                f"Creating index for FASTA: {self._filepath}"
            )
            faidx(str(self._filepath))
        except (OSError, Exception):
            self._loggers.errorlog.log_exception(
                f"Failed indexing for FASTA: {self._filepath}", os.EX_DATAERR
            )
        assert find_fai_index(str(self._filepath))  # now should be available
        return f"{self._filepath}.{FAI}"

    def _search_index(self) -> Path:
        # look for index for the current fasta, if not found compute it
        if find_fai_index(str(self._filepath)):  # index found, store it
            return Path(f"{self._filepath}.{FAI}")
        # index not found -> compute it de novo and store it in the same folder
        # as the input fasta
        self._loggers.verboselog.debug(f"FASTA index not found for {self._filepath}")
        return Path(self._index_fasta())

    def _init_contig_length(self) -> None:
        self.open()  # manually open fasta file
        assert self._fasta_handle  # should be open
        self._ncontigs = len(self._fasta_handle.references)
        if self._ncontigs != 1:
            self._loggers.errorlog.log_raise_exception(
                f"Multiple contigs found in {self._filepath}",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        # we're 100% sure that there is only one contig in this fasta
        self._contig = self._fasta_handle.references[0]
        self._length = self._fasta_handle.lengths[0]
        self.close()  # manually close fasta file

    def open(self) -> "Fasta":
        if self._is_open:
            self._loggers.errorlog.log_raise_exception(
                f"FASTA file {self._filepath} is already open",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        try:  # open fasta, assumes that index is already available
            self._fasta_handle = FastaFile(str(self._filepath))
            self._is_open = True
        except (OSError, Exception) as e:
            self._loggers.errorlog.log_exception(
                f"Failed to open FASTA file {self._filepath}: {str(e)}", os.EX_IOERR
            )
        return self

    def close(self) -> None:
        if self._fasta_handle is not None:
            self._fasta_handle.close()
            self._is_open = False

    def __enter__(self) -> "Fasta":
        return self.open()

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def read(self):
        with open(self._filepath, mode="rb") as fin:
            fin.readline()  # consume header buffer
            return bytearray(fin.read())

    def fetch(
        self, reference: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> ContigSequence:
        if not self._is_open or self._fasta_handle is None:
            self._loggers.errorlog.log_raise_exception(
                "FASTA file must be opened before fetching",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        assert self._fasta_handle  # must not be none
        try:
            if start is None and end is None:  # access string by contig name
                return ContigSequence(
                    self._fasta_handle.fetch(reference),
                    self._contig,
                    0,
                    self._length,
                    self._loggers,
                )
            elif start is not None and end is not None:
                if start < 0 or end < start:
                    self._loggers.errorlog.log_raise_exception(
                        f"Invalid coordinates: start={start}, end={end}",
                        os.EX_DATAERR,
                        Crisprme2SequenceError,
                    )
                return ContigSequence(
                    self._fasta_handle.fetch(reference, start, end),
                    self._contig,
                    start,
                    end,
                    self._loggers,
                )
            else:
                self._loggers.errorlog.log_raise_exception(
                    "Both start and end must be specified or both None",
                    os.EX_DATAERR,
                    Crisprme2SequenceError,
                )
        except KeyError:
            self._loggers.errorlog.log_raise_exception(
                f"Reference '{reference}' not found in FASTA file",
                os.EX_DATAERR,
                Crisprme2FastaError,
            )
        except Exception as e:
            self._loggers.errorlog.log_exception(
                f"Error fetching sequence: {str(e)}", os.EX_DATAERR
            )

    def __contains__(self, reference: str) -> bool:
        return (
            reference in self._fasta_handle.references
            if self._is_open and self._fasta_handle
            else False
        )

    def __repr__(self) -> str:
        status = "open" if self._is_open else "closed"
        return f"<{self.__class__.__name__} object; sequences={self._contig}, status={status}>"

    def __del__(self):
        if self._is_open:
            self.close()

    @property
    def contig(self) -> str:
        return self._contig if self._contig.startswith("chr") else f"chr{self._contig}"

    @property
    def length(self) -> int:
        return self._length  # return contig lengths in fasta

    @property
    def nreferences(self) -> int:
        return self._ncontigs  # return the number of contigs in input fasta


class GuideFasta(Fasta):

    def __init__(self, filepath: str, loggers: CrisprmeLoggers) -> None:
        super().__init__(filepath, loggers)

    def _read_guides(self) -> None:
        f = FastaFile(str(self._filepath), filepath_index=str(self._index))
        gnames = f.references  # retrieve guides seqnames (fasta headers)
        try:  # extract guide sequence
            self._guides = list({f.fetch(gname) for gname in gnames})
        except (SamtoolsError, Exception):
            self._loggers.errorlog.log_exception(
                f"Failed parsing guides from {self._filepath}", os.EX_DATAERR
            )

    @property
    def guides(self) -> List[str]:
        return self._guides
