""" """

from .crisprme2_error import Crisprme2FastaError, Crisprme2FastaFileNotFoundError, Crisprme2SequenceError
from .utils import FAI, warning, find_fai_index
from .sequence import Sequence
from .logger import CrisprmeLoggers

from typing import Optional, List
from pysam import FastaFile, faidx
from pathlib import Path

import sys
import os

# fasta file extensions
FASTAEXTENSIONS = valid_extensions = {"fasta", "fa", "fna", "ffn", "faa", "frn", "fas"}

class Fasta:
    
    def __init__(self, filepath: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._filepath = Path(filepath)  # fasta filename
        self._validate_file()  # validate fasta file structure
        self._index = self._search_index()  # fai index
        self._fasta_handle: Optional[FastaFile] = None
        self._is_open = False  
        self._contig = None      
    
    def _validate_file(self) -> None:
        if not self._filepath.exists():
            self._loggers.errorlog.log_raise_exception(f"FASTA file not found: {self._filepath}", os.EX_DATAERR, Crisprme2FastaFileNotFoundError)
        if not self._filepath.is_file():
            self._loggers.errorlog.log_raise_exception(f"Path is not a file: {self._filepath}", os.EX_DATAERR, Crisprme2FastaFileNotFoundError)        
        # check file extension
        if not any(str(self._filepath).lower().endswith(ext) for ext in FASTAEXTENSIONS):
            self._loggers.errorlog.log_raise_exception(f"File {self._filepath} does not have a standard FASTA extension", os.EX_DATAERR, Crisprme2FastaError)
    
    def _index_fasta(self, pytest: bool = False) -> str:
        if self._index and not pytest:  # launch warning
            warning("FASTA index already present, forcing update", 1)
        try:  # create index in the same folder as the input fasta
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
    
    def open(self) -> "Fasta":
        if self._is_open:
            self._loggers.errorlog.log_raise_exception(f"FASTA file {self._filepath} is already open", os.EX_DATAERR, Crisprme2FastaError)
        try:  # open fasta, assumes that index is already available
            self._fasta_handle = FastaFile(str(self._filepath))
            self._is_open = True
            self._contig = self.references[0]  # contig name
        except (OSError, Exception) as e:
            self._loggers.errorlog.log_exception(f"Failed to open FASTA file {self._filepath}: {str(e)}", os.EX_IOERR)
        return self
        
    
    def close(self) -> None:
        if self._fasta_handle is not None:
            self._fasta_handle.close()
            self._is_open = False
    
    def __enter__(self) -> "Fasta":
        return self.open()
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()
    
    @property
    def references(self) -> List[str]:
        if not self._is_open or self._fasta_handle is None:
            self._loggers.errorlog.log_raise_exception("FASTA file must be opened before accessing references", os.EX_DATAERR, Crisprme2FastaError)
        assert self._fasta_handle  # must not be none
        return list(self._fasta_handle.references)  # return contig names in fasta
    
    @property
    def length(self) -> int:
        if not self._is_open or self._fasta_handle is None:
            self._loggers.errorlog.log_raise_exception("FASTA file must be opened before accessing lengths", os.EX_DATAERR, Crisprme2FastaError)
        assert self._fasta_handle  # must not be none
        return self._fasta_handle.lengths[0]  # return contig lengths in fasta
    
    @property
    def nreferences(self) -> int:
        # return the number of contigs in input fasta
        return len(self.references) if self._is_open else 0
    
    
    def fetch(self, reference: str, start: Optional[int] = None, end: Optional[int] = None) -> Sequence:
        if not self._is_open or self._fasta_handle is None:
            self._loggers.errorlog.log_raise_exception("FASTA file must be opened before fetching", os.EX_DATAERR, Crisprme2FastaError)
        assert self._fasta_handle  # must not be none
        try:
            if start is None and end is None:  # access string by contig name
                return Sequence(self._fasta_handle.fetch(reference), self._loggers)
            elif start is not None and end is not None:
                if start < 0 or end < start:
                    self._loggers.errorlog.log_raise_exception(f"Invalid coordinates: start={start}, end={end}", os.EX_DATAERR, Crisprme2SequenceError)
                return Sequence(self._fasta_handle.fetch(reference, start, end), self._loggers)
            else:
                self._loggers.errorlog.log_raise_exception("Both start and end must be specified or both None", os.EX_DATAERR, Crisprme2SequenceError)
        except KeyError:
            self._loggers.errorlog.log_raise_exception(f"Reference '{reference}' not found in FASTA file", os.EX_DATAERR, Crisprme2FastaError)
        except Exception as e:
            self._loggers.errorlog.log_exception(f"Error fetching sequence: {str(e)}", os.EX_DATAERR)
        sys.exit(1)  # base case

    def __contains__(self, reference: str) -> bool:
        return reference in self.references if self._is_open else False
    
    def __repr__(self) -> str:
        status = "open" if self._is_open else "closed"
        return f"<{self.__class__.__name__} object; sequences={self._contig}, status={status}>"
    
    def __del__(self):
        if self._is_open:
            self.close()


