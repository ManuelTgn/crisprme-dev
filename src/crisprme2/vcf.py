""" """

from .crisprme2_error import (
    Crisprme2VCFError,
    Crisprme2VCFFileNotFoundError,
    Crisprme2VCFFormatError,
    Crisprme2VCFIndexError,
)
from .logger import CrisprmeLoggers
from .utils import find_tbi_index, warning, TBI

from typing import Optional, Iterator, List, Dict, Any
from pysam import TabixFile, TabixIterator, tabix_index
from pathlib import Path

import sys
import os

# vcf extensions
VCFEXTENSIONS = {"vcf", "vcf.gz", "bcf", "bcf.gz"}


class VCF:

    def __init__(self, filepath: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._filepath = Path(filepath)  # vcf filename
        self._validate_file()  # validate vcf file structure
        self._index = self._search_index()  # tbi index
        self._vcf_handle: Optional[TabixFile] = None
        self._is_open = False

    def __repr__(self) -> str:
        status = "open" if self._is_open else "closed"
        return f"<{self.__class__.__name__} object; file={self._filepath}, status={status}>"

    def __del__(self):
        if self._is_open:
            self.close()

    def __enter__(self) -> None:
        return self.open()

    def __exit__(self) -> None:
        self.close()

    def _index_vcf(self, pytest: bool = False) -> str:
        if self._index and not pytest:  # launch warning
            warning("Tabix index already present, forcing update", 1)
        try:  # create index in the same folder as the input vcf
            tabix_index(str(self._filepath), preset="vcf", force=True)
        except (OSError, Exception):
            self._loggers.errorlog.log_exception(
                f"Failed indexing for VCF: {self._filepath}", os.EX_DATAERR
            )
        assert find_tbi_index(str(self._filepath))  # now should be available
        return f"{self._filepath}.{TBI}"

    def _search_index(self) -> Path:
        # look for index for the current vcf, if not found compute it
        if find_tbi_index(str(self._filepath)):  # index found, store it
            return Path(f"{self._filepath}.{TBI}")
        # index not found -> compute it de novo and store it in the same folder
        # as the input vcf
        self._loggers.verboselog.debug(f"Tabix index not found for {self._filepath}")
        return Path(self._index_vcf())

    def _validate_file(self) -> None:
        if not self._filepath.exists():
            self._loggers.errorlog.log_raise_exception(
                f"VCF file not found: {self._filepath}",
                os.EX_DATAERR,
                Crisprme2VCFFileNotFoundError,
            )
        if not self._filepath.is_file():
            self._loggers.errorlog.log_raise_exception(
                f"Path is not a file: {self._filepath}",
                os.EX_DATAERR,
                Crisprme2VCFFileNotFoundError,
            )
        # check file extension
        if not any(str(self._filepath).endswith(ext) for ext in VCFEXTENSIONS):
            self._loggers.errorlog.log_raise_exception(
                f"File {self._filepath} does not have a standard VCF extension",
                os.EX_DATAERR,
                Crisprme2VCFError,
            )

    def open(self) -> None:
        if self._is_open:  # vcf already open
            self._loggers.errorlog.log_raise_exception(
                f"VCF file {self._filepath} is already open",
                os.EX_IOERR,
                Crisprme2VCFError,
            )
        self._loggers.verboselog.debug(f"Opening VCF file {self._filepath}")
        try:  # open vcf, assumes that index is already available
            self._vcf_handle = TabixFile(str(self._filepath), index=str(self._index))
            self._is_open = True
        except (OSError, Exception) as e:
            self._loggers.errorlog.log_exception(
                f"Failed to open VCF file {self._filepath}: {e}", os.EX_IOERR
            )
        self._loggers.verboselog.debug(
            f"Successfully opened VCF file: {self._filepath}"
        )

    def close(self) -> None:
        if self._vcf_handle is not None and self._is_open:
            self._vcf_handle.close()
            self._is_open = False
            self._loggers.verboselog.debug(f"Closed VCF file: {self._filepath}")

    def get_samples(self) -> List[str]:
        if not self._is_open or self._vcf_handle is None:
            self._loggers.errorlog.log_raise_exception(
                f"VCF file not open: {self._filepath}; cannot retrieve samples",
                os.EX_IOERR,
                Crisprme2VCFError,
            )
        assert self._vcf_handle is not None
        return self._vcf_handle.header[-1].strip().split()[9:]

    def fetch(
        self,
        contig: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> TabixIterator:
        if not self._is_open or self._vcf_handle is None:
            self._loggers.errorlog.log_raise_exception(
                "VCF file must be opened before fetching variants",
                os.EX_IOERR,
                Crisprme2VCFError,
            )
        assert self._vcf_handle is not None
        try:  # fetch variants in range
            if contig is not None:
                if start is not None or end is not None:
                    assert self._index  # region-based fetch requires index
                return self._vcf_handle.fetch(contig, start, end)
            else:
                return self._vcf_handle.fetch()
        except Exception as e:
            self._loggers.errorlog.log_exception(
                f"Error fetching variants from {self._filepath}: {str(e)}",
                os.EX_DATAERR,
            )
        sys.exit(os.EX_IOERR)

    def count_variants(
        self,
        contig: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> int:
        print("hello")
        return sum(1 for _ in self.fetch(contig, start, end))

    @property
    def contigs(self) -> List[str]:
        return (
            list(self._vcf_handle.contigs)
            if (self._is_open and self._vcf_handle is not None)
            else []
        )
