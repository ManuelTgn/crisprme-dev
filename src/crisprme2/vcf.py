""" """

from .crisprme2_error import (
    Crisprme2VCFError,
    Crisprme2VCFFileNotFoundError,
    Crisprme2VCFFormatError,
)
from .logger import CrisprmeLoggers
from .variant import VariantRecord
from .utils import find_tbi_index, warning, TBI

from typing import Optional, List
from pysam import TabixFile, TabixIterator, tabix_index
from pathlib import Path

import cyvcf2
import os


# vcf extensions
VCFEXTENSIONS = {"vcf", "vcf.gz", "bcf", "bcf.gz"}


class VCF:

    def __init__(self, filepath: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._filepath = Path(filepath)  # vcf filename
        self._validate_file()  # validate vcf file structure
        self._index = self._search_index()  # tbi index
        self._assess_phasing()  # assess if VCF is phased
        self._contig: str = cyvcf2.VCF(str(self._filepath)).seqnames[0]  # contig name

    def _validate_file(self) -> None:
        # check file existence and extension
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
        if not any(str(self._filepath).endswith(ext) for ext in VCFEXTENSIONS):
            self._loggers.errorlog.log_raise_exception(
                f"File {self._filepath} does not have a standard VCF extension",
                os.EX_DATAERR,
                Crisprme2VCFFormatError,
            )
        # check number of contigs and genotyping availability in VCF
        vcf_ = cyvcf2.VCF(str(self._filepath))
        if len(vcf_.seqnames) > 1:
            self._loggers.errorlog.log_raise_exception(
                f"Multiple contigs in {self._filepath}",
                os.EX_DATAERR,
                Crisprme2VCFFormatError,
            )
        if not vcf_.contains("GT"):
            self._loggers.errorlog.log_raise_exception(
                f"Missing genotype (GT) field in {self._filepath}",
                os.EX_DATAERR,
                Crisprme2VCFFormatError,
            )

    def _index_vcf(self, pytest: bool = False) -> str:
        if not pytest:  # launch warning
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

    def _assess_phasing(self) -> None:
        variant = None
        for v in cyvcf2.VCF(self._filepath):
            variant = v
            break
        assert variant is not None
        self._phasing = all(phase for phase in variant.gt_phases)

    def read(
        self, start: Optional[int], stop: Optional[int], threads: int = 1
    ) -> List[VariantRecord]:
        reader = cyvcf2.VCF(str(self._filepath), lazy=True, threads=threads)  # open vcf
        if start is not None and stop is not None:
            region = f"{self._contig}:{start}-{stop}"
            return [VariantRecord(v, self._loggers) for v in reader(region)]
        return [VariantRecord(v, self._loggers) for v in reader]

    @property
    def contig(self) -> str:
        return self._contig if self._contig.startswith("chr") else f"chr{self._contig}"

    def get_samples(self) -> List[str]:
        return cyvcf2.VCF(str(self._filepath)).samples

    @property
    def filepath(self) -> str:
        return str(self._filepath)

    @property
    def index(self) -> str:
        return str(self._index)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; contigs={self.contig}>"
