""" """

from .crisprme2_error import Crisprme2AnnotationError
from .logger import CrisprmeLoggers
from .utils import find_tbi_index, TBI

from typing import Optional, Literal, TypeAlias, cast
from pysam.utils import SamtoolsError
from pysam import tabix_index, TabixFile

import os


# define tabix preset types
TabixPreset: TypeAlias = Literal["gff", "bed", "sam", "vcf", "psltbl", "pileup"]
_ALLOWED: set[str] = {"gff", "bed", "sam", "vcf", "psltbl", "pileup"}


class Annotation:

    def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # set crisprme loggers
        self._fname = fname  # set annotation file name
        self._is_open = False

    def _search_index(self, preset: str, tbx: Optional[str] = "") -> str:
        # look for index for the current bedfile, if not found compute it
        if not tbx:
            if find_tbi_index(self._fname):  # index found, return it
                return f"{self._fname}.{TBI}"
            # index not found -> compute it de novo and store it
            self._loggers.verboselog.debug(f"Tabix index not found for {self._fname}\n")
            return _tabix_index(self._fname, preset, self._loggers)
        return tbx


class AnnotationBed(Annotation):

    def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
        super().__init__(fname, loggers)  # initialize annotation object
        self._bedidx = self._search_index("bed")  # set bed tabix index
        self._bed: Optional[TabixFile] = None

    def open(self) -> "AnnotationBed":
        if self._is_open:
            self._loggers.errorlog.log_raise_exception(
                f"Annottaion BED file {self._fname} is already open",
                os.EX_IOERR,
                Crisprme2AnnotationError,
            )
        try:  # open bed, assumes that index is already available
            self._bed = TabixFile(self._fname, index=self._bedidx)
            self._is_open = True
        except (IOError, Exception) as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed opening {self._fname}: {e}",
                os.EX_IOERR,
                Crisprme2AnnotationError,
            )
        return self

    def close(self) -> None:
        if self._bed is not None:
            self._bed.close()
            self._is_open = False

    def __enter__(self) -> "AnnotationBed":
        return self.open()

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def __repr__(self) -> str:
        status = "open" if self._is_open else "closed"
        return f"<{self.__class__.__name__} object; status={status}>"

    def __del__(self):
        if self._is_open:
            self.close()

    def fetch_features(self, contig: str, start: int, stop: int) -> Optional[str]:
        if not self._is_open or self._bed is None:
            self._loggers.errorlog.log_raise_exception(
                "Annotation BED file must be opened before fetching",
                os.EX_DATAERR,
                Crisprme2AnnotationError,
            )
        if contig not in self._bed.contigs:
            return None
        return ",".join([e.strip() for e in self._bed.fetch(contig, start, stop)])


class AnnotationGff(Annotation):

    def __init__(self, fname: str, loggers: CrisprmeLoggers) -> None:
        super().__init__(fname, loggers)  # initialize annotation object
        self._gffidx = self._search_index("gff")  # set gff tabix index
        self._gff: Optional[TabixFile] = None

    def open(self) -> "AnnotationGff":
        if self._is_open:
            self._loggers.errorlog.log_raise_exception(
                f"Annotation GFF file {self._fname} is already open",
                os.EX_IOERR,
                Crisprme2AnnotationError,
            )
        try:  # open gff, assumes that index is already available
            self._gff = TabixFile(self._fname, index=self._gffidx)
            self._is_open = True
        except (IOError, Exception) as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed opening {self._fname}: {e}",
                os.EX_IOERR,
                Crisprme2AnnotationError,
            )
        return self

    def close(self) -> None:
        if self._gff is not None:
            self._gff.close()
            self._is_open = False

    def __enter__(self) -> "AnnotationGff":
        return self.open()

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def __repr__(self) -> str:
        status = "open" if self._is_open else "closed"
        return f"<{self.__class__.__name__} object; status={status}>"

    def __del__(self):
        if self._is_open:
            self.close()

    def fetch_features(self, contig: str, start: int, stop: int) -> Optional[str]:
        if not self._is_open or self._gff is None:
            self._loggers.errorlog.log_raise_exception(
                "Annotation GFF file must be opened before fetching",
                os.EX_DATAERR,
                Crisprme2AnnotationError,
            )
        if contig not in self._gff.contigs:
            return None
        return ",".join([e.strip() for e in self._gff.fetch(contig, start, stop)])


def _tabix_index(
    fname: str, preset: str, loggers: CrisprmeLoggers, force: bool = True
) -> str:
    if preset is not None and preset not in _ALLOWED:
        loggers.errorlog.log_raise_exception(
            f"Invalid tabix preset: {preset!r}. Must be one of {sorted(_ALLOWED)}",
            os.EX_DATAERR,
            Crisprme2AnnotationError,
        )
    # compute tabix index if not provided during annotation object initialization
    try:
        tabix_index(fname, preset=cast(TabixPreset, preset), force=force)
    except (SamtoolsError, Exception) as e:
        loggers.errorlog.log_raise_exception(
            f"Failed indexing {fname}: {e}", os.EX_DATAERR, Crisprme2AnnotationError
        )
    assert find_tbi_index(fname)
    return f"{fname}.{TBI}"
