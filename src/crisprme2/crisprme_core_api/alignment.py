""" """

from .crisprme_api_error import Crisprme2AlignmentError
from ..logger import CrisprmeLoggers

from typing import List

import os


class Alignment:

    def __init__(
        self,
        occurrences: List[int],
        strand: List[int],
        guide_al: str,
        target_al: str,
        mm: int,
        bdna: int,
        brna: int,
        edit_d: int,
        cfd_score: float,
        annotations: List[bytearray],
        loggers: CrisprmeLoggers,
    ) -> None:
        self._loggers = loggers  # set crisprme logger (local to python)
        self._occurrences = occurrences  # sequence occurrence positions
        self._strand = strand  # occurrence strandness (1 | 0)
        self._guide_al = guide_al  # aligned guide sequence
        self._target_al = target_al  # aligned target sequence
        self._mm = mm  # number of mismatches in current alignment
        self._bdna = bdna  # number of dna bulges in current alignment
        self._brna = brna  # number of rna bulges in current alignment
        self._edit_d = edit_d  # edit distance (guide -> target)
        # NOTE for each score considered add a different attribute
        # at init time the score at rust level should be nan
        # eventually it will be computed up in python and assigned
        # through property assignment
        self._cfd_score = cfd_score
        # annotations supported as vector v (|v| = 10) of bytearrays (32 bits -> at most 32 different features)
        self._annotations = annotations

    # the following properties are required to access guide and target alignment
    # views to assign the score
    @property
    def guide(self) -> str:
        return self._guide_al

    @property
    def offtarget(self) -> str:
        return self._target_al

    @property
    def cfd_score(self) -> float:
        return self._cfd_score

    @cfd_score.setter
    def cfd_score(self, value: float) -> None:
        if not isinstance(value, float):
            self._loggers.errorlog.log_raise_exception(
                f"CFD score must be a float, got {type(value).__name__}",
                os.EX_DATAERR,
                Crisprme2AlignmentError,
            )
        self._cfd_score = value
