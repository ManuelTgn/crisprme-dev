""" """

from .crisprme2_error import Crisprme2TargetError
from .coordinate import Coordinate
from .logger import CrisprmeLoggers

from typing import List, Set

import os


class Target:

    def __init__(self, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers  # store loggers
        self._coordinates: Set[Coordinate] = set()  # target positions across the genome
        self._alignments: List[bytes] = (
            []
        )  # target alignments stored as raw CIGAR ASCII

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; coordinate={len(self.coordinates)}>"

    def add_target(self, contig: str, position: int, strand: bool) -> None:
        # NOTE on Strand Conversion:
        # The Rust binding uses True=Forward, False=Reverse.
        # It's conventional in bioinformatics/internal code to use 0 for FWD (+) and 1 for REV (-).
        # We invert the boolean (True=0, False=1) to follow this convention.
        # Rust's True (fwd) -> Python's 0 (fwd)
        # Rust's False (rev) -> Python's 1 (rev)
        self._coordinates.add(
            Coordinate(contig, position + 1, int(not strand), self._loggers)
        )

    @property
    def coordinates(self) -> Set[Coordinate]:
        return self._coordinates

    @property
    def alignments(self) -> List[bytes]:
        return self._alignments
