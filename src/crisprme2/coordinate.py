""" """

from .logger import CrisprmeLoggers

import os


class Coordinate:
    def __init__(
        self, contig: str, position: int, strand: int, loggers: CrisprmeLoggers
    ) -> None:
        self._loggers = loggers  # store loggers
        if not isinstance(contig, str) or not contig:
            self._loggers.errorlog.log_raise_exception(
                "Coordinate contig must be a non-empty string", os.EX_DATAERR, TypeError
            )
            raise ValueError("Contig must be a non-empty string.")
        self._contig: str = contig  # contig name
        if not isinstance(position, int):
            self._loggers.errorlog.log_raise_exception(
                f"Coordinate position must be an integer, got {type(position).__name__}",
                os.EX_DATAERR,
                TypeError,
            )
            raise TypeError("Position must be an integer.")
        if position < 0:
            self._loggers.errorlog.log_raise_exception(
                f"Coordinate position must be non-negative, got {position}",
                os.EX_DATAERR,
                TypeError,
            )
            raise ValueError("Position must be non-negative.")
        self._position: int = position
        if not isinstance(strand, int):
            self._loggers.errorlog.log_raise_exception(
                f"Coordinate strand must be an integer, got {type(strand).__name__}",
                os.EX_DATAERR,
                TypeError,
            )
            raise TypeError("Strand must be an integer.")
        # standard convention: 0 for FWD (+), 1 for REV (-)
        if strand not in {0, 1}:
            self._loggers.errorlog.log_raise_exception(
                f"Coordinate strand must be 0 (forward) or 1 (reverse), got {strand}",
                os.EX_DATAERR,
                TypeError,
            )
            raise ValueError("Strand must be 0 or 1.")
        self._strand: int = strand

    def __hash__(self) -> int:
        # Python's built-in hash function on a tuple of immutable components
        # is the standard and most robust way to implement object hashing.
        return hash((self._contig, self._position, self._strand))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Coordinate):
            return NotImplemented
        # compare all immutable fields that contribute to the unique identity
        return (
            self._contig == other._contig
            and self._position == other._position
            and self._strand == other._strand
        )

    @property
    def contig(self) -> str:
        return self._contig

    @property
    def position(self) -> int:
        return self._position

    @property
    def strand(self) -> int:
        return self._strand

    def __repr__(self) -> str:
        strand_char = "+" if self._strand == 0 else "-"
        return (
            f"<{self.__class__.__name__} object; contig='{self.contig}', pos="
            f"{self.position}, strand='{strand_char}')>"
        )

    def __str__(self) -> str:
        strand_char = "+" if self._strand == 0 else "-"
        return f"{self.contig}:{self.position}:{strand_char}"
