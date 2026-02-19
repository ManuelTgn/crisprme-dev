""" """

from .crisprme2_error import (
    Crisprme2SequenceError,
    Crisprme2ContigSequenceError,
    Crisprme2ReverseComplementError,
    Crisprme2DnaRnaError,
)
from .logger import CrisprmeLoggers
from .utils import RC

from typing import Optional, List

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
            self._loggers.errorlog.log_raise_exception(
                f"Invalid coordinates: start={start}, end={end}, length={len(self)}",
                os.EX_DATAERR,
                Crisprme2SequenceError,
            )
        return "".join(self._sequence[start:end])

    def reverse_complement(self) -> List[str]:
        try:
            return list(reverse_complement(self.sequence, self._loggers))
        except (KeyError, Exception) as e:
            self._loggers.errorlog.log_raise_exception(
                f"Error computing reverse complement: {str(e)}",
                os.EX_DATAERR,
                Crisprme2SequenceError,
            )

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

    def __init__(
        self,
        sequence: str,
        contig: str,
        start: int,
        stop: int,
        loggers: CrisprmeLoggers,
    ) -> None:
        super().__init__(sequence, loggers)
        self._contig = contig  # store sequence contig name
        self._start = start  # start sequence start position
        self._stop = stop  # store sequence stop position

    def chunk(self, size: int, overlap: int) -> List[str]:
        # validate size and overlap values
        if size <= 0:
            self._loggers.errorlog.log_raise_exception(
                f"Chunk size must be > 0 (got {size})",
                os.EX_DATAERR,
                Crisprme2ContigSequenceError,
            )
        if overlap < 0:
            self._loggers.errorlog.log_raise_exception(
                f"Overlap must be >= 0 (got {overlap})",
                os.EX_DATAERR,
                Crisprme2ContigSequenceError,
            )
        if overlap >= size and self._length > size:
            self._loggers.errorlog.log_raise_exception(
                f"Overlap size ({overlap}) must be less than the chunk size ({size})",
                os.EX_DATAERR,
                Crisprme2ContigSequenceError,
            )
        # handle simple cases: empty sequence or short sequence
        if self._length == 0:  # empty sequence
            return []
        # small contig: single chunk, no need to do anything fancy
        if size >= self._length:
            return [self.sequence]
        sequence = self.sequence  # materialize the full constig string once
        # number of core chunks
        n_chunks = (self._length + size - 1) // size  # ceil(n / size)
        # initialize chunk sequences list + preallocation
        chunks: List[str] = [None] * n_chunks  # type: ignore
        write_idx = 0
        for start in range(0, self._length, size):  # construct chunks
            if (stop := start + size) > self._length:
                stop = self._length
            ext_start = 0 if start == 0 else start - overlap
            ext_stop = stop  # ext_stop is stop (left-overlap)
            # slice directly from materialized string
            chunks[write_idx] = sequence[ext_start:ext_stop]
            write_idx += 1
            if stop == self._length:
                break  # nothing to do from here
        if write_idx != n_chunks:
            chunks = chunks[:write_idx]  # trash preallocated unused cells
        return chunks

    @property
    def contig(self) -> str:
        return self._contig

    @property
    def start(self) -> int:
        return self._start

    @property
    def stop(self) -> int:
        return self._stop


def reverse_complement(sequence: str, loggers: CrisprmeLoggers) -> str:
    try:
        return "".join([RC[nt] for nt in sequence[::-1]])
    except (KeyError, Exception) as e:
        loggers.errorlog.log_raise_exception(
            f"Failed reverse complement on {sequence}: {e}",
            os.EX_DATAERR,
            Crisprme2ReverseComplementError,
        )


def dna2rna(sequence: str, loggers: CrisprmeLoggers) -> str:
    try:
        return sequence.replace("T", "U").replace("t", "u")
    except ValueError as e:
        loggers.errorlog.log_raise_exception(
            f"Failed translating DNA to RNA on {sequence}: {e}",
            os.EX_DATAERR,
            Crisprme2DnaRnaError,
        )
