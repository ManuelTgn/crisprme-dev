""" """

from .crisprme2_error import Crisprme2SequenceError, Crisprme2ContigSequenceError
from .logger import CrisprmeLoggers
from .encoder import BitSequence
from .utils import RC

from typing import Optional, Union, List, Generator

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
            return [RC[nt] for nt in self._sequence[::-1]]
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

    def chunk(self, size: int, overlap: int):
        if overlap >= size:
            self._loggers.errorlog.log_raise_exception(
                f"Overlap size ({overlap}) must be less than the chunk size ({size})",
                os.EX_DATAERR,
                Crisprme2ContigSequenceError,
            )

        if size >= self._length:  # small contig
            yield ContigSequence(
                self.subsequence(0, self._length),
                self._contig,
                0,
                self._length,
                self._loggers,
            )
            return

        step = size - overlap  # ← FIX: move forward by size−overlap

        for i in range(0, self._length, step):
            start = i
            stop = i + size

            # clip to contig boundaries
            if stop > self._length:
                stop = self._length

            yield ContigSequence(
                self.subsequence(start, stop), 
                self._contig, 
                start, 
                stop, 
                self._loggers
            )

    def encode(self) -> bytearray:
        # encode contig sequence as byte array
        return BitSequence(self.sequence, self._loggers).data


    @property
    def contig(self) -> str:
        return self._contig

    @property
    def start(self) -> int:
        return self._start

    @property
    def stop(self) -> int:
        return self._stop


