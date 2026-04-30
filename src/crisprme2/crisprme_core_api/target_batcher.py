""" """

from .crisprme2_api_error import Crisprme2BatcherError
from ..logger import CrisprmeLoggers

try:  # import rust API modules
    from .._crisprme2_native import TargetBatcher as RustTargetBatcher
    from .._crisprme2_native import BatcherStats, FeedStatus
except ImportError:
    # fallback for development/testing
    RustTargetBatcher = None
    BatcherStats = None
    FeedStatus = None

from typing import Optional, Tuple, List, Iterator
from dataclasses import dataclass

import os


# define number of batch targets collected and sent to aligner
BATCHITS = 1_000_000

# define maximum number of unique targets before flushing
BATCHMAXUNQ = 500_000


@dataclass
class BatchStats:
    hits_in_batch: int
    unique_windows: int

    @classmethod
    def from_rust(cls, stats: "BatcherStats") -> "BatchStats":
        return cls(
            hits_in_batch=stats.hits_in_batch, unique_windows=stats.unique_windows
        )


@dataclass
class FeedResult:
    flushed: bool
    stats: BatchStats

    @classmethod
    def from_rust(cls, status: "FeedStatus") -> "FeedResult":
        return cls(flushed=status.flushed, stats=status.stats)


class TargetBatcher:
    def __init__(
        self,
        pam: str,
        guide: str,
        size: int,
        right: bool,
        overlap: int,
        threads: int,
        loggers: CrisprmeLoggers,
    ) -> None:
        self._loggers = loggers  # store loggers
        if RustTargetBatcher is None:
            self._loggers.errorlog.log_raise_exception(
                "Rust TargetBatcher module not exposed to python",
                os.EX_CANTCREAT,
                ValueError,
            )
        _validate_overlap(size, overlap, self._loggers)  # validate overlap
        self._batcher = RustTargetBatcher(
            pam, guide, size, right, threads, BATCHITS, BATCHMAXUNQ, overlap
        )
        # initialize chunks stats
        self._total_chunks_fed = 0
        self._total_flushes = 0

    def __repr__(self) -> str:
        stats = self._get_stats()
        return (
            f"<{self.__class__.__name__} object; id={self.id}, "
            f"hits={stats.hits_in_batch}, unique={stats.unique_windows}>"
        )
    
    @property
    def batcher(self) -> RustTargetBatcher:
        return self._batcher

    @property
    def id(self) -> int:
        return self._batcher.id

    def _get_stats(self) -> BatchStats:
        stats = self._batcher.stats()
        return BatchStats.from_rust(stats)

    def feed_chunk(
        self, contig_id: int, chunk_start: int, chunk_seq: str, valid_len: int
    ) -> FeedResult:
        self._total_chunks_fed += 1  # increase number of actual chunks fed
        # feed current sequence chunk to batcher to collect targets
        result = self._batcher.feed_chunk(contig_id, chunk_start, chunk_seq, valid_len)
        # read batcher status
        feed_result = FeedResult.from_rust(result)
        if feed_result.flushed:
            self._loggers.verboselog.debug(
                f"Batch flush triggered: {feed_result.stats.hits_in_batch} hits, {feed_result.stats.unique_windows} unique windows"
            )
        return feed_result

    def flush_and_align(
        self, max_mismatches: int, max_bdna: int, max_brna: int
    ) -> None:
        # NOTE This method is currently just a placeholder
        # TODO: replace current placeholder with actual flush_and_align
        self._total_flushes += 1  # increase number of executed flushes
        # write logging details before flushing
        stats = self._get_stats()
        self._loggers.verboselog.debug(
            f"Flushing batch #{self._total_flushes}: {stats.hits_in_batch} hits, {stats.unique_windows} unique windows"
        )
        # batcher filled, send off-target candidates to aligner
        self._batcher.flush_and_align(max_mismatches, max_bdna, max_brna)

    def finalize(self) -> BatchStats:
        # read final batcher stats and write details on logger
        stats = self._batcher.finalize()
        stats = BatchStats.from_rust(stats)
        self._loggers.basiclog.info(
            f"Finalized {self.__class__.__name__} (id={self.id}): processed {self._total_chunks_fed} chunks, {self._total_flushes} flushes"
        )
        return stats


def _validate_overlap(size: int, overlap: int, loggers: CrisprmeLoggers) -> int:
    # validate overlap size
    if size > 0 and overlap < size - 1:
        loggers.errorlog.log_raise_exception(
            f"Overlap is less than size - 1 ({overlap} < {size - 1})",
            os.EX_DATAERR,
            Crisprme2BatcherError,
        )
    return overlap
