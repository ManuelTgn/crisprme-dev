"""
target_batcher.py
-----------------
Python wrapper for the Rust ``TargetBatcher`` struct exposed via PyO3.

Responsibilities
~~~~~~~~~~~~~~~~
- Validate every argument before crossing the FFI boundary.
- Expose typed, documented accessors for the fields that matter at the
  Python orchestration level (id, stats, flush signal).
- Translate ``FeedStatus.flushed`` into a clear return value so
  ``scanner.py`` can decide whether to call ``pipeline.submit()``
  without inspecting internal Rust state.
- Map all Rust ``PyValueError`` / ``PyResult`` failures to the typed
  ``Crisprme2BatcherError`` so callers get predictable exception types.

What this class does NOT do
~~~~~~~~~~~~~~~~~~~~~~~~~~~
- It does not own the scanning loop — that lives in ``scanner.py``.
- It does not decide alignment thresholds — those come from
  :class:`~crisprme2.crisprme_core_api.Thresholds` and are passed
  through at flush time.
- It does not hold a reference to the pipeline — submission is
  orchestrated externally.

Rust -> Python type mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    BatcherStats  ->  BatchStats   (dataclass, immutable)
    FeedStatus    ->  FeedResult   (dataclass, immutable)
    TargetBatcher ->  TargetBatcher (this wrapper class)

Typical usage (from scanner.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    from crisprme2.crisprme_core_api import TargetBatcher, Thresholds

    thresholds = Thresholds(max_mm=4, max_bdna=1, max_brna=1, loggers=loggers)
    batcher = TargetBatcher.create(
        pam=pam, guide=guide, size=30, upstream=True,
        overlap=29, threads=8, loggers=loggers,
    )

    result = batcher.feed_chunk(contig_id=0, chunk_start=0,
                                chunk_seq=seq, valid_len=len(seq))
    if result.flushed:
        pipeline.submit(batcher)          # drains internal map
        batcher.flush_and_align(thresholds)  # not normally called directly

    batcher.finalize()
"""

from __future__ import annotations

from typing import Any

from .crisprme2_api_error import Crisprme2BatcherError
from ..guide import Guide
from ..logger import CrisprmeLoggers
from ..pam import PAM

try:  # import rust API modules
    from .._crisprme2_native import TargetBatcher as RustTargetBatcher
    from .._crisprme2_native import BatcherStats, FeedStatus
except ImportError:
    # fallback for development/testing
    RustTargetBatcher = None
    BatcherStats = None
    FeedStatus = None

from dataclasses import dataclass

import os


# ------------------------------------------------------------------------------
# constants (mirror Rust defaults; override via TargetBatcher.create kwargs)
# ------------------------------------------------------------------------------

# Maximum total hits accumulated before an automatic flush is triggered.
_DEFAULT_BATCH_HITS: int = 100_000

# Maximum number of *unique* windows before an automatic flush is triggered.
_DEFAULT_MAX_UNIQUE: int = 500_000

# ------------------------------------------------------------------------------
# immutable result dataclasses
# ------------------------------------------------------------------------------


@dataclass(frozen=True)
class BatchStats:
    """
    Snapshot of a :class:`TargetBatcher`'s internal counters.

    Attributes
    ----------
    hits_in_batch : int
        Total occurrence count accumulated since the last flush (including
        duplicate windows at different genomic positions).
    unique_windows : int
        Number of distinct window sequences currently held in the map.
    """

    hits_in_batch: int
    unique_windows: int

    @classmethod
    def from_rust(cls, rust_stats: Any) -> "BatchStats":
        """Build from an opaque Rust ``BatcherStats`` object"""
        return cls(
            hits_in_batch=rust_stats.hits_in_batch,
            unique_windows=rust_stats.unique_windows,
        )

    def __repr__(self) -> str:
        return (
            f"BatchStats(hits={self.hits_in_batch}, "
            f"unique_windows={self.unique_windows})"
        )


@dataclass(frozen=True)
class FeedResult:
    """
    Result returned by :meth:`TargetBatcher.feed_chunk`.

    Attributes
    ----------
    flushed : bool
        ``True`` when the Rust batcher has signalled that the batch is full
        and should be submitted to the pipeline.  The caller **must** call
        ``pipeline.submit(batcher)`` before feeding further chunks.
    stats : BatchStats
        Counter snapshot *after* the chunk was processed.
    """

    flushed: bool
    stats: BatchStats

    @classmethod
    def from_rust(cls, rust_status: Any) -> "FeedResult":
        """Build from an opaque Rust ``FeedStatus`` object"""
        return cls(flushed=rust_status.flushed, stats=rust_status.stats)

    def __repr__(self) -> str:
        return f"FeedResult(flushed={self.flushed}, stats={self.stats})"


# ------------------------------------------------------------------------------
# internal validators
# ------------------------------------------------------------------------------


def _require_native(loggers: CrisprmeLoggers) -> None:
    if RustTargetBatcher is None:
        loggers.errorlog.log_raise_exception(
            "Rust TargetBatcher module not exposed to Python. "
            "Ensure the native extension (_crisprme2_native) is compiled and installed.",
            os.EX_CANTCREAT,
            Crisprme2BatcherError,
        )


def _validate_overlap(size: int, overlap: int, loggers: CrisprmeLoggers) -> None:
    """
    Mirror the Rust constructor guard:
    ``overlap_left`` must be >= ``size - 1`` when ``size > 0``.
    """
    if size > 0 and overlap < size - 1:
        loggers.errorlog.log_raise_exception(
            f"'overlap' ({overlap}) must be >= size - 1 ({size - 1}) "
            "to avoid losing k-mers at chunk boundaries.",
            os.EX_DATAERR,
            Crisprme2BatcherError,
        )


# ------------------------------------------------------------------------------
# public wrapper
# ------------------------------------------------------------------------------


class TargetBatcher:
    """
    Python wrapper around the Rust ``TargetBatcher`` struct.

    The batcher accumulates candidate off-target windows as genome chunks
    are fed in.  When the internal map reaches capacity (measured by either
    total hit count or unique-window count), it signals that the batch
    should be flushed to the alignment pipeline.

    Construction
    ~~~~~~~~~~~~
    Always use `create` - do **not** call ``__init__`` directly.

    ::

        batcher = TargetBatcher.create(
            pam=pam, guide=guide, size=30, upstream=True,
            overlap=29, threads=8, loggers=loggers,
        )

    Parameters
    ----------
    pam : PAM
        Parsed PAM object.  Its ``.pam`` attribute (str) is forwarded to Rust.
    guide : Guide
        Guide RNA object.  Its ``.sequence`` attribute (str) is forwarded.
    size : int
        Window width (guide length + PAM length + any bulge offset).
        Must be strictly positive.
    upstream : bool
        ``True`` if the PAM lies upstream of the protospacer.
    overlap : int
        Left-overlap kept between consecutive chunks (>= ``size - 1``).
    threads : int
        Number of threads for the parallel scanner inside the batcher.
    loggers : CrisprmeLoggers
        Shared logger bundle.
    batch_hits : int
        Flush threshold on total hit count (default: 1 000 000).
    max_unique : int
        Flush threshold on unique-window count (default: 500 000).

    Raises
    ------
    Crisprme2BatcherError
        On invalid arguments or if the native extension is missing.
    """

    def __init__(self, _rust_handle: Any, loggers: CrisprmeLoggers) -> None:
        # DO NOT CALL DIRECTLY -> use TargetBatcher.create()
        self._batcher = _rust_handle
        self._loggers = loggers  # store loggers
        self._total_chunks_fed: int = 0
        self._total_flushes: int = 0

    @classmethod
    def create(
        cls,
        pam: PAM,
        guide: Guide,
        size: int,
        upstream: bool,
        overlap: int,
        threads: int,
        loggers: CrisprmeLoggers,
        batch_hits: int = _DEFAULT_BATCH_HITS,
        max_unique: int = _DEFAULT_MAX_UNIQUE,
    ) -> "TargetBatcher":
        """
        Build and return a new :class:`TargetBatcher`.

        Parameters
        ----------
        pam : PAM
            Parsed PAM object; ``.pam`` str is forwarded to Rust.
        guide : Guide
            Guide RNA object; ``.sequence`` str is forwarded to Rust.
        size : int
            Window extraction width (> 0).
        upstream : bool
            Whether the PAM is upstream of the protospacer.
        overlap : int
            Left overlap between consecutive FASTA chunks (>= size - 1).
        threads : int
            Scanner thread count (> 0).
        loggers : CrisprmeLoggers
            Shared logger bundle.
        batch_hits : int
            Total-hit flush threshold (default 100 000).
        max_unique : int
            Unique-window flush threshold (default 500 000).

        Returns
        -------
        TargetBatcher
            Ready-to-use batcher instance.

        Raises
        ------
        Crisprme2BatcherError
            On invalid arguments or native extension unavailability.
        """
        _require_native(loggers)
        _validate_overlap(size, overlap, loggers)
        pam_seq: str = pam.pam
        guide_seq: str = guide.sequence
        loggers.verboselog.debug(
            f"Constructing TargetBatcher("
            f"pam={pam_seq!r}, guide={guide_seq!r}, size={size}, "
            f"upstream={upstream}, overlap={overlap}, threads={threads}, "
            f"batch_hits={batch_hits}, max_unique={max_unique}"
        )
        try:
            rust_handle = RustTargetBatcher(pam_seq, guide_seq, size, upstream, threads, batch_hits, max_unique, overlap)  # type: ignore
        except Exception as e:
            loggers.errorlog.log_raise_exception(
                f"Rust TargetBatcher construction failed: {e}",
                os.EX_UNAVAILABLE,
                Crisprme2BatcherError,
            )
        loggers.verboselog.debug(f"TargetBatcher created (id={rust_handle.id})")
        return cls(rust_handle, loggers)

    # --------------------------------------------------------------------------
    # properties
    # --------------------------------------------------------------------------

    @property
    def id(self) -> int:
        """Monotonically increasing batcher id assigned by Rust"""
        return self._batcher.id

    @property
    def batcher(self) -> Any:
        """Rust-implemented batcher"""
        return self._batcher

    @property
    def total_chunks_fed(self) -> int:
        """Number of sequence chunks fed since construction"""
        return self._total_chunks_fed

    @property
    def total_flushes(self) -> int:
        """Number of times the batch has been flushed to the pipeline"""
        return self._total_flushes

    # --------------------------------------------------------------------------
    # core API
    # --------------------------------------------------------------------------

    def feed_chunk(
        self, contig_id: int, chunk_start: int, chunk_seq: str, valid_len: int
    ) -> FeedResult:
        """
        Feed a sequence chunk to the batcher.

        The batcher encodes the chunk into IUPAC bitmasks, runs the parallel
        PAM/target scanner, filters positions to the valid core window, and
        accumulates (window -> occurrences) entries into its internal map.

        Parameters
        ----------
        contig_id : int
            Index of the current contig (0-based, fits in u32).
        chunk_start : int
            Absolute genomic start of this chunk in the contig (fits in u32).
            Pass ``0`` for the first chunk.
        chunk_seq : str
            Raw nucleotide string for this chunk (may include overlap region).
        valid_len : int
            Length of the *core* (non-overlap) region within ``chunk_seq``.
            Positions outside this core are discarded by the batcher.

        Returns
        -------
        FeedResult
            ``result.flushed`` is ``True`` when the batch has reached its
            capacity threshold.  The caller **must** submit the batcher to
            the pipeline before feeding the next chunk.

        Raises
        ------
        Crisprme2BatcherError
            On argument validation failure or Rust-side errors (e.g.
            position overflow, encoding failure).
        """
        self._total_chunks_fed += 1  # increase number of actual chunks fed
        try:
            rust_status = self._batcher.feed_chunk(
                contig_id, chunk_start, chunk_seq, valid_len
            )
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"feed_chunk() failed (contig={contig_id}, start={chunk_start}): {e}",
                os.EX_DATAERR,
                Crisprme2BatcherError,
            )
        # feed current sequence chunk to batcher to collect targets
        result = FeedResult.from_rust(rust_status)
        if result.flushed:
            self._loggers.verboselog.debug(
                f"Batcher {self.id}: flush threshold reached - "
                f"{result.stats.hits_in_batch} hits, "
                f"{result.stats.unique_windows} unique windows"
            )
        return result

    def finalize(self) -> BatchStats:
        """
        Flush any remaining state and return final batch statistics.

        Must be called once after all chunks have been fed (including any
        tail chunks after the last automatic flush).  Clears the internal
        map; the batcher should not be used after this call.

        Returns
        -------
        BatchStats
            Counter snapshot of whatever remained in the batch at finalize
            time (before the clear).

        Raises
        ------
        Crisprme2BatcherError
            If the Rust finalize call fails.
        """
        self._loggers.basiclog.info(
            f"Finalizing TargetBatcher (id={self.id}): "
            f"{self._total_chunks_fed} chunks fed, {self._total_flushes} flushes"
        )
        try:
            rust_stats = self._batcher.finalize()
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"finalize() failed: {e}",
                os.EX_IOERR,
                Crisprme2BatcherError,
            )
        return BatchStats.from_rust(rust_stats)

    def stats(self) -> BatchStats:
        """
        Return a live counter snapshot without modifying internal state.

        Returns
        -------
        BatchStats
            Current hit count and unique-window count.

        Raises
        ------
        Crisprme2BatcherError
            If the Rust stats call fails.
        """
        try:
            rust_stats = self._batcher.stats()
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"stats() failed: {e}",
                os.EX_IOERR,
                Crisprme2BatcherError,
            )
        return BatchStats.from_rust(rust_stats)

    # --------------------------------------------------------------------------
    # other helpers
    # --------------------------------------------------------------------------

    def __repr__(self) -> str:
        try:
            s = self.stats()
            stats_str = f"hits={s.hits_in_batch}, unique={s.unique_windows}"
        except Exception:
            stats_str = "stats=unavailable"
        return (
            f"{self.__class__.__name__}(id={self.id}, {stats_str}, "
            f"chunks_fed={self._total_chunks_fed}, flushes={self._total_flushes})"
        )
