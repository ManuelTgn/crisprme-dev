"""
pipeline.py
-----------
Python wrapper for the Rust ``PyPipeline`` struct exposed via PyO3.

The public surface is intentionally minimal: the Rust side owns all
performance-critical state; this wrapper is responsible for:

- Argument validation before any Rust FFI call is made.
- Lifecycle enforcement (the pipeline is a strict context manager).
- Mapping Rust panics / PyO3 errors to typed Python exceptions.
- Providing typed, documented Python entry-points for scanner.py and
  any future callers.

Typical usage
~~~~~~~~~~~~~
::

    from crisprme2.crisprme_core_api import Pipeline, Thresholds

    thresholds = Thresholds(max_mm=4, max_bdna=1, max_brna=1, loggers=loggers)

    with Pipeline.create(chunks=8, thresholds=thresholds, transforms=[my_transform], loggers=loggers) as pipeline:
        pipeline.submit(batcher)   # called once per full batcher flush

Notes
~~~~~
- ``close()`` is **not** part of the public API; it is called exclusively
  inside ``__exit__``.  This prevents double-close bugs and makes the
  usage contract unambiguous.
- Transform validation happens entirely in Python before any Rust call,
  so error messages are clear and do not reference internal Rust symbols.
"""

from __future__ import annotations

from types import TracebackType
from typing import Any, Dict, List, Optional, Type

import os


from ..logger import CrisprmeLoggers
from ..pam import PAM

from .crisprme2_api_error import (
    Crisprme2PipelineConfigError,
    Crisprme2PipelineLifecycleError,
    Crisprme2PipelineSubmitError,
)
from .target_batcher import TargetBatcher
from .thresholds import Thresholds

try:
    from .._crisprme2_native import pipeline as _rust_pipeline_factory
except ImportError:
    # allow the module to be imported in development / doc-build environments
    # where the native extension has not been compiled yet
    _rust_pipeline_factory = None


# ==============================================================================
# internal helpers
# ==============================================================================


def _require_native(loggers: CrisprmeLoggers) -> None:
    """Raise a configuration error if the native extension is unavailable"""
    if _rust_pipeline_factory is None:
        loggers.errorlog.log_raise_exception(
            "Rust pipeline factory not exposed to Python. Ensure the native "
            "extension (_crisprme2_native) is compiled and installed",
            os.EX_CANTCREAT,
            Crisprme2PipelineConfigError,
        )


def _validate_transforms(transforms: List[Any], loggers: CrisprmeLoggers) -> None:
    """
    Validate the *transforms* list.

    Rules:
    - Must be a non-empty list
    - Every element must be callable (``__call__`` present)
    """
    if not isinstance(transforms, list):
        loggers.errorlog.log_raise_exception(
            f"'transforms' must be a list, got {type(transforms).__name__!r}",
            os.EX_DATAERR,
            Crisprme2PipelineConfigError,
        )
    if not transforms:
        loggers.errorlog.log_raise_exception(
            "'transforms' must contain at least one element",
            os.EX_DATAERR,
            Crisprme2PipelineConfigError,
        )
    for i, t in enumerate(transforms):
        if not callable(t):
            loggers.errorlog.log_raise_exception(
                f"Transform at index {i} is not callable "
                f"(type={type(t).__name__!r}). Each transform must "
                "implement __call__",
                os.EX_DATAERR,
                Crisprme2PipelineConfigError,
            )


def _contig_names_in_id_order(
    contig_ids: Dict[str, int], loggers: CrisprmeLoggers
) -> List[str]:
    inverted: Dict[int, str] = {i: n for n, i in contig_ids.items()}
    assert len(inverted) == len(contig_ids)
    expected = set(range(len(contig_ids)))
    if set(inverted) != expected:
        missing = sorted(expected - set(inverted))
        loggers.errorlog.log_raise_exception(
            f"'contig_ids' must be dense 0..{len(contig_ids) - 1}; missing {missing}",
            os.EX_DATAERR,
            Crisprme2PipelineConfigError,
        )
    names = [inverted[i] for i in range(len(inverted))]
    # mirror the Rust guard so the failure names the offending contig
    for name in names:
        bad = next((c for c in ',"\n\r' if c in name), None)
        if bad is not None:
            loggers.errorlog.log_raise_exception(
                f"Contig name {name} contains {bad}, which would break the report structure",
                os.EX_DATAERR,
                Crisprme2PipelineConfigError,
            )
    return names


# ==============================================================================
# public wrapper
# ==============================================================================


class Pipeline:
    """
    Python wrapper around the Rust ``PyPipeline`` processing pipeline.

    The pipeline owns a pool of GPU/CPU memory and a chain of worker
    stages (GPU miner -> resolver -> broadcast -> user transforms -> CSV sink).
    It is **not** reusable: once the context manager exits the underlying
    Rust object is consumed and the instance must be discarded.

    Construction
    ~~~~~~~~~~~~
    Always use the :meth:`create` classmethod — do **not** call ``__init__``
    directly.

    ::

        with Pipeline.create(chunks=8, thresholds=t, transforms=[f]) as p:
            p.submit(batcher)

    Parameters
    ----------
    chunks : int
        Number of memory-pool chunks to pre-allocate.  Each chunk maps to
        ``CHUNK_SIZE`` bytes (defined in the Rust ``columnar`` crate).
        Must be a positive integer.
    thresholds : RustThresholds
        Alignment thresholds (max mismatches, max DNA bulges, max RNA bulges).
        Construct via ``crisprme2._crisprme2_native.Thresholds``.
    transforms : list[callable]
        Ordered list of Python callables that form the transform stage chain.
        Every element must implement ``__call__``; the list must be non-empty.
    loggers : CrisprmeLoggers
        Logger bundle used for structured logging and error propagation.

    Raises
    ------
    Crisprme2PipelineConfigError
        If any argument is invalid, or if the native extension is unavailable.
    Crisprme2PipelineLifecycleError
        If the pipeline is used after it has been closed.
    Crisprme2PipelineSubmitError
        If a batch submission to the Rust pipeline fails.
    """

    def __init__(self, _rust_handle: Any, loggers: CrisprmeLoggers) -> None:
        # _rust_handle is the opaque PyPipeline object returned by the Rust
        # function. Callers should NEVER construct this directly: use
        # Pipeline.create() instead
        self._pipeline = _rust_handle
        self._loggers = loggers
        self._closed: bool = False

    # ==========================================================================
    # construction
    # ==========================================================================

    @classmethod
    def create(
        cls,
        chunks: int,
        thresholds: Thresholds,
        transforms: List[Any],
        pam: PAM,
        upstream: bool,
        outpath: str,
        contig_ids: Dict[str, int],
        loggers: CrisprmeLoggers,
    ) -> "Pipeline":
        """
        Build and return a new :class:`Pipeline` instance.

        This is the only supported constructor.  All arguments are validated
        in Python before the Rust factory is invoked, so errors carry
        descriptive messages without Rust symbol noise.

        Parameters
        ----------
        chunks : int
            Positive integer controlling the size of the pre-allocated memory
            pool (``chunks x CHUNK_SIZE`` bytes).
        thresholds : Thresholds
            ``crisprme2._crisprme2_native.Thresholds`` instance carrying
            ``max_mm``, ``max_bdna``, and ``max_brna`` limits.
        transforms : list[callable]
            Non-empty list of Python callables forming the transform chain.
        pam : PAM
            Parsed PAM; ``.pam`` (str) is forwarded to Rust and rendered into
            the guide column of the CSV report.
        upstream : bool
            ``True``  -> guide column is ``<PAM><aligned-guide>`` (Cas12a TTTV)
            ``False`` -> guide column is ``<aligned-guide><PAM>`` (SpCas9 NGG)
        outpath : str
            Path of the CSV report. Truncated on open.
        contig_ids : dict[str, int]
            Contig ids mapping. Ids must be dense ``0..N-1``; the report
            resolves each ``Occurence``'s contig id through this table.
        loggers : CrisprmeLoggers
            Shared logger bundle.

        Returns
        -------
        Pipeline
            A ready-to-use pipeline, not yet entered as a context manager.

        Raises
        ------
        Crisprme2PipelineConfigError
            On any invalid argument or if the native extension is missing.
        """
        _require_native(loggers)  # ensure native rust api is installed
        _validate_transforms(transforms, loggers)  # ensure transforms are callable
        contigs = _contig_names_in_id_order(contig_ids, loggers)
        loggers.verboselog.debug(
            f"Constructing Pipeline (chunks={chunks}). num_transforms={len(transforms)} "
            f"pam={pam.pam!r}, upstream={upstream}, output={outpath!r}, contigs={len(contigs)}))"
        )
        try:
            rust_handle = _rust_pipeline_factory(chunks, thresholds.rust_handle, transforms, pam.pam, upstream, outpath, contigs)  # type: ignore
        except Exception as e:
            loggers.errorlog.log_raise_exception(
                f"Rust pipeline initialization failed: {e}",
                os.EX_UNAVAILABLE,
                Crisprme2PipelineConfigError,
            )
        loggers.basiclog.info(
            f"Pipeline created (chunks={chunks}, transforms={len(transforms)} "
            f"pam={pam.pam!r}, upstream={upstream}, contigs={len(contigs)}))"
        )
        return cls(rust_handle, loggers)

    # ==========================================================================
    # lifecycle helpers
    # ==========================================================================

    def _assert_open(self) -> None:
        """
        Guard: raise :exec:`Crisprme2PipelineLifeCycleError` if the pipeline has
        already been closed.
        """
        if self._closed:
            self._loggers.errorlog.log_raise_exception(
                "Pipeline has already been closed. A pipeline instance cannot "
                "be reused after its context manager exits",
                os.EX_UNAVAILABLE,
                Crisprme2PipelineLifecycleError,
            )

    def _close(self) -> None:
        """
        Internal shutdown: signal EOF to the Rust pipeline and join all
        worker threads.  Idempotent — safe to call more than once but the
        second call is a no-op.
        """
        if self._closed:
            return
        self._closed = True
        self._loggers.verboselog.debug("Closing pipeline end joining worker threads")
        try:
            self._pipeline.close()
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Error while closing pipeline: {e}",
                os.EX_UNAVAILABLE,
                Crisprme2PipelineLifecycleError,
            )
        self._loggers.basiclog.info("Pipeline closed")

    # ==========================================================================
    # Context manager protocol
    # ==========================================================================

    def __enter__(self) -> "Pipeline":
        """
        Enter the pipeline context.

        Raises
        ------
        Crisprme2PipelineLifecycleError
            If the pipeline has already been closed (re-entry guard).
        """
        self._assert_open()
        self._loggers.verboselog.debug("Entering Pipeline context")
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> bool:
        """
        Exit the pipeline context, always closing the underlying Rust object.

        The pipeline is closed regardless of whether an exception occurred.
        Any exception propagates normally (return value is ``False``).
        If :meth:`_close` itself raises, that secondary exception is logged
        and re-raised, masking the original only when both fail simultaneously.
        """
        self._loggers.verboselog.debug(
            f"Exiting Pipeline context (exc_type={exc_type})"
        )
        self._close()
        return False  # do not suppress exceptions raised inside the with-block

    # ==========================================================================
    # Public API
    # ==========================================================================

    def submit(self, batcher: TargetBatcher) -> None:
        """
        Submit the contents of a :class:`~crisprme2.crisprme_core_api.TargetBatcher`
        to the pipeline for alignment and scoring.

        The batcher is drained (``flush_to_batch`` is called on the Rust
        side) and its windows + occurrences are transferred into pipeline
        memory frames.  The GIL is released while the data is sent across
        the channel so worker threads can make progress concurrently.

        Parameters
        ----------
        batcher : TargetBatcher
            A populated batcher instance.  Must be a ``RustTargetBatcher``
            (i.e. ``crisprme2._crisprme2_native.TargetBatcher``) wrapped or
            unwrapped.

        Raises
        ------
        Crisprme2PipelineLifecycleError
            If the pipeline has already been closed.
        Crisprme2PipelineSubmitError
            If the batcher type is wrong, or if the Rust ``submit()`` call
            raises (e.g. channel disconnected, sequence too long).
        """
        self._assert_open()
        self._loggers.verboselog.debug("Submitting TargetBatcher batch to pipeline")
        try:
            self._pipeline.submit(batcher.batcher)
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Pipeline batch submission failed: {e}",
                os.EX_IOERR,
                Crisprme2PipelineSubmitError,
            )
        self._loggers.verboselog.debug("Batch submitted successfully")

    # ==========================================================================
    # other helpers
    # ==========================================================================

    @property
    def is_closed(self) -> bool:
        """``True`` once the context manager has exited"""
        return self._closed

    def __repr__(self) -> str:
        status = "closed" if self._closed else "open"
        return f"Pipeline(status={status!r})"
