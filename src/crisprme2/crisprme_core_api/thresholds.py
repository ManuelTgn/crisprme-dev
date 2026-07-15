"""
thresholds.py
-------------
Python wrapper around the Rust ``Thresholds`` value type exposed via PyO3.

Because the Rust object is opaque after construction (no readable fields),
this wrapper stores the three threshold values in Python and acts as the
single source of truth for inspection, hashing, equality, and repr.

The validated Python values are forwarded to ``RustThresholds`` exactly
once at construction time; after that the Rust handle is carried opaquely
and handed to :class:`~crisprme2.crisprme_core_api.Pipeline` when needed.

Typical usage
~~~~~~~~~~~~~
::

    from crisprme2.crisprme_core_api import Thresholds

    # Exact-match only (zero is valid)
    exact = Thresholds(max_mm=0, max_bdna=0, max_brna=0)

    # Up to 4 mismatches, 1 DNA bulge, 1 RNA bulge
    t = Thresholds(max_mm=4, max_bdna=1, max_brna=1)

    with Pipeline.create(chunks=8, thresholds=t, transforms=[f], loggers=log) as p:
        p.submit(batcher)

Notes
~~~~~
- All three parameters must be non-negative integers (``>= 0``).
- ``max_mm=0, max_bdna=0, max_brna=0`` activates exact-match mode.
- The wrapper is **immutable**: there are no setters. Build a new instance
  if different thresholds are needed.
- The wrapper is **hashable** and implements ``__eq__`` so it can be used
  as a dict key or in sets.
"""

from __future__ import annotations

from typing import Any

import os

from .crisprme2_api_error import Crisprme2PipelineConfigError
from ..logger import CrisprmeLoggers

try:
    from .._crisprme2_native import Thresholds as RustThresholds
except ImportError:
    RustThresholds = None


# ==============================================================================
# internal helpers
# ==============================================================================


def _require_native(loggers: CrisprmeLoggers) -> None:
    """Raise if the native extension has not been compiled"""
    if RustThresholds is None:
        loggers.errorlog.log_raise_exception(
            "Rust Thresholds type not exposed to Python. Ensure the native "
            "extension (_crisprme2_native) is compiled and installed.",
            os.EX_CANTCREAT,
            Crisprme2PipelineConfigError,
        )


# ==============================================================================
# public python wrapper
# ==============================================================================


class Thresholds:
    """
    Immutable alignment-threshold value type.

    Wraps the Rust ``Thresholds`` struct and stores the three constituent
    values in Python for inspection, equality testing, and hashing.  The
    underlying Rust handle is opaque and is accessed only by pipeline
    internals via :attr:`rust_handle`.

    Parameters
    ----------
    max_mm : int
        Maximum number of mismatches allowed per alignment (``>= 0``).
        Pass ``0`` for exact-match-only mode.
    max_bdna : int
        Maximum number of DNA-strand bulges allowed (``>= 0``).
    max_brna : int
        Maximum number of RNA-strand bulges allowed (``>= 0``).
    loggers : CrisprmeLoggers
        Shared logger bundle used for validation errors and info messages.

    Raises
    ------
    Crisprme2PipelineConfigError
        If any value is not a non-negative integer, or if the native
        extension is unavailable.

    Examples
    --------
    Exact-match mode::

        t = Thresholds(max_mm=0, max_bdna=0, max_brna=0, loggers=loggers)

    Tolerant search::

        t = Thresholds(max_mm=4, max_bdna=1, max_brna=1, loggers=loggers)
    """

    __slots__ = ("_max_mm", "_max_bdna", "_max_brna", "_loggers", "_rust_handle")

    def __init__(
        self, max_mm: int, max_bdna: int, max_brna: int, loggers: CrisprmeLoggers
    ) -> None:
        _require_native(loggers)  # check that rust extension is available
        self._loggers = loggers  # set loggers
        # set Thresholds class attributes
        self._max_mm = max_mm
        self._max_bdna = max_bdna
        self._max_brna = max_brna
        loggers.verboselog.debug(f"Constructing object {repr(self)}")
        try:  # initialize rust Thresholds object
            self._rust_handle: Any = RustThresholds(max_brna, max_bdna, max_mm)  # type: ignore[assignment]
        except Exception as e:
            loggers.errorlog.log_raise_exception(
                f"Rust Thresholds construction failed: {e}",
                os.EX_DATAERR,
                Crisprme2PipelineConfigError,
            )

    # ==========================================================================
    # read-only fields accessors
    # (values are stored in Python; the rust handle is opaque)
    # ==========================================================================

    @property
    def mm(self) -> int:
        """Maximum mismatches allowed (``>= 0``)"""
        return self._max_mm

    @property
    def bdna(self) -> int:
        """Maximum DNA bulges allowed (``>= 0``)"""
        return self._max_bdna

    @property
    def brna(self) -> int:
        """Maximum RNA bulges allowed (``>= 0``)"""
        return self._max_brna

    @property
    def is_exact_match(self) -> bool:
        """``True`` when all thresholds are zero (exact-match mode)"""
        return self._max_mm == 0 and self._max_bdna == 0 and self._max_brna == 0

    # ==========================================================================
    # internal: raw rust handle for pipeline internals
    # ==========================================================================

    @property
    def rust_handle(self) -> Any:
        """
        The opaque Rust ``Thresholds`` object.

        This property is intended for use by pipeline internals
        (e.g. :class:`~crisprme2.crisprme_core_api.Pipeline`) that need to
        pass the raw handle across the FFI boundary.  DO NOT rely on its
        type or attributes in application code
        """
        return self._rust_handle

    # ==========================================================================
    # value-type protocol
    # ==========================================================================

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Thresholds):
            return NotImplemented
        return (
            self._max_mm == other._max_mm
            and self._max_bdna == other._max_bdna
            and self._max_brna == other._max_brna
        )

    def __hash__(self) -> int:
        return hash((self._max_mm, self._max_bdna, self._max_brna))

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(max_mm={self._max_mm}, "
            f"max_bdna={self._max_bdna}, max_brna={self._max_brna})"
        )

    def __str__(self) -> str:
        return f"MM={self._max_mm}, BDNA={self._max_bdna}, BRNA={self._max_brna}"
