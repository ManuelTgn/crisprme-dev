"""
feature_registry.py
--------------------
Python-side handle for one annotation BED file's term <-> bit mapping.

Role
~~~~
Functional annotation encodes, per off-target, the set of BED column-4 *terms*
it overlaps as a 32-bit mask — one bit per distinct term (mirrored by the Rust
``FeatureRegistry``). This class is the **encode** half of that scheme on the
Python side:

* It builds the native ``PyRegistry`` from a BED path and pulls the whole
  ``term -> bitmask`` mapping across the FFI boundary **once**, caching it as a
  plain ``dict[str, int]``.
* The annotation transform then turns each pysam-fetched term into bits with
  :meth:`accumulate` / :meth:`mask_for` and ORs them into the alignment's
  ``u32`` annotation slot — with no per-term calls back into Rust.

The **decode** half — expanding an accumulated ``u32`` back into a comma-joined
list of term names for the report — lives Rust-side in the TSV sink. This class
never renders names.

Notes
~~~~~
* One instance per BED. The cached mapping is small (<= 32 entries) and
  immutable; :attr:`masks` hands out a read-only view of it.
* A term fetched at annotate time that is absent from the registry contributes
  no bits, so a stray label degrades to "unannotated" rather than raising in the
  middle of a batch.
* All initialization failures route through ``loggers.errorlog`` (which logs the
  traceback and halts with the given exit code) rather than surfacing a raw
  exception to the caller.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from types import MappingProxyType

import os


from ..logger import CrisprmeLoggers

from .crisprme2_api_error import Crisprme2AnnotationError

try:  # native extension, built via maturin
    from .._crisprme2_native import PyRegistry as _PyRegistry
except ImportError:  # extension not built yet (dev / test)
    _PyRegistry = None


class FeatureRegistry:
    """Cached ``term -> bitmask`` mapping for a single annotation BED.

    Parameters
    ----------
    path : str
        Path to the annotation BED (plain or bgzip-compressed; compression is
        detected by magic bytes, so the extension is irrelevant).
    loggers : CrisprmeLoggers
        Shared logger bundle; every error is raised through it.
    """

    __slots__ = ("_loggers", "_path", "_registry", "_masks")

    def __init__(self, path: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers
        self._path = path
        if _PyRegistry is None:
            loggers.errorlog.log_raise_exception(
                "Native PyRegistry is unavailable; build the extension with "
                "`maturin develop` before running annotation.",
                os.EX_CANTCREAT,
                Crisprme2AnnotationError,
            )
        loggers.basiclog.info(f"Loading annotation features from {path}")
        try:
            self._registry = _PyRegistry(path)
            # Cross the FFI boundary exactly once for the full mapping; every
            # later lookup is then a pure-Python dict hit.
            self._masks = self._registry.masks()
        except Exception as e:  # native ValueError / IOError, malformed BED, ...
            loggers.errorlog.log_raise_exception(
                f"Feature registry initialization failed on {path}: {e}",
                os.EX_IOERR,
                Crisprme2AnnotationError,
            )
        loggers.verboselog.info(
            f"FeatureRegistry ready: {len(self._masks)} terms from {path}"
        )

    # -- introspection ---------------------------------------------------------

    @property
    def path(self) -> str:
        """Path of the BED this registry was built from."""
        return self._path

    @property
    def num_features(self) -> int:
        """Number of distinct annotation terms (bits in use); ``0..=32``."""
        return len(self._masks)

    @property
    def masks(self) -> Mapping[str, int]:
        """Read-only view of the ``term -> bitmask`` mapping.

        Each value is a pre-shifted mask (``1 << bit``). The view shares the
        cached dict, so it is cheap to obtain; it must not be mutated (and can't
        be, being a ``MappingProxyType``).
        """
        return MappingProxyType(self._masks)

    # -- encoding --------------------------------------------------------------

    def mask_for(self, term: str) -> int:
        """Bitmask for a single term, or ``0`` if the term is not registered."""
        return self._masks.get(term, 0)

    def accumulate(self, terms: Iterable[str]) -> int:
        """OR the bitmasks of ``terms`` into one ``u32`` annotation value.

        This is the value the transform writes into an alignment's annotation
        slot. Unregistered terms contribute nothing, and the result is always
        <= 32 bits wide by construction.
        """
        acc = 0
        get = self._masks.get  # bind once for the hot loop
        for term in terms:
            acc |= get(term, 0)
        return acc

    # -- dunder ----------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._masks)

    def __repr__(self) -> str:
        return f"<FeatureRegistry path={self._path!r} terms={len(self._masks)}>"
