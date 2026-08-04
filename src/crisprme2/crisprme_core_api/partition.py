"""
partition.py
------------
Python wrapper for the Rust report partitioner exposed via PyO3
(``_crisprme2_native.partition_report``).

Splits the single intermediate off-target report into two reports:

- **primary**     - one representative alignment per cluster
- **alternative** - every other alignment in the cluster

A *cluster* groups hits mapped within ``n`` bases, keyed on
``(contig, strand, start)`` (strand-aware; only ``n`` is user-tunable). The
primary within each cluster is chosen by a configurable ranking (``criteria``)
whose fixed ``(strand, target_aligned)`` tiebreak is applied last on the Rust
side. Output lines are byte-identical to the intermediate, so all score and
annotation columns are preserved verbatim.

Responsibilities
~~~~~~~~~~~~~~~~
- Validate ``n`` and the three report paths before crossing the FFI boundary.
- Normalise the optional ``criteria`` spec into ``(str, str)`` pairs.
- Map every native failure to the typed :class:`Crisprme2PartitionError` so
  callers get a predictable exception type.

What this does NOT do
~~~~~~~~~~~~~~~~~~~~~
- It does not choose the report paths - the composition root
  (``search.py``) derives them from the intermediate report name.
- It does not run the search - it is a finalization step invoked only
  after the pipeline context has closed and the intermediate is flushed.

Typical usage (from search.py, after the pipeline context exits)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    partition_offtargets(
        intermediate=outpath,
        primary_out=primary_path,
        alternative_out=alternative_path,
        n=3,
        loggers=loggers,
        criteria=None,  # default: edit distance, dna/rna bulges, mismatches (asc)
    )
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import os

from .crisprme2_api_error import Crisprme2PartitionError
from ..logger import CrisprmeLoggers

try:  # import rust API function
    from .._crisprme2_native import partition_report as _native_partition_report
except ImportError:
    # fallback for development/testing before the native extension is built
    _native_partition_report = None


def _require_native(loggers: CrisprmeLoggers) -> None:
    """Halt the run if the native extension has not been compiled."""
    if _native_partition_report is None:
        loggers.errorlog.log_raise_exception(
            "Rust partition_report not exposed to Python. Ensure the native "
            "extension (_crisprme2_native) is compiled and installed.",
            os.EX_CANTCREAT,
            Crisprme2PartitionError,
        )


def partition_offtargets(
    intermediate: str,
    primary_out: str,
    alternative_out: str,
    criteria: List[str],
    n: int,
    loggers: CrisprmeLoggers,
) -> Tuple[int, int, int]:
    """
    Split ``intermediate`` into ``primary_out`` and ``alternative_out``.

    Parameters
    ----------
    intermediate : str
        Path of the combined report written by the pipeline. Must exist and be
        fully flushed (call only after the ``Pipeline`` context has exited).
    primary_out, alternative_out : str
        Destination paths. Truncated on open; must differ from each other and
        from ``intermediate``.
    n : int
        Clustering window in bases (non-negative). The only user-tunable
        clustering parameter.
    loggers : CrisprmeLoggers
        Shared logger bundle.
    criteria : sequence of str, optional
        Ordered primary-selection fields, each drawn from :data:`CRITERIA`. The
        order sets ranking priority; direction is inherent per field
        (edit-family ascending, scores descending). ``None`` uses the default
        (``edit-dist``, ``bdna``, ``brna``, ``mm``).

    Returns
    -------
    (clusters, primary, alternative) : tuple[int, int, int]
        Counts, for logging/verification.

    Raises
    ------
    Crisprme2PartitionError
        On invalid arguments, a missing native extension, or any native
        parsing / I/O failure (the underlying message is preserved).
    """
    _require_native(loggers)
    assert _native_partition_report is not None
    if not isinstance(n, int) or n < 0:
        loggers.errorlog.log_raise_exception(
            f"Clustering window n must be a non-negative integer, got {n!r}",
            os.EX_USAGE,
            Crisprme2PartitionError,
        )
    # Outputs are truncate-on-create: refuse to clobber the intermediate or
    # collapse the two reports onto one path.
    if len({Path(intermediate), Path(primary_out), Path(alternative_out)}) != 3:
        loggers.errorlog.log_raise_exception(
            "intermediate, primary, and alternative report paths must be "
            f"distinct (got {intermediate!r}, {primary_out!r}, {alternative_out!r})",
            os.EX_USAGE,
            Crisprme2PartitionError,
        )
    try:
        clusters, primary, alternative = _native_partition_report(
            intermediate, primary_out, alternative_out, n, criteria
        )
    except Exception as e:  # native ValueError / OSError -> typed API error
        loggers.errorlog.log_raise_exception(
            f"Report partitioning failed: {e}",
            os.EX_DATAERR,
            Crisprme2PartitionError,
        )
    loggers.basiclog.info(
        f"Partitioned report: {clusters} clusters -> "
        f"{primary} primary, {alternative} alternative"
    )
    return clusters, primary, alternative
