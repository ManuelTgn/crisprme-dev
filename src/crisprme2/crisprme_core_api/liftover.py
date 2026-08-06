"""
liftover.py
-----------
Python wrapper for the Rust report liftover exposed via PyO3
(``_crisprme2_native.lift_report``).

Lifts one per-haplotype intermediate off-target report from assembly
coordinates to reference (e.g. GRCh38) coordinates, using that haplotype's own
chain file. Coordinates are added *additively*: the four columns
``hg38_chr, hg38_pos, hg38_flipped, map_status`` are appended and native
columns are preserved verbatim, so every row is kept — including
``unmapped`` (assembly-specific) rows, which carry ``NA`` coordinates.

Each row is lifted at its leftmost-forward ``start`` anchor (the same anchor
the hg38 clustering key uses). A row is ``ambiguous`` when a lower-score chain
covering the same position maps it more than ``tolerance`` bp away from the
primary chain's mapping.

Responsibilities
~~~~~~~~~~~~~~~~
- Validate ``tolerance`` and the paths before crossing the FFI boundary.
- Map every native failure to the typed :class:`Crisprme2LiftoverError`.

What this does NOT do
~~~~~~~~~~~~~~~~~~~~~
- It does not choose paths or merge haplotypes — the finalization step
  (``search.py``) owns the ``.assemblies`` layout and the downstream merge.
"""

from __future__ import annotations

from typing import Tuple

import os

from .crisprme2_api_error import Crisprme2LiftoverError
from ..logger import CrisprmeLoggers

try:  # import rust API function
    from .._crisprme2_native import lift_report as _native_lift_report
except ImportError:
    # fallback for development/testing before the native extension is built
    _native_lift_report = None


def _require_native(loggers: CrisprmeLoggers) -> None:
    """Halt the run if the native extension has not been compiled."""
    if _native_lift_report is None:
        loggers.errorlog.log_raise_exception(
            "Rust lift_report not exposed to Python. Ensure the native "
            "extension (_crisprme2_native) is compiled and installed.",
            os.EX_CANTCREAT,
            Crisprme2LiftoverError,
        )


def lift_offtargets(
    report: str,
    chain: str,
    out_path: str,
    tolerance: int,
    loggers: CrisprmeLoggers,
    contig_col: str = "chromosome",
    start_col: str = "start",
) -> Tuple[int, int, int]:
    """
    Lift ``report`` to reference coordinates, writing ``out_path``.

    Parameters
    ----------
    report : str
        Per-haplotype intermediate report (native coordinates). Must exist.
    chain : str
        Chain file lifting this haplotype to the reference. Must exist.
    out_path : str
        Destination for the hg38-augmented report. Truncated on open; must
        differ from ``report``.
    tolerance : int
        Ambiguity threshold in bases (non-negative).
    loggers : CrisprmeLoggers
        Shared logger bundle.
    contig_col, start_col : str
        Column names resolved by the native side (order-independent).

    Returns
    -------
    (mapped, ambiguous, unmapped) : tuple[int, int, int]
        Row counts, for logging/verification.

    Raises
    ------
    Crisprme2LiftoverError
        On invalid arguments, a missing native extension, or any native
        parsing / I/O failure (the underlying message is preserved).
    """
    _require_native(loggers)
    assert _native_lift_report is not None
    if not isinstance(tolerance, int) or tolerance < 0:
        loggers.errorlog.log_raise_exception(
            f"Ambiguity tolerance must be a non-negative integer, got {tolerance!r}",
            os.EX_USAGE,
            Crisprme2LiftoverError,
        )
    if os.path.abspath(out_path) == os.path.abspath(report):
        loggers.errorlog.log_raise_exception(
            "Liftover output path must differ from the input report",
            os.EX_USAGE,
            Crisprme2LiftoverError,
        )
    try:
        mapped, ambiguous, unmapped = _native_lift_report(
            report, chain, out_path, tolerance, contig_col, start_col
        )
    except Exception as e:  # native raises builtin ValueError / IOError
        loggers.errorlog.log_raise_exception(
            f"Liftover failed for {report}: {e}",
            os.EX_DATAERR,
            Crisprme2LiftoverError,
        )
    return mapped, ambiguous, unmapped