"""
merge.py
--------
Python wrapper for the Rust assembly merge exposed via PyO3
(``_crisprme2_native.merge_assemblies``).

Merges every sample's per-haplotype *lifted* reports into one final report in
reference coordinates. Within a sample, haplotype hits at the same
sequence-aware locus collapse (copy-wise OR -> homozygous); across samples,
shared sites collapse and their carriers union. The ``chromosome``/``start``/
``strand`` columns carry reference coordinates for mapped rows and native
assembly coordinates for ``unmapped`` (assembly-specific) rows; two columns
trail: ``samples`` (e.g. ``HG1:1|1,HG2:1|0``) then ``map_status``.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

import os

from ..logger import CrisprmeLoggers

try:
    from .._crisprme2_native import merge_assemblies as _native_merge_assemblies
except ImportError:
    _native_merge_assemblies = None


def _require_native(loggers: CrisprmeLoggers) -> None:
    if _native_merge_assemblies is None:
        loggers.errorlog.log_raise_exception(
            "Rust merge_assemblies not exposed to Python. Ensure the native "
            "extension (_crisprme2_native) is compiled and installed.",
            os.EX_CANTCREAT,
            ValueError,
        )


def merge_assemblies(
    sample_names: Sequence[str],
    hap_layout: Sequence[Sequence[int]],
    reports: Sequence[Tuple[str, int, int]],
    header: str,
    out_path: str,
    merge_bp: int,
    criteria: List[str],
    loggers: CrisprmeLoggers,
) -> int:
    """
    Merge lifted haplotype reports into ``out_path``; return rows written.

    Parameters
    ----------
    sample_names : sequence of str
        Sample names in u32-id order (from the run's ``SampleTable``).
    hap_layout : sequence of sequence of int
        Each sample's sorted PanSN hap_ids; length == that sample's ploidy.
    reports : sequence of (path, sample_index, hap_id)
        One lifted report per haplotype.
    header : str
        The reference ``REPORT_HEADER`` to reuse verbatim.
    out_path : str
        Final merged report path (truncated on open).
    merge_bp : int
        Single-linkage clustering window (the search's ``--merge``).
    criteria : list of str
        Primary-selection field tokens; empty uses the default order.
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Raises
    ------
    Crisprme2ReportError
        On a missing native extension or any native parse / I/O failure.
    """
    _require_native(loggers)
    assert _native_merge_assemblies is not None
    if not isinstance(merge_bp, int) or merge_bp < 0:
        loggers.errorlog.log_raise_exception(
            f"merge_bp must be a non-negative integer, got {merge_bp!r}",
            os.EX_USAGE,
            ValueError,
        )
    try:
        return _native_merge_assemblies(
            list(sample_names),
            [list(h) for h in hap_layout],
            [tuple(r) for r in reports],
            header,
            out_path,
            merge_bp,
            list(criteria),
        )
    except Exception as e:
        loggers.errorlog.log_raise_exception(
            f"Assembly merge failed: {e}",
            os.EX_DATAERR,
            ValueError,
        )
