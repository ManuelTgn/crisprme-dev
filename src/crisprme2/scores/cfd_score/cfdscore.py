"""
cfdscore.py
-----------
CFD (Cutting Frequency Determination) score model for CRISPR off-target
assessment (Doench et al., 2016, Nature Biotechnology).

Module structure
~~~~~~~~~~~~~~~~
- :func:`load_models`     — load mismatch and PAM score tables from disk once.
- :class:`CfdScorer`      — callable transform that scores an entire
                            :class:`~crisprme2.crisprme_core_api.AlignmentBatch`
                            in-place, writing results to score slot 0.

Scoring logic
~~~~~~~~~~~~~
The CFD score for a guide/off-target pair is::

    cfd = PAM_weight(concrete_pam) x prod_i mismatch_weight(rX:dY, position_i)

where the product runs over the protospacer positions where guide != target
(gaps from bulge alignments are skipped), positions are 1-indexed from the
**PAM-distal** end, and ``concrete_pam`` is the actual PAM the off-target
matched (see "PAM weighting" below).

Sequence encoding (IMPORTANT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``rguide`` and ``rseq`` reach this transform as **ASCII** bytes written by the
Rust ``Resolver`` — ``A``/``C``/``G``/``T`` for bases, ``-`` for a bulge gap,
and ``0x00`` padding past the resolved length. They are *not* IUPAC bitmasks.
The kernel decodes each byte to an index 0..3 (A/C/G/T; U folds onto T); every
other byte (``-``, ``0x00``, ambiguous codes) decodes to ``-1`` and is skipped
with weight ``1.0``. ``rseq`` carries the protospacer only — the PAM is never
in it (it lives in ``pam_id``).

PAM weighting
~~~~~~~~~~~~~
The CFD PAM weight depends on the two 3'-most PAM bases and therefore varies
per off-target whenever the run PAM is degenerate in those positions (e.g.
``NRG`` -> the middle base is ``A`` or ``G`` per hit, with different weights).
So the weight is resolved **per row** from ``batch.pam_id``: at construction we
enumerate the run PAM's concrete variants in ``pam_id`` order (via the native
``pam_variants_ascii``, the same enumeration ``Occurence::pam()`` indexes) and
precompute a weight per id. At call time this is a single vectorised gather.
The ``PAM_ID_NONE`` sentinel (and any out-of-range id) is weighted
``_PAM_NONE_WEIGHT``.

Vectorised implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``__call__`` path avoids Python loops over rows. ``rguide``/``rseq`` are
read as flat ``uint8`` views, reshaped to ``(N, seq_len)``, and processed
column-by-column so NumPy does the heavy lifting; the PAM weight is a single
fancy-index gather.

Score slot convention
~~~~~~~~~~~~~~~~~~~~~
CFD scores are written to **score slot 0** (``batch.score(0)``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

import os
import pickle

from ...crisprme_core_api import AlignmentBatch
from ...dna_alphabet import RNA, RC, dna2rna_nt
from ...logger import CrisprmeLoggers
from ...protocol import Transformer

from ..crisprme2_scores_error import Crisprme2CfdScoreError


try:  # Concrete-PAM-variant enumeration, in the exact order `Occurence::pam()`
    from ..._crisprme2_native import pam_variants_ascii
except ImportError:
    pam_variants_ascii = None


# ==============================================================================
# paths
# ==============================================================================

_MODELS_DIR: Path = Path(__file__).parent / "models"
_MM_SCORES_FILE: Path = _MODELS_DIR / "mismatch_score.pkl"
_PAM_SCORES_FILE: Path = _MODELS_DIR / "pam_scores.pkl"


# ==============================================================================
# constants
# ==============================================================================

#: score slot index for CFD within PyAlignmentBatch
CFD_SCORE_SLOT: int = 0

#: number of positions scored by the CFD model (PAM-distal position 1..20)
_CFD_POSITIONS: int = 20

#: sentinel written by Rust for "no concrete PAM" (matches PAM_ID_NONE = u16::MAX)
PAM_ID_NONE: int = 0xFFFF

#: weight applied to rows whose PAM could not be resolved to a concrete variant.
#: 0.0 => an off-target without a confirmed cutting-competent PAM predicts no cut.
_PAM_NONE_WEIGHT: float = 0.0


# ==============================================================================
# model types
# ==============================================================================

MismatchScores = Dict[str, float]
PamScores = Dict[str, float]


# ==============================================================================
# ASCII -> base-index decode table (module-level, built once)
# ==============================================================================

# Maps an ASCII byte to a base index 0..3 (A/C/G/T; U folds onto T to match the
# RNA-alphabet convention used by `_build_mm_lookup`). Everything else -- '-'
# (bulge gap), 0x00 (padding), and any ambiguous IUPAC letter -- stays -1 and is
# skipped during scoring (weight 1.0). Both cases are mapped so lowercase
# soft-masked bases (if ever emitted) still score.
_ASCII_TO_IDX: np.ndarray = np.full(256, -1, dtype=np.int8)
for _base, _idx in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("U", 3)):
    _ASCII_TO_IDX[ord(_base)] = _idx
    _ASCII_TO_IDX[ord(_base.lower())] = _idx


# ==============================================================================
# model loading
# ==============================================================================


def load_models(loggers: CrisprmeLoggers) -> Tuple[MismatchScores, PamScores]:
    """
    Load CFD mismatch and PAM score tables from the bundled pickle files.

    The models are the original tables from Doench et al. (2016). Both files
    are shipped with the package.

    Parameters
    ----------
    loggers : CrisprmeLoggers
        Shared logger bundle used for error reporting.

    Returns
    -------
    tuple[MismatchScores, PamScores]
        ``(mm_scores, pam_scores)`` where each is a ``dict[str, float]``.

    Raises
    ------
    Crisprme2CfdScoreError
        If either pickle file cannot be found or unpickled.
    """
    try:  # load mismatch score model
        with open(_MM_SCORES_FILE, mode="rb") as fh:
            mm_scores: MismatchScores = pickle.load(fh)
    except Exception as e:
        loggers.errorlog.log_raise_exception(
            f"Failed to load CFD mismatch scores from {_MM_SCORES_FILE}: {e}",
            os.EX_IOERR,
            Crisprme2CfdScoreError,
        )
    try:  # load pam score model
        with open(_PAM_SCORES_FILE, mode="rb") as fh:
            pam_scores: PamScores = pickle.load(fh)
    except Exception as e:
        loggers.errorlog.log_raise_exception(
            f"Failed to load CFD PAM scores from {_PAM_SCORES_FILE}: {e}",
            os.EX_IOERR,
            Crisprme2CfdScoreError,
        )
    return mm_scores, pam_scores


# ==============================================================================
# pre-computed lookup tables
# ==============================================================================


def _build_mm_lookup(mm_scores: MismatchScores, loggers: CrisprmeLoggers) -> np.ndarray:
    """
    Build a dense ``(20, 4, 4)`` float32 lookup table from the sparse mismatch
    score dict.

    Axes: ``[position-1 (0..19), guide_base (A/C/G/T), target_base (A/C/G/T)]``,
    with the base index following the RNA-alphabet order (A=0, C=1, G=2, U/T=3),
    which matches the kernel's ASCII decode.

    Missing entries default to ``1.0`` (perfect-match weight), which is safe:
    guide == target positions are never keyed in the model anyway.

    The CFD key format is ``"rX:dY,pos"`` where ``X`` is the RNA (guide) base and
    ``Y`` is the DNA (target) base **reverse-complemented** — so we RC ``Y`` back
    to recover the actual target base before indexing.

    Parameters
    ----------
    mm_scores : MismatchScores
        Raw dict loaded from ``mismatch_score.pkl``.
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Returns
    -------
    np.ndarray
        Shape ``(20, 4, 4)``, dtype ``float32``.
    """
    base_idx = {nt: i for i, nt in enumerate(RNA)}  # A=0, C=1, G=2, U=3
    table = np.ones((20, 4, 4), dtype=np.float32)
    for key, val in mm_scores.items():
        try:  # key format: "rX:dY,pos" e.g. "rU:dT,12"
            pair, pos_str = key.split(",")
            pos = int(pos_str) - 1  # convert to 0-based (position 1 -> index 0)
            rna_nt = pair[1]  # 'r' + guide base
            dna_nt = pair[4]  # 'd' + RC(target base)
            # undo the RC stored in the key to recover the actual target base
            target_rna_nt = dna2rna_nt(RC.get(dna_nt, dna_nt))
            g_idx, t_idx = base_idx.get(rna_nt), base_idx.get(target_rna_nt)
            if g_idx is None or t_idx is None or not (0 <= pos < 20):
                continue
            table[pos, g_idx, t_idx] = float(val)
        except (ValueError, IndexError) as e:
            loggers.verboselog.debug(
                f"CFD mismatch table: skipping malformed key {key!r}: {e}"
            )
    return table


def _build_pam_weight_by_id(
    variants: List[str],
    pam_scores: PamScores,
    loggers: CrisprmeLoggers,
) -> np.ndarray:
    """
    Precompute a per-``pam_id`` PAM weight vector.

    ``variants[i]`` is the concrete PAM whose ``pam_id == i`` (as produced by the
    native ``pam_variants_ascii``, matching the id packed into each
    ``Occurence``). The CFD PAM weight keys on the two 3'-most PAM bases, so the
    weight for id ``i`` is ``pam_scores[variants[i][-2:]]``.

    Parameters
    ----------
    variants : list[str]
        Concrete PAM strings in ``pam_id`` order.
    pam_scores : PamScores
        Raw dict loaded from ``pam_scores.pkl`` (16 two-character keys).
    loggers : CrisprmeLoggers
        Shared logger bundle.

    Returns
    -------
    np.ndarray
        Shape ``(len(variants),)``, dtype ``float32``. Indexed directly by
        ``pam_id``.

    Raises
    ------
    Crisprme2CfdScoreError
        If a variant's 2-base key is absent from the PAM score table.
    """
    pam_scores_u = {k.upper(): float(v) for k, v in pam_scores.items()}
    weights = np.empty(len(variants), dtype=np.float32)
    for i, variant in enumerate(variants):
        key = variant[-2:].upper()  # two 3'-most PAM bases
        try:
            weights[i] = pam_scores_u[key]
        except KeyError:
            loggers.errorlog.log_raise_exception(
                f"CFD PAM table has no entry for {key!r} (from variant "
                f"{variant!r}, pam_id {i}). Known keys: {sorted(pam_scores_u)}",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
    return weights


# ==============================================================================
# vectorised scoring kernel
# ==============================================================================


def _score_batch_vectorized(
    rguide: np.ndarray,
    rseq: np.ndarray,
    mm_table: np.ndarray,
    seq_len: int,
    n_rows: int,
) -> np.ndarray:
    """
    Compute the mismatch component of CFD for all rows simultaneously.

    ``rguide``/``rseq`` are **ASCII** (see module docstring). Each column is
    decoded to base indices via :data:`_ASCII_TO_IDX`; a position contributes a
    mismatch weight only where both bases are known (A/C/G/T) and differ. Gaps
    (``-``), padding (``0x00``), and ambiguous bytes decode to ``-1`` and are
    skipped with weight ``1.0``.

    Only the first ``min(_CFD_POSITIONS, seq_len)`` columns are scored. This
    assumes the protospacer is laid out 5'->3' with the PAM-distal base at
    column 0 (matching the ``<protospacer>NGG`` orientation), so column ``c``
    corresponds to CFD position ``c + 1``.

    Parameters
    ----------
    rguide, rseq : np.ndarray
        Shape ``(n_rows, seq_len)``, dtype ``uint8``. ASCII bytes.
    mm_table : np.ndarray
        Shape ``(20, 4, 4)`` pre-built mismatch weight table.
    seq_len : int
        Width of the resolved buffer (columns beyond the protospacer are padding).
    n_rows : int
        Number of alignment rows.

    Returns
    -------
    np.ndarray
        Shape ``(n_rows,)`` float32 — product of mismatch weights per row.
    """
    scores = np.ones(n_rows, dtype=np.float32)
    n_pos = min(_CFD_POSITIONS, seq_len)
    for col in range(n_pos):
        g_idx = _ASCII_TO_IDX[rguide[:, col]]  # (n_rows,) int8, -1 = skip
        t_idx = _ASCII_TO_IDX[rseq[:, col]]
        # score positions where both bases are known and mismatched
        mism = (g_idx >= 0) & (t_idx >= 0) & (g_idx != t_idx)
        if not mism.any():
            continue
        scores[mism] *= mm_table[col, g_idx[mism], t_idx[mism]]
    return scores


# ==============================================================================
# callable scorer class
# ==============================================================================


class CfdScorer(Transformer):
    """
    Callable transform that computes CFD scores for a batch of alignments and
    writes them to score slot 0 in-place.

    Implements the ``Transformer`` protocol: the pipeline delivers a **raw**
    ``PyAlignmentBatch``; ``__call__`` wraps it in an
    :class:`~crisprme2.crisprme_core_api.AlignmentBatch`, reads ``rguide`` /
    ``rseq`` / ``pam_id``, computes CFD vectorised over all rows, and fills
    ``batch.score(CFD_SCORE_SLOT)`` without allocating extra buffers.

    The run's (possibly degenerate) PAM motif is supplied at construction. Its
    concrete variants are enumerated once — in ``pam_id`` order — so the PAM
    weight can be gathered per row at call time.

    Parameters
    ----------
    pam : str
        The run's IUPAC PAM motif (e.g. ``"NGG"``, ``"NRG"``). Not the scored
        2-base region — the full motif, matching what the batcher was built with.
    loggers : CrisprmeLoggers
        Shared logger bundle for error propagation and debug logging.

    Raises
    ------
    Crisprme2CfdScoreError
        If the native extension is unavailable, the model files cannot be
        loaded, or a PAM variant maps to a key absent from the PAM score table.

    Examples
    --------
    ::

        from crisprme2.scores.cfd_score import CfdScorer

        scorer = CfdScorer(pam="NGG", loggers=loggers)
        # passed as a transform: Pipeline.create(transforms=[scorer])
    """

    def __init__(self, pam: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers
        self._pam_motif: str = pam.upper()
        if pam_variants_ascii is None:
            loggers.errorlog.log_raise_exception(
                "Native function 'pam_variants_ascii' is unavailable; the "
                "compiled extension is required for CFD PAM scoring. Rebuild "
                "with `maturin develop`.",
                os.EX_UNAVAILABLE,
                Crisprme2CfdScoreError,
            )
        # load + pre-process models at construction time, not per call
        mm_scores, pam_scores = load_models(loggers)
        self._mm_table: np.ndarray = _build_mm_lookup(mm_scores, loggers)
        # enumerate concrete PAM variants in pam_id order (Rust is the source of
        # truth), then fold each into a per-id CFD PAM weight -- computed once.
        try:
            variants: List[str] = list(pam_variants_ascii(self._pam_motif))
        except Exception as e:
            loggers.errorlog.log_raise_exception(
                f"Failed enumerating PAM variants for {self._pam_motif!r}: {e}",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        self._pam_weight_by_id: np.ndarray = _build_pam_weight_by_id(
            variants, pam_scores, loggers
        )
        loggers.verboselog.debug(
            f"{self.__class__.__name__} initialized (pam={self._pam_motif!r}, "
            f"{self._pam_weight_by_id.shape[0]} variant(s))"
        )

    # --------------------------------------------------------------------------
    # internal helpers
    # --------------------------------------------------------------------------

    def _pam_weight_for_batch(self, pam_id: np.ndarray) -> np.ndarray:
        """
        Gather the per-row PAM weight from ``pam_id``.

        ``PAM_ID_NONE`` and any out-of-range id fall back to
        :data:`_PAM_NONE_WEIGHT`.

        Parameters
        ----------
        pam_id : np.ndarray
            Shape ``(N,)`` uint16, the batch's ``pam_id`` column.

        Returns
        -------
        np.ndarray
            Shape ``(N,)`` float32 PAM weights.
        """
        n_variants = self._pam_weight_by_id.shape[0]
        valid = pam_id < n_variants  # excludes PAM_ID_NONE (0xFFFF) and OOB ids
        safe = np.where(valid, pam_id, 0).astype(np.intp)
        pam_w = self._pam_weight_by_id[safe]  # fancy index -> owned copy
        if not valid.all():
            pam_w[~valid] = _PAM_NONE_WEIGHT
        return pam_w

    # --------------------------------------------------------------------------
    # transformer protocol
    # --------------------------------------------------------------------------

    def __call__(self, raw_batch: Any) -> None:
        """
        Score all alignments in *raw_batch* using the CFD model.

        Wraps the raw ``PyAlignmentBatch`` in an
        :class:`~crisprme2.crisprme_core_api.AlignmentBatch`, computes CFD
        vectorised over all rows, and writes results to
        ``score[CFD_SCORE_SLOT]`` in-place.

        Parameters
        ----------
        raw_batch : PyAlignmentBatch
            Raw Rust alignment batch delivered by the pipeline stage.

        Raises
        ------
        Crisprme2CfdScoreError
            If the batch cannot be wrapped or scoring fails.
        """
        batch = AlignmentBatch(raw_batch, self._loggers)
        n_rows = batch.n_rows
        if n_rows == 0:
            return  # nothing to score
        # rguide/rseq arrive as 2-D (n_rows, W) arrays: the Rust `[u8; 32]`
        # columns are exposed via `PyBuffer::from_array`, which presents a 2-D
        # buffer (n_rows rows x 32 bytes), NOT a flat vector. Normalise by total
        # size
        #
        # reshape(n_rows, -1) yields (n_rows, W) whether the buffer
        # arrives 2-D or flat, and W is then the resolved buffer width (32).
        rguide_raw = batch.rguide  # (n_rows, W) uint8, read-only
        rseq_raw = batch.rseq  # (n_rows, W) uint8, read-only
        total_bytes = rguide_raw.size
        if total_bytes % n_rows != 0:
            self._loggers.errorlog.log_raise_exception(
                f"{self.__class__.__name__}: rguide element count ({total_bytes}) "
                f"is not divisible by n_rows ({n_rows}). Buffer layout unexpected.",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        try:
            rguide = rguide_raw.reshape(n_rows, -1)
            rseq = rseq_raw.reshape(n_rows, -1)
        except ValueError as e:
            self._loggers.errorlog.log_raise_exception(
                f"{self.__class__.__name__}: failed to reshape sequence buffers: {e}",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        seq_len = rguide.shape[1]
        try:  # mismatch component, vectorised
            mm_scores = _score_batch_vectorized(
                rguide, rseq, self._mm_table, seq_len, n_rows
            )
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"{self.__class__.__name__}: vectorised mismatch scoring failed: {e}",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        # per-row PAM weight, gathered from pam_id
        pam_weights = self._pam_weight_for_batch(batch.pam_id)
        cfd_scores = (mm_scores * pam_weights).astype(np.float32, copy=False)
        # write in-place to the Rust-owned score buffer (no copy)
        out = batch.score(CFD_SCORE_SLOT)
        out[:] = cfd_scores
        self._loggers.verboselog.debug(
            f"{self.__class__.__name__}: scored {n_rows} alignments "
            f"(mean={float(cfd_scores.mean()):.4f})"
        )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(pam={self._pam_motif!r}, "
            f"variants={self._pam_weight_by_id.shape[0]}, "
            f"slot={CFD_SCORE_SLOT})"
        )
