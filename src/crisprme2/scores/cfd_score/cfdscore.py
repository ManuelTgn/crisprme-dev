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

    cfd = PAM_score(pam) x mismatch_score(rU:dX, position)

where the product runs over all positions where guide ≠ target (gaps from
bulge alignments are skipped), and positions are 1-indexed from the PAM-
proximal end.

Vectorised implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``__call__`` path avoids Python loops over rows entirely.  IUPAC
bitmask arrays (``rguide``, ``rseq``) are read as flat ``uint8`` views,
reshaped to ``(N, seq_len)``, and processed column-by-column so that NumPy
SIMD paths do the heavy lifting.  The PAM score look-up is a single
vectorised index into a pre-built array.

Score slot convention
~~~~~~~~~~~~~~~~~~~~~
CFD scores are written to **score slot 0** (``batch.score(0)``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np

import os
import pickle

from ...crisprme_core_api import AlignmentBatch
from ..crisprme2_scores_error import Crisprme2CfdScoreError
from ...dna_alphabet import RNA, RC, dna2rna_nt
from ...logger import CrisprmeLoggers


# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------

_MODELS_DIR: Path = Path(__file__).parent / "models"
_MM_SCORES_FILE: Path = _MODELS_DIR / "mismatch_score.pkl"
_PAM_SCORES_FILE: Path = _MODELS_DIR / "pam_scores.pkl"

# score slot index for CFD within PyAlignmentBatch
CFD_SCORE_SLOT: int = 0

# number of positions scored by the CFD model (PAM-distal position 1 to 20)
_CFD_POSITIONS: int = 20


# ------------------------------------------------------------------------------
# model types
# ------------------------------------------------------------------------------

MismatchScores = Dict[str, float]
PamScores = Dict[str, float]


# ------------------------------------------------------------------------------
# model loading
# ------------------------------------------------------------------------------


def load_models(loggers: CrisprmeLoggers) -> Tuple[MismatchScores, PamScores]:
    """
    Load CFD mismatch and PAM score tables from the bundled pickle files.

    The models are the original tables from Doench et al. (2016).  Both
    files are shipped with the package.
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
    except (FileNotFoundError, pickle.UnpicklingError, Exception) as e:
        loggers.errorlog.log_raise_exception(
            f"Failed to load CFD mismatch scores from {_MM_SCORES_FILE}: {e}",
            os.EX_IOERR,
            Crisprme2CfdScoreError,
        )
    try:  # load pam score model
        with open(_PAM_SCORES_FILE, mode="rb") as fh:
            pam_scores: PamScores = pickle.load(fh)
    except (FileNotFoundError, pickle.UnpicklingError, Exception) as e:
        loggers.errorlog.log_raise_exception(
            f"Failed to load CFD PAM scores from {_PAM_SCORES_FILE}: {e}",
            os.EX_IOERR,
            Crisprme2CfdScoreError,
        )
    return mm_scores, pam_scores


# ------------------------------------------------------------------------------
# pre-computed lookup tables
# ------------------------------------------------------------------------------


def _build_mm_lookup(mm_scores: MismatchScores, loggers: CrisprmeLoggers) -> np.ndarray:
    """
    Build a dense ``(20, 4, 4)`` float32 lookup table from the sparse
    mismatch score dict.

    Axes: ``[position-1 (0..19), guide_base (A/C/G/T), target_base (A/C/G/T)]``

    Missing entries default to ``1.0`` (perfect-match weight), which is
    safe — guide == target positions are never keyed in the model anyway.

    The CFD key format is ``"rX:dY,pos"`` where X is the RNA (guide) base
    and Y is the DNA (target) base, using the convention that the DNA base
    is reverse-complemented before look-up.

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
    # base -> column index (RNA alphabet after dna2rna)
    base_idx = {nt: i for i, nt in enumerate(RNA)}
    table = np.ones((20, 4, 4), dtype=np.float32)
    for key, val in mm_scores.items():
        try:  # key format: "rX:dY,pos" e.g. "rU:dT,12"
            pair, pos_str = key.split(",")
            pos = int(pos_str) - 1  # convert to 0-based
            rna_nt = pair[1]  # 'r' + base
            dna_nt = pair[4]  # 'd' + base (before RC in key construction)
            # the key stores the RC of the target base
            # so we reverse-complement dna_nt back to get the actual target base
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


def _build_pam_lookup(pam_scores: PamScores) -> Dict[str, float]:
    """
    Return the PAM score dict with keys normalised to uppercase.

    No structural transformation is needed, the dict is small (16 entries)
    and look-ups happen once per batch, not once per row.

    Parameters
    ----------
    pam_scores : PamScores
        Raw dict loaded from ``pam_scores.pkl``.

    Returns
    -------
    dict[str, float]
        Uppercase-keyed PAM score dict.
    """
    return {k.upper(): v for k, v in pam_scores.items()}


# ------------------------------------------------------------------------------
# vectorised scoring kernels
# ------------------------------------------------------------------------------


def _score_batch_vectorized(
    rguide: np.ndarray,
    rseq: np.ndarray,
    mm_table: np.ndarray,
    seq_len: int,
    n_rows: int,
) -> np.ndarray:
    """
    Compute the mismatch component of CFD for all rows simultaneously.

    Parameters
    ----------
    rguide : np.ndarray
        Shape ``(n_rows, seq_len)``, dtype ``uint8``. IUPAC bitmask bytes.
    rseq : np.ndarray
        Shape ``(n_rows, seq_len)``, dtype ``uint8``. IUPAC bitmask bytes.
    mm_table : np.ndarray
        Shape ``(20, 4, 4)`` pre-built mismatch weight table.
    seq_len : int
        Total sequence length including PAM; only positions 0..19 are
        scored (PAM-distal protospacer).
    n_rows : int
        Number of alignment rows.

    Returns
    -------
    np.ndarray
        Shape ``(n_rows,)`` float32 — product of mismatch weights per row.
    """
    scores = np.ones(n_rows, dtype=np.float32)
    # IUPAC bitmask → base index mapping (A=1,C=2,G=4,T=8 in the Rust encoder)
    # We map each bitmask byte to an index 0..3 for the lookup table.
    # Single-base bitmasks: A=0x01, C=0x02, G=0x04, T=0x08 (or U after RNA conv)
    bitmask_to_idx = np.full(256, -1, dtype=np.int8)
    bitmask_to_idx[0b0001] = 0  # A
    bitmask_to_idx[0b0010] = 1  # C
    bitmask_to_idx[0b0100] = 2  # G
    bitmask_to_idx[0b1000] = 3  # T / U
    for col in range(min(_CFD_POSITIONS, seq_len)):
        g_bytes = rguide[:, col]  # (n_rows,) uint8
        t_bytes = rseq[:, col]  # (n_rows,) uint8
        # rows where guide == target → weight 1.0, skip
        mismatch_mask = g_bytes != t_bytes
        if not mismatch_mask.any():
            continue
        # rows with gap (0b0000) in either sequence -> bulge, skip (weight 1.0)
        valid_mask = mismatch_mask & (g_bytes != 0) & (t_bytes != 0)
        if not valid_mask.any():
            continue
        g_idx = bitmask_to_idx[g_bytes[valid_mask].astype(np.intp)]
        t_idx = bitmask_to_idx[t_bytes[valid_mask].astype(np.intp)]
        # rows where either base is ambiguous (index == -1) -> skip
        known = (g_idx >= 0) & (t_idx >= 0)
        if not known.any():
            continue
        rows_known = np.where(valid_mask)[0][known]
        weights = mm_table[col, g_idx[known], t_idx[known]]
        scores[rows_known] *= weights
    return scores


# ------------------------------------------------------------------------------
# callable scorer class
# ------------------------------------------------------------------------------


class CfdScorer:
    """
    Callable transform that computes CFD scores for a batch of alignments
    and writes them to score slot 0 in-place.

    Implements the :class:`~crisprme2.protocol.Transformer` protocol:
    accepts a raw ``PyAlignmentBatch``, wraps it in
    :class:`~crisprme2.crisprme_core_api.AlignmentBatch`, and fills
    ``batch.score(CFD_SCORE_SLOT)`` without allocating extra memory.

    The PAM sequence must be provided at construction time so the PAM
    weight can be applied per batch.  The PAM is the two nucleotides
    immediately 3′ of the protospacer for SpCas9 (e.g. ``"GG"`` for NGG).

    Parameters
    ----------
    pam : str
        Two-character PAM string (uppercase) used for PAM score look-up.
        Must exist in the CFD PAM score table.
    loggers : CrisprmeLoggers
        Shared logger bundle for error propagation and debug logging.

    Raises
    ------
    Crisprme2CfdScoreError
        If the model files cannot be loaded, or if *pam* is not found in
        the PAM score table.

    Examples
    --------
    ::

        from crisprme2.scores.cfd_score import CfdScorer

        scorer = CfdScorer(pam="GG", loggers=loggers)
        # scorer is passed as a transform to Pipeline.create(transforms=[scorer])
    """

    def __init__(self, pam: str, loggers: CrisprmeLoggers) -> None:
        self._loggers = loggers
        self._pam = pam
        # load and pre-process models at construction time, not at call time
        mm_scores, pam_scores = load_models(loggers)
        self._mm_table: np.ndarray = _build_mm_lookup(mm_scores, loggers)
        self._pam_scores: Dict[str, float] = _build_pam_lookup(pam_scores)
        # validate PAM at construction so failure is immediate, not buried
        # inside the first batch call
        if self._pam not in self._pam_scores:
            loggers.errorlog.log_raise_exception(
                f"PAM {self._pam!r} not found in CFD PAM score table. "
                f"Available PAMs: {sorted(self._pam_scores)}",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        self._pam_weight: float = self._pam_scores[self._pam]
        loggers.verboselog.debug(
            f"{self.__class__.__name__} initialized (pam={self._pam!r}, "
            f"pam_weight={self._pam_weight:.4f})"
        )

    # --------------------------------------------------------------------------
    # transformer protocol
    # --------------------------------------------------------------------------

    def __call__(self, raw_batch: object) -> Any:
        """
        Score all alignments in *raw_batch* using the CFD model.

        Wraps the raw ``PyAlignmentBatch`` in an
        :class:`~crisprme2.crisprme_core_api.AlignmentBatch` view, reads
        ``rguide`` and ``rseq`` as read-only NumPy arrays, computes CFD
        scores vectorised over all rows, and writes results to
        ``score[CFD_SCORE_SLOT]`` in-place.

        Parameters
        ----------
        raw_batch : PyAlignmentBatch
            Opaque Rust alignment batch delivered by the pipeline stage.

        Raises
        ------
        Crisprme2CfdScoreError
            If the batch cannot be wrapped or scoring fails.
        """
        try:
            batch = AlignmentBatch(raw_batch, self._loggers)
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"CfdScorer: failed to wrap PyAlignmentBatch: {e}",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        n_rows = batch.n_rows
        if n_rows == 0:
            return  # nothing to score
        # flat IUPAC bitmask arrays -> reshape to (N, seq_len)
        rguide_flat = batch.rguide  # (N * seq_len,) uint8, read-only
        rseq_flat = batch.rseq  # (N * seq_len,) uint8, read-only
        total_bytes = rguide_flat.shape[0]
        if total_bytes % n_rows != 0:
            self._loggers.errorlog.log_raise_exception(
                f"{self.__class__.__name__}: rguide byte count ({total_bytes}) "
                f"is not divisible by n_rows ({n_rows}). Buffer layout unexpected",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        seq_len = total_bytes // n_rows
        try:
            rguide = rguide_flat.reshape(n_rows, seq_len)
            rseq = rseq_flat.reshape(n_rows, seq_len)
        except ValueError as e:
            self._loggers.errorlog.log_raise_exception(
                f"CfdScorer: failed to reshape sequence buffers: {e}",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        try:  # compute mismatch component vectorised
            mm_scores = _score_batch_vectorized(
                rguide, rseq, self._mm_table, seq_len, n_rows
            )
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"CfdScorer: vectorised mismatch scoring failed: {e}",
                os.EX_DATAERR,
                Crisprme2CfdScoreError,
            )
        # apply PAM weight (scalar broadcast)
        cfd_scores = mm_scores * self._pam_weight
        # write in-place to the Rust-owned score buffer (no copy)
        out = batch.score(CFD_SCORE_SLOT)
        out[:] = cfd_scores
        self._loggers.verboselog.debug(
            f"{self.__class__.__name__}: scored {n_rows} alignments "
            f"(mean={float(cfd_scores.mean()):.4f}, "
            f"pam_weight={self._pam_weight:.4f})"
        )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(pam={self._pam}, "
            f"pam_weight=f{self._pam_weight:.4f}, "
            f"slot={CFD_SCORE_SLOT})"
        )
