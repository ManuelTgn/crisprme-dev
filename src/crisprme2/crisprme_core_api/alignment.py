"""
alignment.py
------------
Python wrapper around the Rust ``PyAlignmentBatch`` struct exposed via PyO3.

Memory model
~~~~~~~~~~~~
Every field on ``PyAlignmentBatch`` is a ``PyBuffer`` — a zero-copy view into a
contiguous, Rust-owned memory region.  The wrapper turns each buffer into a
NumPy array via ``np.asarray`` **without copying**:

- Read-only fields (``seq_id``, ``offset``, ``rguide``, ``rseq``) are wrapped
  and their ``writeable`` flag is immediately cleared, so a transform cannot
  corrupt the alignment records it is scoring.
- The mutable field (``score``) is wrapped as a writeable array: transforms are
  expected to fill it in-place.

Buffer shapes
~~~~~~~~~~~~~
::

    seq_id    : np.ndarray[uint32,  1-D, read-only]  - (N,)               window id per row
    offset    : np.ndarray[uint32,  1-D, read-only]  - (N,)               genomic offset per row
    rguide    : np.ndarray[uint8,   1-D, read-only]  - (N * SEQ_MAX_LEN,) IUPAC bitmasks, flat
    rseq      : np.ndarray[uint8,   1-D, read-only]  - (N * SEQ_MAX_LEN,) IUPAC bitmasks, flat
    score(i)  : np.ndarray[float32, 1-D, writeable]  - (N,)               score slot i ∈ [0, N_SCORE_SLOTS)

``rguide`` and ``rseq`` are flat byte arrays.  To recover per-row sequences,
view them as fixed-width byte strings::

    seq_len = 32          # or whatever SEQ_MAX_LEN is for this run
    seqs = batch.rseq.view(f"S{seq_len}")   # shape (N,), dtype '|S32'

.. note::
    The native ``PyAlignmentBatch`` also carries per-row *feature/annotation*
    slots.  These are deliberately **not** surfaced by this wrapper yet; the
    annotation API lands in a later commit.

.. warning::
    Every array returned by this class is valid **only** for the lifetime of
    the ``PyAlignmentBatch`` object delivered to the owning transform's
    ``__call__``.  Never keep a reference to a returned array beyond a single
    transform invocation — the Rust pipeline may reclaim or overwrite the
    underlying memory the moment ``__call__`` returns.

Error handling
~~~~~~~~~~~~~~
All failures are routed through ``loggers.errorlog.log_raise_exception`` rather
than being raised directly.  That path logs the error (with traceback) to
``errors.log``, flushes the halt banner to ``stderr``, and terminates the run
with the supplied ``os.EX_*`` exit code.  Because it does not return, the
wrapper never leaves a partially-constructed array in the caller's hands.

Typical usage (inside a scoring transform)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    from crisprme2.crisprme_core_api import AlignmentBatch

    class CfdScorer:
        def __init__(self, loggers):
            self._loggers = loggers

        def __call__(self, raw_batch) -> None:
            # the pipeline hands the transform a *raw* PyAlignmentBatch;
            # wrap it before touching any buffer.
            batch  = AlignmentBatch(raw_batch, self._loggers)
            guide  = batch.rguide             # (N * SEQ_MAX_LEN,), uint8, read-only
            target = batch.rseq               # (N * SEQ_MAX_LEN,), uint8, read-only
            out    = batch.score(0)           # (N,), float32, writeable
            out[:] = _compute_cfd(guide, target)
"""

from __future__ import annotations

from typing import Any

import numpy as np

import os


from ..logger import CrisprmeLoggers

from .crisprme2_api_error import Crisprme2AlignmentBatchError

try:
    from .._crisprme2_native import PyAlignmentBatch as RustAlignmentBatch
except ImportError:
    # fallback for development/testing before the native extension is built
    RustAlignmentBatch = None


# ==============================================================================
# fixed slot counts - must match the native transform stage (transform.rs)
# ==============================================================================

#: number of score slots on PyAlignmentBatch (Rust: ``scores: [PyBuffer; 4]``)
N_SCORE_SLOTS: int = 4


# ==============================================================================
# dtype constants
# ==============================================================================

_DTYPE_U32 = np.dtype(np.uint32)
_DTYPE_U16 = np.dtype(np.uint16)
_DTYPE_U8 = np.dtype(np.uint8)
_DTYPE_F32 = np.dtype(np.float32)


# ==============================================================================
# internal helpers
# ==============================================================================


def _require_native(loggers: CrisprmeLoggers) -> None:
    """Halt the run if the native extension has not been compiled."""
    if RustAlignmentBatch is None:
        loggers.errorlog.log_raise_exception(
            "Rust PyAlignmentBatch type not exposed to Python. Ensure the "
            "native extension (_crisprme2_native) is compiled and installed.",
            os.EX_CANTCREAT,
            Crisprme2AlignmentBatchError,
        )


def _validate_raw_batch(raw: Any, loggers: CrisprmeLoggers) -> None:
    """
    Confirm that *raw* is an instance of the Rust ``PyAlignmentBatch`` type.

    Runs before any buffer access so a misuse surfaces as a typed, logged
    Python exception instead of an opaque PyO3 error deep inside an accessor.
    Skipped when the native type is unavailable (development builds), in which
    case :func:`_require_native` has already halted.
    """
    if RustAlignmentBatch is not None and not isinstance(raw, RustAlignmentBatch):
        loggers.errorlog.log_raise_exception(
            "'raw_batch' must be a PyAlignmentBatch instance, got "
            f"{type(raw).__name__!r}. Instances are created exclusively by the "
            "Rust pipeline stage.",
            os.EX_DATAERR,
            Crisprme2AlignmentBatchError,
        )


def _validate_score_idx(idx: int, loggers: CrisprmeLoggers) -> None:
    """Halt if *idx* is not an ``int`` in ``[0, N_SCORE_SLOTS)``."""
    # ``bool`` is a subclass of ``int``; reject it explicitly to avoid
    # ``score(True)`` silently selecting slot 1.
    if isinstance(idx, bool) or not isinstance(idx, int):
        loggers.errorlog.log_raise_exception(
            f"Score index must be an int, got {type(idx).__name__!r}.",
            os.EX_DATAERR,
            Crisprme2AlignmentBatchError,
        )
    if not (0 <= idx < N_SCORE_SLOTS):
        loggers.errorlog.log_raise_exception(
            f"Score index {idx} out of range - valid range is [0, {N_SCORE_SLOTS}).",
            os.EX_DATAERR,
            Crisprme2AlignmentBatchError,
        )


def _buf_to_readonly(buf: Any, dtype: np.dtype) -> np.ndarray:
    """
    Convert a ``PyBuffer`` to a read-only NumPy array without copying.

    Consumes the buffer protocol via ``np.asarray`` (sharing memory with the
    Rust allocation), then clears the ``writeable`` flag so downstream
    transforms cannot mutate immutable alignment records.

    Parameters
    ----------
    buf : PyBuffer
        Zero-copy buffer returned by a ``PyAlignmentBatch`` accessor method.
    dtype : np.dtype
        Element dtype to apply.  The buffer's byte length must be an integer
        multiple of ``dtype.itemsize``.

    Returns
    -------
    np.ndarray
        A 1-D read-only array sharing memory with the Rust allocation.
    """
    arr = np.asarray(buf, dtype=dtype)
    arr.flags.writeable = False
    return arr


def _buf_to_writable(buf: Any, dtype: np.dtype) -> np.ndarray:
    """
    Convert a ``PyBuffer`` to a writeable NumPy array without copying.

    The buffer must originate from a ``&mut`` slice on the Rust side (true for
    every ``score`` buffer on ``PyAlignmentBatch``).  Forcing ``writeable`` to
    ``True`` doubles as an invariant check: if the buffer is unexpectedly
    read-only, NumPy raises here and the caller routes the failure through
    ``errorlog``.

    Parameters
    ----------
    buf : PyBuffer
        Mutable zero-copy buffer returned by a ``PyAlignmentBatch`` accessor.
    dtype : np.dtype
        Element dtype to apply.

    Returns
    -------
    np.ndarray
        A 1-D writeable array sharing memory with the Rust allocation.
    """
    arr = np.asarray(buf, dtype=dtype)
    arr.flags.writeable = True
    return arr


# ==============================================================================
# public wrapper
# ==============================================================================


class AlignmentBatch:
    """
    Zero-copy, transform-facing view over a native ``PyAlignmentBatch``.

    Wraps the raw PyO3 object handed to a transform's ``__call__`` and exposes
    its columns as NumPy arrays that share memory with the Rust pipeline.  The
    read-only columns (``seq_id``, ``offset``, ``rguide``, ``rseq``) describe
    each alignment; the writeable ``score`` slots are where a transform records
    its results in-place.

    Instances are cheap and stateless beyond the wrapped handle: construct one
    per ``__call__`` and let it fall out of scope when the transform returns.

    Parameters
    ----------
    raw_batch : PyAlignmentBatch
        The native batch delivered by the ``PyTransform`` pipeline stage.
    loggers : CrisprmeLoggers
        Shared logger bundle used for all error propagation.
    """

    __slots__ = ("_raw", "_loggers")

    def __init__(self, raw_batch: Any, loggers: CrisprmeLoggers) -> None:
        _require_native(loggers)
        _validate_raw_batch(raw_batch, loggers)
        self._raw = raw_batch
        self._loggers = loggers

    # ==========================================================================
    # read-only fields
    # ==========================================================================

    @property
    def seq_id(self) -> np.ndarray:
        """
        Window id for each alignment row.

        Shape : ``(N,)``
        Dtype : ``uint32``
        Access: read-only

        Maps each row back to its originating unique window in the
        ``TargetBatcher`` map, used for occurrence look-up after scoring.
        """
        try:
            return _buf_to_readonly(self._raw.seq_id(), _DTYPE_U32)
        except Crisprme2AlignmentBatchError:
            raise
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed accessing seq_id buffer: {e}",
                os.EX_IOERR,
                Crisprme2AlignmentBatchError,
            )

    @property
    def offset(self) -> np.ndarray:
        """
        Genomic offset (absolute position within the contig) for each row.

        Shape : ``(N,)``
        Dtype : ``uint32``
        Access: read-only

        Corresponds to the ``pos`` field unpacked from the Rust ``Occ``
        u64 occurrence record.
        """
        try:
            return _buf_to_readonly(self._raw.offset(), _DTYPE_U32)
        except Crisprme2AlignmentBatchError:
            raise
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed accessing offset buffer: {e}",
                os.EX_IOERR,
                Crisprme2AlignmentBatchError,
            )

    @property
    def rguide(self) -> np.ndarray:
        """
        Aligned guide sequence for all rows, encoded as IUPAC bitmasks.

        Shape : ``(N * SEQ_MAX_LEN,)``
        Dtype : ``uint8``
        Access: read-only

        Flat byte array.  To recover per-row sequences, view it as fixed-width
        byte strings::

            seq_len = 32
            rows = batch.rguide.view(f"S{seq_len}")  # shape (N,)

        Gaps introduced by bulge alignment are encoded as ``0x00``.
        """
        try:
            return _buf_to_readonly(self._raw.rguide(), _DTYPE_U8)
        except Crisprme2AlignmentBatchError:
            raise
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed accessing rguide buffer: {e}",
                os.EX_IOERR,
                Crisprme2AlignmentBatchError,
            )

    @property
    def rseq(self) -> np.ndarray:
        """
        Aligned off-target sequence for all rows, encoded as IUPAC bitmasks.

        Shape : ``(N * SEQ_MAX_LEN,)``
        Dtype : ``uint8``
        Access: read-only

        Parallel to :attr:`rguide`.  To recover per-row sequences::

            seq_len = 32
            rows = batch.rseq.view(f"S{seq_len}")   # shape (N,)

        Gaps introduced by bulge alignment are encoded as ``0x00``.
        """
        try:
            return _buf_to_readonly(self._raw.rseq(), _DTYPE_U8)
        except Crisprme2AlignmentBatchError:
            raise
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed accessing rseq buffer: {e}",
                os.EX_IOERR,
                Crisprme2AlignmentBatchError,
            )

    @property
    def pam_id(self) -> np.ndarray:
        """(N,) uint16, read-only. Concrete PAM-variant index per row; indexes the
        run PAM's variant table (PAM::pam_variant_ascii order). 0xFFFF (PAM_ID_NONE)
        marks rows with no concrete PAM."""
        try:
            return _buf_to_readonly(self._raw.pam_id(), _DTYPE_U16)
        except Crisprme2AlignmentBatchError:
            raise
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed accessing pam_id buffer: {e}",
                os.EX_IOERR,
                Crisprme2AlignmentBatchError,
            )

    # ==========================================================================
    # writeable fields
    # ==========================================================================

    def score(self, idx: int) -> np.ndarray:
        """
        Return the writeable score array for slot *idx*.

        Shape : ``(N,)``
        Dtype : ``float32``
        Access: writeable — assign results in-place.

        Score slot assignments (by convention):

        ===  ==================
        idx  Score model
        ===  ==================
        0    CFD
        1    (reserved)
        2    (reserved)
        3    (reserved)
        ===  ==================

        Parameters
        ----------
        idx : int
            Score slot index in ``[0, N_SCORE_SLOTS)``.

        Returns
        -------
        np.ndarray
            Shape ``(N,)`` float32, writeable, sharing memory with Rust.

        Raises
        ------
        Crisprme2AlignmentBatchError
            If *idx* is out of range or the buffer cannot be accessed.  (Raised
            via ``errorlog``, which logs and halts the run.)

        Examples
        --------
        ::

            scores = batch.score(0)   # CFD slot
            scores[:] = cfd_values    # in-place assignment
        """
        _validate_score_idx(idx, self._loggers)
        try:
            return _buf_to_writable(self._raw.score(idx), _DTYPE_F32)
        except Crisprme2AlignmentBatchError:
            raise
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed accessing score[{idx}] buffer: {e}",
                os.EX_IOERR,
                Crisprme2AlignmentBatchError,
            )

    # ==========================================================================
    # convenience helpers
    # ==========================================================================

    @property
    def n_rows(self) -> int:
        """
        Number of alignment rows in this batch.

        Derived from ``seq_id.shape[0]``.  An empty batch legitimately reports
        ``0``.  A genuine failure to read ``seq_id`` is **not** swallowed: it is
        surfaced by the :attr:`seq_id` accessor, which logs the error and halts
        the run after flushing stderr.
        """
        seq_id = self.seq_id  # halts via errorlog on a real buffer fault
        # a well-formed (possibly empty) buffer is always 1-D
        if seq_id.ndim != 1:
            self._loggers.errorlog.log_raise_exception(
                f"seq_id buffer is {seq_id.ndim}-D, expected 1-D; "
                "cannot derive n_rows.",
                os.EX_DATAERR,
                Crisprme2AlignmentBatchError,
            )
        return int(seq_id.shape[0])

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(rows={self.n_rows}, "
            f"score_slots={N_SCORE_SLOTS})"
        )
