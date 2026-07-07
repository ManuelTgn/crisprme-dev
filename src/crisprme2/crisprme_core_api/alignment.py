"""
alignment.py
------------
Python wrapper around the Rust ``PyAlignmentBatch`` struct exposed via PyO3.

Memory model
~~~~~~~~~~~~
Every field in ``PyAlignmentBatch`` is a ``PyBuffer`` — a zero-copy view into
a contiguous Rust-owned memory region.  The wrapper converts each buffer into
a NumPy array via ``np.asarray`` **without copying**:

- Read-only fields (``seq_id``, ``offset``, ``rguide``, ``rseq``) are wrapped
  and their ``writeable`` flag is immediately cleared so transforms cannot
  corrupt alignment records.
- Mutable fields (``score``, ``feature``) are wrapped as writeable arrays:
  transforms are expected to fill these in-place.

Buffer shapes
~~~~~~~~~~~~~
::

    seq_id     : np.ndarray[uint32, 1-D, read-only]   - (N,)               window id per row
    offset     : np.ndarray[uint32, 1-D, read-only]   - (N,)               genomic offset per row
    rguide     : np.ndarray[uint8,  1-D, read-only]   - (N * SEQ_MAX_LEN,) IUPAC bitmasks, flat
    rseq       : np.ndarray[uint8,  1-D, read-only]   - (N * SEQ_MAX_LEN,) IUPAC bitmasks, flat
    score(i)   : np.ndarray[float32, 1-D, writeable]  - (N,)               score slot i ∈ [0, 4)
    feature(i) : np.ndarray[uint32,  1-D, writeable]  - (N,)               feature slot i ∈ [0, 10)

``rguide`` and ``rseq`` are flat byte arrays.  To recover per-row sequences,
view them as fixed-width byte strings::

    seq_len = 32          # or whatever SEQ_MAX_LEN is for this run
    seqs = batch.rseq.view(f'S{seq_len}')   # shape (N,), dtype '|S32'

.. warning::
    Every array returned by this class is only valid for the lifetime of the
    ``PyAlignmentBatch`` object delivered to the owning transform's
    ``__call__``.  Never store a reference to a returned array beyond a single
    transform invocation — the Rust pipeline may reclaim or overwrite the
    underlying memory immediately after ``__call__`` returns.

Typical usage (inside a scoring transform)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    from crisprme2.crisprme_core_api import AlignmentBatch

    class CfdScorer:
        def __init__(self, loggers):
            self._loggers = loggers

        def __call__(self, raw_batch) -> None:
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
    RustAlignmentBatch = None


# ==============================================================================
# fixed slot counts - match transform.rs constants
# ==============================================================================

#: number of score slots in PyAlignmentBatch (scores: [PyBuffer: 4])
N_SCORE_SLOTS: int = 4

#: number of annotation slots in PyAlignmentBatch (features: [PyBuffer: 4])
N_ANNOTATION_SLOTS: int = 4


# ==============================================================================
# dtype constants
# ==============================================================================

_DTYPE_U32 = np.dtype(np.uint32)
_DTYPE_U8 = np.dtype(np.uint8)
_DTYPE_F32 = np.dtype(np.float32)


# ==============================================================================
# internal helpers
# ==============================================================================


def _require_native(loggers: CrisprmeLoggers) -> None:
    """Raise if the native extension has not been compiled"""
    if RustAlignmentBatch is None:
        loggers.errorlog.log_raise_exception(
            "Rust PyAlignmentBatch type not exposed to Python. Ensure the "
            "native extension (_crisprme2_native) is compiled and installed",
            os.EX_CANTCREAT,
            Crisprme2AlignmentBatchError,
        )


def _validate_raw_batch(raw: Any, loggers: CrisprmeLoggers) -> None:
    """
    Confirm that *raw* is an instance of the Rust ``PyAlignmentBatch`` type.

    This guard runs before any buffer access so errors surface as typed
    Python exceptions rather than opaque PyO3 panics.
    """
    if RustAlignmentBatch is not None and not isinstance(raw, RustAlignmentBatch):
        loggers.errorlog.log_raise_exception(
            "'raw_batch' must be a PyAlignmentBatch instance, got "
            f"{type(raw).__name__!r}. Instances are created exclusively by the "
            "Rust pipeline stage",
            os.EX_DATAERR,
            Crisprme2AlignmentBatchError,
        )


def _validate_score_idx(idx: int, loggers: CrisprmeLoggers) -> None:
    """Raise if *idx* is outside ``[0, N_SCORE_SLOTS)``"""
    if isinstance(idx, bool) or not isinstance(idx, int):
        loggers.errorlog.log_raise_exception(
            f"Score index must be an int, got {type(idx).__name__!r}",
            os.EX_DATAERR,
            Crisprme2AlignmentBatchError,
        )
    if not (0 <= idx < N_SCORE_SLOTS):
        loggers.errorlog.log_raise_exception(
            f"Score index {idx} out of range - valid range is [0, {N_SCORE_SLOTS}]",
            os.EX_DATAERR,
            Crisprme2AlignmentBatchError,
        )


def _validate_annotation_idx(idx: int, loggers: CrisprmeLoggers) -> None:
    """Raise if *idx* is outside ``[0, N_ANNOTATION_SLOTS)``"""
    if isinstance(idx, bool) or not isinstance(idx, int):
        loggers.errorlog.log_raise_exception(
            f"Annotation index must be an int, got {type(idx).__name__!r}",
            os.EX_DATAERR,
            Crisprme2AlignmentBatchError,
        )
    if not (0 <= idx < N_ANNOTATION_SLOTS):
        loggers.errorlog.log_raise_exception(
            f"Annotation index {idx} out of range - valid range is [0, {N_ANNOTATION_SLOTS}]",
            os.EX_DATAERR,
            Crisprme2AlignmentBatchError,
        )


def _buf_to_readonly(buf: Any, dtype: np.dtype) -> np.ndarray:
    """
    Convert a ``PyBuffer`` to a read-only NumPy array without copying.

    Uses ``np.asarray`` to consume the buffer protocol, then immediately
    clears the ``writeable`` flag.  The array shares memory with the
    Rust allocation (no copy is made).

    Parameters
    ----------
    buf : PyBuffer
        Zero-copy buffer returned by a ``PyAlignmentBatch`` accessor method.
    dtype : np.dtype
        Element dtype to apply.  The buffer's byte length must be an
        integer multiple of ``dtype.itemsize``.

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

    The buffer must originate from a ``&mut`` slice on the Rust side
    (which is the case for all ``score`` and ``feature`` buffers in
    ``PyAlignmentBatch``).

    Parameters
    ----------
    buf : PyBuffer
        Mutable zero-copy buffer returned by a ``PyAlignmentBatch``
        accessor method.
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

        The buffer is a flat byte array.  To recover per-row sequences,
        view it as fixed-width byte strings::

            seq_len = 32
            rows = batch.rguide.view(f'S{seq_len}')  # shape (N,)

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
            rows = batch.rseq.view(f'S{seq_len}')   # shape (N,)

        Gaps introduced by bulge alignment are encoded as ``0x00``.
        """
        try:
            return _buf_to_readonly(self._raw, _DTYPE_U8)
        except Crisprme2AlignmentBatchError:
            raise
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed accessing rseq buffer: {e}",
                os.EX_IOERR,
                Crisprme2AlignmentBatchError,
            )

    # ==========================================================================
    # writable mutable fields
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
        1    ?
        2    ?
        3    ?
        ===  ==================

        Parameters
        ----------
        idx : int
            Score slot index in ``[0, N_SCORE_SLOTS)`` i.e. ``[0, 4)``.

        Returns
        -------
        np.ndarray
            Shape ``(N,)`` float32, writeable, sharing memory with Rust.

        Raises
        ------
        Crisprme2AlignmentError
            If *idx* is out of range or the buffer cannot be accessed.

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

    def annotation(self, idx: int) -> np.ndarray:
        """
        Return the writeable feature/annotation bitmask array for slot *idx*.

        Shape : ``(N,)``
        Dtype : ``uint32``
        Access: writeable — set annotation bits in-place.

        Each element is a bitmask where individual bits correspond to
        genomic annotation features (e.g. exon, promoter, repeat region).
        The bit-to-feature mapping is defined by the
        :class:`~crisprme2.crisprme_core_api.FeatureRegistry`.

        Parameters
        ----------
        idx : int
            Feature slot index in ``[0, N_FEATURE_SLOTS)`` i.e. ``[0, 10)``.

        Returns
        -------
        np.ndarray
            Shape ``(N,)`` uint32, writeable, sharing memory with Rust.

        Raises
        ------
        Crisprme2AlignmentError
            If *idx* is out of range or the buffer cannot be accessed.

        Examples
        --------
        ::

            feat = batch.feature(0)
            feat[:] = annotation_bitmasks   # in-place assignment
        """
        _validate_annotation_idx(idx, self._loggers)
        try:
            return _buf_to_writable(self._raw.feature(idx), _DTYPE_U32)
        except Crisprme2AlignmentBatchError:
            raise
        except Exception as e:
            self._loggers.errorlog.log_raise_exception(
                f"Failed accessing feature[{idx}] buffer: {e}",
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

        Derived from ``seq_id.shape[0]``.  Returns ``0`` if the buffer
        cannot be read (e.g. empty batch).
        """
        try:
            return int(self.seq_id.shape[0])
        except Exception:
            return 0

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(rows={self.n_rows}, "
            f"score_slots={N_SCORE_SLOTS}, "
            f"annotation_slots={N_ANNOTATION_SLOTS})"
        )
