//! Zero-copy Python views over ring buffer batches.
//!
//! Both [`AlignmentBatchView`] and [`SequenceBatchView`] follow the same
//! pattern: they wrap an `Option<RingBatch>`, expose the underlying memory
//! to Python via the buffer protocol, and return the ring slot to the pool
//! on drop.
//!
//! ## Buffer protocol and slot lifetime
//!
//! When Python calls `memoryview(view)`, CPython invokes `__getbuffer__`,
//! which increments the object's reference count via `into_ptr()`.  The ring
//! slot is therefore kept alive for as long as any `memoryview` derived from
//! it exists — even if the [`AlignmentBatchView`] Python object itself is
//! garbage-collected first.  The slot is returned to the ring only when the
//! last `memoryview` is released **and** the view is dropped.
//!
//! ## Validity
//!
//! Both views can be in an *empty* state (constructed via `::empty()`), which
//! represents a non-available result from a non-blocking receive.  All
//! `#[pymethods]` that require data check for this state and raise
//! [`PyValueError`] / [`PyBufferError`] accordingly.

use crate::alignment::alignment::Alignment;
use crate::memory::batch::{AlignmentRingBatch, SequenceRingBatch};
use pyo3::exceptions::{PyBufferError, PyValueError};
use pyo3::ffi::Py_buffer;
use pyo3::{ffi, pyclass, pymethods, PyRef, PyResult};
use std::ffi::{c_int, c_void};
use std::ptr;
use tracing::trace;

/// A Python-visible, zero-copy view over an [`AlignmentRingBatch`].
///
/// Produced by [`HybridEngine::receive`] / [`HybridEngine::receive_blocking`]
/// and handed directly to Python.  Exposes the batch's [`Alignment`] array
/// as a read-only [`memoryview`] so NumPy or `struct.unpack_from` can consume
/// results without any copying.
///
/// The view may be *empty* when a non-blocking receive found no data; callers
/// should check [`is_valid`][AlignmentBatchView::is_valid] before use.
///
/// ## Memory layout
///
/// Each element is a 24-byte [`Alignment`] with the following struct format
/// (compatible with Python's `struct` module):
///
/// | Format | Field               | Bytes |
/// |--------|---------------------|-------|
/// | `Q`    | `cigarx.storage`    | 8     |
/// | `B`    | `cigarx.bits`       | 1     |
/// | `xxxxxxx` | padding          | 7     |
/// | `I`    | `id`                | 4     |
/// | `B`    | `offset`            | 1     |
/// | `B`    | `strand`            | 1     |
/// | `xx`   | padding             | 2     |
///
/// Full format string: `"QBxxxxxxxIBBxx"` (24 bytes total).
///
/// ## Python usage
/// ```python
/// import struct
///
/// batch = engine.receive_blocking()
/// if batch.is_valid():
///     view = memoryview(batch)
///     fmt, item = "QBxxxxxxxIBBxx", struct.calcsize("QBxxxxxxxIBBxx")
///     for i in range(len(view) // item):
///         storage, bits, id_, offset, strand = struct.unpack_from(fmt, view, i * item)
/// ```
#[pyclass]
pub struct AlignmentBatchView {
    /// The underlying ring buffer slot. `None` when the view is empty.
    inner: Option<AlignmentRingBatch>,
}

impl AlignmentBatchView {
    /// Wrap a ring buffer slot in a view. The slot will be returned to the
    /// ring when this view is dropped.
    pub fn new(batch: AlignmentRingBatch) -> Self {
        Self { inner: Some(batch) }
    }

    /// Create an empty sentinel view, indicating that no batch is currently
    /// available (e.g. from a non-blocking receive that found nothing).
    pub fn empty() -> Self {
        Self { inner: None }
    }
}

impl Drop for AlignmentBatchView {
    /// Returns the ring buffer slot to the pool.
    fn drop(&mut self) {
        if let Some(inner) = self.inner.take() {
            trace!("Dropping python alignment-view");
            drop(inner);
        }
    }
}

#[pymethods]
impl AlignmentBatchView {
    /// Return the ID of the [`TargetBatcher`] whose sequences produced this
    /// batch.
    ///
    /// Useful for correlating alignment results back to the original query
    /// on the Python side.
    ///
    /// # Errors
    ///
    /// Raises [`PyValueError`] if the view is empty.
    fn batcher_id(&self) -> PyResult<usize> {
        match &self.inner {
            None => Err(PyValueError::new_err("View is empty")),
            Some(inner) => Ok(inner.id()),
        }
    }

    /// Return `True` if the view holds a live alignment batch.
    ///
    /// Always check this before calling [`memoryview`] or
    /// [`batcher_id`][AlignmentBatchView::batcher_id] on a view returned by
    /// the non-blocking [`HybridEngine::receive`].
    fn is_valid(&self) -> PyResult<bool> {
        Ok(self.inner.is_some())
    }

    /// Return the number of alignments in this batch.
    fn size(&self) -> PyResult<usize> {
        match &self.inner {
            None => Err(PyValueError::new_err("View is empty")),
            Some(inner) => Ok(inner.len()),
        }
    }

    /// Implement the Python buffer protocol, exposing the alignment array as
    /// a flat, read-only `memoryview`.
    ///
    /// CPython calls this when the user writes `memoryview(batch)`.
    /// On success, `(*view).obj` is set to a new strong reference to `slf`
    /// (via [`PyRef::into_ptr`]), keeping the ring slot alive until
    /// [`PyBuffer_Release`] is called.
    ///
    /// # Layout
    ///
    /// The buffer is a 1-D array of 24-byte items with format `"QBxxxxxxxIBBxx"`.
    /// See the [`AlignmentBatchView`] type-level docs for the full field table.
    ///
    /// # Errors
    ///
    /// | Condition | Error |
    /// |-----------|-------|
    /// | `view` pointer is null | [`PyBufferError`] |
    /// | Caller requested a writable buffer | [`PyBufferError`] |
    /// | View is empty (no batch attached) | [`PyBufferError`] |
    ///
    /// # Safety
    ///
    /// Unsafe because it writes to a raw `*mut Py_buffer` supplied by CPython.
    /// All pointer fields (`shape`, `strides`, `suboffsets`, `internal`) are
    /// set to well-defined values; `shape` and `strides` point into the view
    /// struct itself, which CPython guarantees outlives this call.
    unsafe fn __getbuffer__(
        slf: PyRef<'_, Self>,
        view: *mut Py_buffer,
        flags: c_int,
    ) -> PyResult<()> {
        if view.is_null() {
            return Err(PyBufferError::new_err("view is null"));
        }

        if (flags & ffi::PyBUF_WRITABLE) == ffi::PyBUF_WRITABLE {
            return Err(PyBufferError::new_err("buffer is read-only"));
        }

        if slf.inner.is_none() {
            return Err(PyBufferError::new_err("view is empty"));
        }

        let inner = slf.inner.as_ref().unwrap();

        // Compile-time guard: if Alignment's size changes (e.g. fields are
        // added or repr changes), this will fail to compile rather than
        // silently producing a corrupt memoryview.
        const _: () = assert!(size_of::<Alignment>() == 24);

        unsafe {
            (*view).buf = inner.alignments().as_ptr() as *mut c_void;
            (*view).len = (inner.len() * size_of::<Alignment>()) as isize;
            (*view).readonly = 1;
            (*view).itemsize = size_of::<Alignment>() as isize;

            // TODO: 7 bytes of padding between cigarx.bits and id is wasteful.
            // Consider flattening Cigarx<u64> fields directly into Alignment
            // (storage: u64, bits: u8, id: u32, offset: u8, strand: u8) and
            // re-packing with #[repr(C, packed)] or reordering fields to
            // reduce padding to 3 bytes (total 16 bytes).
            (*view).format = c"QBxxxxxxxIBBxx".as_ptr() as *mut _;

            (*view).ndim = 1;
            (*view).shape = &mut (*view).len;
            (*view).strides = &mut (*view).itemsize;
            (*view).suboffsets = ptr::null_mut();
            (*view).internal = ptr::null_mut();

            // Must be set last — consumes slf
            (*view).obj = slf.into_ptr();
        }

        Ok(())
    }
}

/// A Python-visible, zero-copy view over a [`SequenceRingBatch`].
///
/// The symmetric counterpart to [`AlignmentBatchView`] for the input side of
/// the pipeline. Currently unimplemented — the buffer protocol stub is
/// present to reserve the Python API surface.
///
/// ## Memory layout
///
/// Each element is a raw [`Iupac`][crate::sequence::iupac::Iupac] byte.
/// The exact format string will be `"B"` (unsigned char) once implemented.
#[pyclass]
pub struct SequenceBatchView {
    /// The underlying ring buffer slot. `None` when the view is empty.
    inner: Option<SequenceRingBatch>,
}

impl SequenceBatchView {
    /// Wrap a sequence ring buffer slot in a view.
    pub fn new(batch: SequenceRingBatch) -> Self {
        Self { inner: Some(batch) }
    }

    /// Create an empty sentinel view.
    pub fn empty() -> Self {
        Self { inner: None }
    }
}

impl Drop for SequenceBatchView {
    /// Returns the ring buffer slot to the pool.
    fn drop(&mut self) {
        if let Some(inner) = self.inner.take() {
            trace!("Dropping python sequence view");
            drop(inner);
        }
    }
}

#[pymethods]
impl SequenceBatchView {
    /// Expose the sequence batch as a read-only `memoryview`.
    ///
    /// # Not yet implemented
    ///
    /// Will expose a 1-D array of `u8` IUPAC bytes with format `"B"` once
    /// the sequence ring batch layout is finalised.
    unsafe fn __getbuffer__(
        _slf: PyRef<'_, Self>,
        _view: *mut Py_buffer,
        _flags: c_int,
    ) -> PyResult<()> {
        unimplemented!()
    }
}
