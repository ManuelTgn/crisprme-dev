use std::ffi::{c_int, c_void};
use std::ptr;
use pyo3::{ffi, pyclass, pymethods, PyRef, PyResult};
use pyo3::exceptions::{PyBufferError, PyValueError};
use pyo3::ffi::Py_buffer;
use tracing::trace;
use crate::alignment::alignment::Alignment;
use crate::memory::batch::{AlignmentRingBatch, SequenceRingBatch};

/// Wrapper for an alignment batch
#[pyclass]
pub struct AlignmentBatchView {
    inner: Option<AlignmentRingBatch>
}

impl AlignmentBatchView {
    pub fn new(batch: AlignmentRingBatch) -> Self {
        Self { inner: Some(batch) }
    }

    pub fn empty() -> Self {
        Self { inner: None }
    }
}

impl Drop for AlignmentBatchView {
    fn drop(&mut self) {
        if let Some(inner) = self.inner.take() {
            trace!("Dropping python alignment-view");
            drop(inner);
        }
    }
}

#[pymethods]
impl AlignmentBatchView {

    /// Get id of the source batcher
    fn batcher_id(&self) -> PyResult<usize> {
        match &self.inner {
            None => Err(PyValueError::new_err("View is empty")),
            Some(inner) => Ok(inner.id())
        }
    }

    /// Returns true if the view is attached to an alignment batch
    fn is_valid(&self) -> PyResult<bool> {
        Ok(self.inner.is_some())
    }

    /// Returns a memoryview of the underling alignment batch
    unsafe fn __getbuffer__(slf: PyRef<'_, Self>, view: *mut Py_buffer, flags: c_int) -> PyResult<()> {

        if view.is_null() {
            return Err(PyBufferError::new_err("view is null"));
        }

        if (flags & ffi::PyBUF_WRITABLE) == ffi::PyBUF_WRITABLE {
            return Err(PyBufferError::new_err("buffer is read-only"));
        }

        if (slf.inner.is_none()) {
            return Err(PyBufferError::new_err("view is empty"));
        }

        let inner = slf.inner.as_ref().unwrap();
        const _: () = assert!(size_of::<Alignment>() == 24);
        unsafe {
            (*view).buf = inner.alignments().as_ptr() as *mut c_void;
            (*view).len = (inner.len() * size_of::<Alignment>()) as isize;
            (*view).readonly = 1;
            (*view).itemsize = size_of::<Alignment>() as isize;

            // TODO: 7 padding?! We have to do something about it, maybe incorporate everything inside Alignment
            // "QBxxxxxxxIBBxx": u64 + u8 + 7 pad + u32 + u8 + u8 + 2 pad = 24 bytes
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

#[pyclass]
pub struct SequenceBatchView {
    inner: Option<SequenceRingBatch>
}

impl SequenceBatchView {
    pub fn new(batch: SequenceRingBatch) -> Self {
        Self { inner: Some(batch) }
    }

    pub fn empty() -> Self {
        Self { inner: None }
    }
}

impl Drop for SequenceBatchView {
    fn drop(&mut self) {
        if let Some(inner) = self.inner.take() {
            trace!("Dropping python sequence view");
            drop(inner);
        }
    }
}

#[pymethods]
impl SequenceBatchView {

    /// Returns a memoryview of the underling alignment batch
    unsafe fn __getbuffer__(slf: PyRef<'_, Self>, view: *mut Py_buffer, flags: c_int) -> PyResult<()> {
        unimplemented!()
    }
}