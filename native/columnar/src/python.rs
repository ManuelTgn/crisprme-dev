
use std::ffi::{CStr, c_int, c_void};
use pyo3::{PyResult, ffi::Py_buffer, pyclass, pymethods};
use bytemuck::Pod;

pub trait PyBufferFormat: Pod {
    /// Null-terminated PEP 3118 format string, e.g. `b"B\0"` for `u8`.
    const FORMAT: &'static CStr;
}

macro_rules! impl_format {
    ($ty:ty, $fmt:literal) => {
        impl PyBufferFormat for $ty {
            const FORMAT: &'static CStr = unsafe {
                CStr::from_bytes_with_nul_unchecked(concat!($fmt, "\0").as_bytes())
            };
        }
    };
}

impl_format!(u8,  "B");
impl_format!(u16, "H");
impl_format!(u32, "I");
impl_format!(u64, "Q");
impl_format!(i8,  "b");
impl_format!(i16, "h");
impl_format!(i32, "i");
impl_format!(i64, "q");
impl_format!(f32, "f");
impl_format!(f64, "d");

// =============================================================================
// PyChunk
// =============================================================================

/// A zero-copy view into a type T (usually a &mut [T] or &mut [[T; N]]).
///
/// Implements the PEP 3118 buffer protocol so that Python code can obtain a
/// `memoryview` (or a numpy array via `numpy.frombuffer`) backed directly by
/// the Rust-owned memory, no data is copied.
#[derive(Debug, Clone, Copy)]
#[pyclass(unsendable, skip_from_py_object)]
pub struct PyBuffer {

    /// Raw pointer to the first element of this chunk's data.
    data: *mut u8,

    /// Size of a single element in bytes.
    item_size: isize,

    /// PEP 3118 format string (null-terminated).
    format: &'static CStr,

    /// Number of dimensions: 1 for scalar columns, 2 for `[T; N]` array columns.
    ndim: c_int,

    // Buffer protocol requires stable addresses for shape/strides, so we
    // store them inline and hand out pointers to these fields.
    // For 1D: shape[0] = rows, strides[0] = item_size.
    // For 2D: shape = [rows, N], strides = [N * item_size, item_size].
    strides: [isize; 2],
    shape:   [isize; 2],
}

#[pymethods]
impl PyBuffer {

    unsafe fn __getbuffer__(&mut self, view: *mut Py_buffer, _flags: c_int) -> PyResult<()> {

        if view.is_null() {
            return Err(pyo3::exceptions::PyBufferError::new_err(
                "null Py_buffer pointer",
            ));
        }

        let view = unsafe { &mut *view };
        view.buf = self.data as *mut c_void;
        view.readonly = 0;

        view.ndim     = self.ndim;
        view.itemsize = self.item_size;
        view.format   = self.format.as_ptr() as *mut _;
        view.shape    = self.shape.as_ptr() as *mut _;
        view.strides  = self.strides.as_ptr() as *mut _;

        // Total byte length = product of shape dimensions * item_size
        view.len = match self.ndim {
            1 => self.shape[0] * self.item_size,
            _ => self.shape[0] * self.shape[1] * self.item_size,
        };

        view.suboffsets = std::ptr::null_mut();
        view.internal   = std::ptr::null_mut();

        Ok(())
    }

    unsafe fn __releasebuffer__(&self, _view: *mut Py_buffer) {
        // No-op: we don't allocate anything extra in __getbuffer__.
    }
}

impl PyBuffer {

    /// Create a PyBuffer from a slice
    pub unsafe fn from_slice<T>(value: &[T]) -> Self 
    where 
        T: PyBufferFormat
    {
        let item_size = std::mem::size_of::<T>() as isize;
        let len = value.len() as isize;
        PyBuffer {
            data: value.as_ptr() as _, 
            item_size, 
            format: T::FORMAT, 
            ndim: 1, 
            strides: [item_size, 0], 
            shape: [len, 0] 
        }
    }

    /// Create a PyBuffer from a mutable slice
    pub unsafe fn from_slice_mut<T>(value: &mut [T]) -> Self 
    where 
        T: PyBufferFormat
    {
        let item_size = std::mem::size_of::<T>() as isize;
        let len = value.len() as isize;
        PyBuffer {
            data: value.as_mut_ptr() as _, 
            item_size, 
            format: T::FORMAT, 
            ndim: 1, 
            strides: [item_size, 0], 
            shape: [len, 0] 
        }
    }

    /// Create a PyBuffer from a slice of arrays
    pub unsafe fn from_array<T, const N: usize>(value: &[[T; N]]) -> Self
    where
        T: PyBufferFormat
    {
        let item_size = std::mem::size_of::<T>() as isize;
        let len = value.len() as isize;
        let n = N as isize;
        PyBuffer {
            data: value.as_ptr() as _, 
            item_size, 
            format: T::FORMAT, 
            ndim: 2, 
            strides: [item_size * n, item_size], 
            shape: [len, n] 
        }
    }

    /// Create a PyBuffer from a mutable slice of arrays
    pub unsafe fn from_array_mut<T, const N: usize>(value: &mut [[T; N]]) -> Self
    where
        T: PyBufferFormat
    {
        let item_size = std::mem::size_of::<T>() as isize;
        let len = value.len() as isize;
        let n = N as isize;
        PyBuffer {
            data: value.as_mut_ptr() as _, 
            item_size, 
            format: T::FORMAT, 
            ndim: 2, 
            strides: [item_size * n, item_size], 
            shape: [len, n] 
        }
    }
}

#[cfg(test)]
mod test {
    use columnar_derive::Columnar;
    use crate::memory::region_min_size;
    use crate::Schema;
    use super::*;

    #[derive(Columnar)]
    struct Sequences {
        pub id: u32,
        pub sequence: [u8; 2],
    }

    #[pyclass(unsendable)]
    struct PySequences {
        pub ids: PyBuffer,
        pub sequences: PyBuffer,
    }

    #[pymethods]
    impl PySequences {

        fn ids(&self) -> PyResult<PyBuffer> { Ok(self.ids) }
        fn sequences(&self) -> PyResult<PyBuffer> { Ok(self.sequences) }
    }
}