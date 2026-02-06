//! CUDA FFI helpers.
//!
//! This module provides thin wrappers around CUDA-side allocation and memory copy
//! routines exposed via C++ (`api.cuh`) and bridged through `cxx`.
//!
//! # Safety
//! The raw-pointer APIs assume:
//! - `cpu` and `gpu` pointers are valid for `len` elements of `T`
//! - `T` is plain-old-data (POD) compatible with the C++ side (no drop glue, no references)
//! - the active CUDA device/context matches the allocation/copy calls
//! - `len * size_of::<T>()` does not overflow and fits in `u64`

use std::{mem, ptr::NonNull, time::Instant};
use tracing::{debug};

#[cxx::bridge(namespace = "cuda")]
mod ffi {
    unsafe extern "C++" {
        include!("crisprme-core/include/api.cuh");
        
        unsafe fn malloc(bytes: u64) -> *mut u8;
        unsafe fn free(memory: *mut u8);
        unsafe fn memcpy_to_gpu(gpu: *mut u8, cpu: *const u8, bytes: u64);
        unsafe fn memcpy_to_cpu(gpu: *const u8, cpu: *mut u8, bytes: u64);
    }
}


/// Allocate `len` elements of `T` on the GPU
/// 
/// Returns a non-null pointer on success, or None if allocation failed.
pub fn malloc<T>(len: usize) -> *mut T {
    let bytes = mem::size_of::<T>()
        .checked_mul(len)?;
    let ptr = unsafe { ffi::malloc(bytes as u64) };
    ptr as *mut T
}


/// Free GPU memory previously allocated via [`malloc`].
///
/// # Safety
/// - `memory` must have been returned by [`malloc`] (same allocator / device).
pub unsafe fn free<T>(memory: NonNull<T>) {
        ffi::free(memory.as_ptr() as *mut u8);
}


/// Copy `len` elements from CPU to GPU.
///
/// # Safety
/// - `cpu` must be valid for reads of `len` elements
/// - `gpu` must be valid for writes of `len` elements
pub unsafe fn memcpy_to_gpu<T>(cpu: *const T, gpu: *mut T, len: usize) {
    let bytes = mem::size_of::<T>()
        .checked_mul(len)
        .expect("memcpy_to_gpu: size overflow");
    let now = Instant::now();
    ffi::memcpy_to_gpu(gpu as *mut u8, cpu as *const u8, bytes as u64);
    debug!("memcpy CPU -> GPU [{}ms]", now.elapsed().as_millis());
}


/// Copy `len` elements from GPU to CPU.
///
/// # Safety
/// - `gpu` must be valid for reads of `len` elements
/// - `cpu` must be valid for writes of `len` elements
pub unsafe fn memcpy_to_cpu<T>(cpu: *mut T, gpu: *const T, len: usize) {
    let bytes = mem::size_of::<T>()
        .checked_mul(len)
        .expect("memcpy_to_cpu: size overflow");
    let now = Instant::now();
    ffi::memcpy_to_cpu(gpu as *const u8, cpu as *mut u8, bytes as u64);
    debug!("memcpy CPU <- GPU [{}ms]", now.elapsed().as_millis());
}




