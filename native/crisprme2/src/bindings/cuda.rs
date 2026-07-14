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

use std::time::Instant;

use tracing::info;

#[cxx::bridge(namespace = "cuda")]
mod ffi {
    unsafe extern "C++" {
        include!("api.cuh");

        unsafe fn malloc(bytes: u64) -> *mut u8;
        unsafe fn free(memory: *mut u8);
        unsafe fn memcpy_to_gpu(gpu: *mut u8, cpu: *const u8, bytes: u64);
        unsafe fn memcpy_to_cpu(gpu: *const u8, cpu: *mut u8, bytes: u64);
        unsafe fn pin(ptry: *const u8, bytes: u64);
        unsafe fn unpin(ptry: *const u8);
    }
}

pub fn malloc<T>(len: usize) -> *mut T {
    let bytes = std::mem::size_of::<T>() * len;
    let ptr = unsafe { ffi::malloc(bytes as u64) };
    ptr as *mut T
}

pub fn free<T>(memory: *mut T) {
    unsafe {
        ffi::free(memory as *mut u8);
    }
}

#[tracing::instrument(name = "cuda", skip_all)]
pub fn memcpy_to_gpu<T>(cpu: *const T, gpu: *mut T, len: usize) {
    let now = Instant::now();
    let bytes = std::mem::size_of::<T>() * len;
    unsafe {
        ffi::memcpy_to_gpu(gpu as *mut u8, cpu as *const u8, bytes as u64);
    }
    tracing::trace!("memcpy CPU -> GPU [{}ms]", now.elapsed().as_millis());
}

#[tracing::instrument(name = "cuda", skip_all)]
pub fn memcpy_to_cpu<T>(cpu: *mut T, gpu: *const T, len: usize) {
    let now = Instant::now();
    let bytes = std::mem::size_of::<T>() * len;
    unsafe {
        ffi::memcpy_to_cpu(gpu as *const u8, cpu as *mut u8, bytes as u64);
    }
    tracing::trace!("memcpy CPU <- GPU [{}ms]", now.elapsed().as_millis());
}

#[tracing::instrument(name = "cuda")]
pub fn pin(memory: *const u8, bytes: usize) {
    unsafe {
        ffi::pin(memory, bytes as u64);
    }
}

#[tracing::instrument(name = "cuda")]
pub fn unpin(memory: *const u8) {
    unsafe {
        ffi::unpin(memory);
    }
}
