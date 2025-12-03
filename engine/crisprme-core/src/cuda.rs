use std::ptr;
use std::time::{Duration, Instant};
use tracing::{info};
use crate::{
    common::{guide::Guide, thresholds::Thresholds},
    memory::{batch::{AlignmentBatch, AlignmentRingBatch, SequenceBatch, SequenceRingBatch}, ring::RingAdapter},
    utils::IUPAC,
};

#[cxx::bridge]
mod ffi {

    unsafe extern "C++" {
        include!("crisprme-core/include/scores.cuh");
        unsafe fn scores(
            query: *const u8,
            strings: *const u8,
            result: *mut u8,
            qlen: i32,
            slen: i32,
            n: i32,
        );

        unsafe fn cuda_malloc(bytes: u64) -> *mut u8;
        unsafe fn cuda_free(memory: *mut u8);

        unsafe fn cuda_memcpy_to_gpu(gpu: *mut u8, cpu: *const u8, bytes: u64);
        unsafe fn cuda_memcpy_to_cpu(gpu: *const u8, cpu: *mut u8, bytes: u64);

        unsafe fn cuda_mine_prepare(
            guide: *const u8,
            qlen: u32,
            slen: u32,
            ggap: u32,
            sgap: u32,
            mism: u32,
            alignments: u32,
        );
        unsafe fn cuda_mine_next(
            batch: *const u8,
            batch_offset: u32,
            batch_size: u32,
            result: *mut u8,
            result_size: u32,
        ) -> u32;

    }
}

pub fn cuda_mine_prepare(guide: &Guide, slen: usize, alignments: usize, thresholds: &Thresholds) {
    unsafe {
        ffi::cuda_mine_prepare(
            guide.as_ptr() as *const u8,
            guide.len() as u32,
            slen as u32,
            thresholds.qgap,
            thresholds.tgap,
            thresholds.mism,
            alignments as u32,
        );
    }
}

pub fn cuda_mine_next(batch: &SequenceRingBatch, batch_offset: usize, result: &mut AlignmentRingBatch) {
    let mined = unsafe {
        ffi::cuda_mine_next(
            batch.gpu_ptr(),
            batch_offset as u32,
            batch.len() as u32,
            result.gpu_ptr_mut(),
            result.capacity() as u32
        )
    };

    result.set_len(mined as usize);
}

pub fn scores(query: &[IUPAC], strings: &[IUPAC], slen: usize, n: usize) -> Vec<u8> {
    assert!(query.len() <= slen);

    let mut result = vec![255; n];
    unsafe {
        ffi::scores(
            query.as_ptr() as *const u8,
            strings.as_ptr() as *const u8,
            result.as_mut_ptr(),
            query.len() as i32,
            slen as i32,
            n as i32,
        );
    }

    result
}

pub fn scores_with_arena(
    query: &[crate::common::iupac::Iupac],
    strings: &[crate::common::iupac::Iupac],
    result: &mut [u8],
    slen: usize,
    n: usize,
) {
    assert!(query.len() <= slen);
    unsafe {
        ffi::scores(
            query.as_ptr() as *const u8,
            strings.as_ptr() as *const u8,
            result.as_mut_ptr(),
            query.len() as i32,
            slen as i32,
            n as i32,
        );
    }
}

pub fn malloc<T>(len: usize) -> *mut T {
    let bytes = std::mem::size_of::<T>() * len;
    let ptr = unsafe { ffi::cuda_malloc(bytes as u64) };
    ptr as *mut T
}

pub fn free<T>(memory: *mut T) {
    unsafe {
        ffi::cuda_free(memory as *mut u8);
    }
}

#[tracing::instrument(name = "cuda", skip_all)]
pub fn memcpy_to_gpu<T>(cpu: *const T, gpu: *mut T, len: usize) {
    let now = Instant::now();
    let bytes = std::mem::size_of::<T>() * len;
    unsafe {
        ffi::cuda_memcpy_to_gpu(gpu as *mut u8, cpu as *const u8, bytes as u64);
    }
    info!("memcpy CPU -> GPU [{} ms]", 
        now.elapsed().as_millis());
}

#[tracing::instrument(name = "cuda", skip_all)]
pub fn memcpy_to_cpu<T>(cpu: *mut T, gpu: *const T, len: usize) {
    let now = Instant::now();
    let bytes = std::mem::size_of::<T>() * len;
    unsafe {
        ffi::cuda_memcpy_to_cpu(gpu as *const u8, cpu as *mut u8, bytes as u64);
    }
    info!("memcpy CPU <- GPU [{} ms]", 
        now.elapsed().as_millis());
}
