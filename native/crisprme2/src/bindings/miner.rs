//! CUDA-backed alignment ("mining") bindings.
//!
//! This module bridges Rust batch buffers to the CUDA implementation via `cxx`.
//!
//! # Call protocol
//! 1. [`initialize`] once at program start (per device).
//! 2. For each guide / parameterization:
//!    - [`pre_mine`] before processing batches.
//! 3. Repeatedly call [`mine`] for each [`SequenceRingBatch`] to fill an [`AlignmentRingBatch`].
//!    - [`mine`] returns `finish = true` when no more work remains for the current stream.
//! 4. [`post_mine`] once after finishing the batches for a guide.
//! 5. [`shutdown`] once at program end (per device).
//!
//! # Safety
//! Although the public APIs are safe, they assume:
//! - `SequenceRingBatch::gpu_ptr()` points to a valid GPU buffer for `len()` bytes
//! - `AlignmentRingBatch::gpu_ptr_mut()` points to a valid GPU buffer for `capacity()` bytes
//! - `Guide` bytes match the CUDA side layout/encoding

use crate::{
    crispr::guide::Guide,
    alignment::thresholds::Thresholds,
    memory::batch::{AlignmentRingBatch, SequenceRingBatch},
};


#[cxx::bridge(namespace = "cuda::miner")]
mod ffi {
    /// Output of the mining operation.
    struct MinerOutput {
        /// Number of alignments written into the provided output buffer.
        alignments_count: usize,
        /// Whether the miner has completed processing for the current stream/batch sequence.
        finish: bool,
    }

    unsafe extern "C++" {
        include!("api.cuh");

        fn initialize(device: u32);

        unsafe fn pre_mine(
            guide: *const u8,
            glen: u32,
            slen: u32,
            ggap: u32,
            sgap: u32,
            mism: u32,
            strand: u8,
        );

        unsafe fn mine(
            batch: *const u8,
            batch_size: u32,
            alignments: *mut u8,
            capacity: u32,
        ) -> MinerOutput;

        fn post_mine();
        fn shutdown(device: u32);
    }
}


/// Initialize CUDA miner state for a given device.
pub fn initialize(device: u32) {
    ffi::initialize(device);
}


/// Configure miner parameters for the next sequence batches.
///
/// `seq_len` is the per-target sequence length in the batch buffer.
pub fn pre_mine(guide: &Guide, seq_len: usize, thresholds: &Thresholds, strand: u8) {
    // If these can ever exceed u32, fail loudly.
    assert!(guide.len() <= u32::MAX as usize);
    assert!(seq_len <= u32::MAX as usize);

    unsafe {
        ffi::pre_mine(
            guide.as_ptr() as *const u8,
            guide.len() as u32,
            seq_len as u32,
            thresholds.qgap,
            thresholds.tgap,
            thresholds.mism,
            strand,
        );
    }
}


/// Mine a batch of sequences and write alignments into `alignments`.
///
/// Returns `true` when the CUDA miner signals completion for the current stream.
pub fn mine(batch: &SequenceRingBatch, alignments: &mut AlignmentRingBatch) -> bool {
    assert!(batch.len() <= u32::MAX as usize);
    assert!(alignments.capacity() <= u32::MAX as usize);

    let output = unsafe {
        ffi::mine(
            batch.gpu_ptr(),
            batch.len() as u32,
            alignments.gpu_ptr_mut(),
            alignments.capacity() as u32,
        )
    };

    alignments.set_len(output.alignments_count);
    output.finish
}


/// Finalize miner state after finishing all batches for the current guide.
pub fn post_mine() {
    ffi::post_mine();
}


/// Shutdown CUDA miner state for a given device.
pub fn shutdown(device: u32) {
    ffi::shutdown(device);
}
