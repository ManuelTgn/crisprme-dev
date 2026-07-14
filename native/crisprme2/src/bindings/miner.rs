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
    alignment::thresholds::Thresholds,
    bindings::miner::ffi::{MinerConfig, MinerInput},
    crispr::guide::Guide,
    memory::batch::{AlignmentRingBatch, SequenceRingBatch},
};

#[cxx::bridge(namespace = "cuda::miner")]
mod ffi {

    /// Output of the mining operation.
    struct MinerOutput {
        /// Number of alignments written into the provided output buffer.
        pub alignments_count: usize,
        /// Whether the miner has completed processing for the current stream/batch sequence.
        pub finish: bool,
    }

    /// Configuration for the miner
    struct MinerConfig {
        // Guide to use
        pub guide: *const u8,
        // Size of the guide and sequences
        pub glen: u32,
        pub slen: u32,
        // Thresholds
        pub ggap: u32,
        pub sgap: u32,
        pub mism: u32,
    }

    /// Input for the miner
    struct MinerInput {
        // Sequences to be mined, as [[Iupac; N]]
        pub sequences: *const u8,
        // Number of sequences
        pub seq_count: u32,
        // Output columns (Cigarx64, SeqRowIdx, u8)
        pub cigarx: *mut u64,
        pub index: *mut u32,
        pub offset: *mut u8,
        // Result buffer capacity
        pub capacity: u32,
    }

    unsafe extern "C++" {
        include!("api.cuh");

        fn initialize(device: u32);

        fn prepare(config: MinerConfig);

        fn launch(input: MinerInput) -> MinerOutput;

        fn post_mine();

        fn shutdown(device: u32);
    }
}

/// Initialize CUDA miner state for a given device.
pub fn initialize(device: u32) {
    ffi::initialize(device);
}

/// Configure miner parameters for the next sequence batches.
pub fn prepare(guide: &Guide, seq_len: usize, thresholds: &Thresholds) {
    ffi::prepare(MinerConfig {
        guide: guide.as_ptr() as *const u8,
        glen: guide.len() as u32,
        slen: seq_len as u32,
        ggap: thresholds.qgap,
        sgap: thresholds.tgap,
        mism: thresholds.mism,
    });
}

/// Mina a batch of sequences and write alignments into `alignments`.
pub fn mine(
    sequences: *const u8,
    seq_count: u32,
    cigarx: *mut u64,
    index: *mut u32,
    offset: *mut u8,
    capacity: u32,
) -> (bool, usize) {
    let output = ffi::launch(MinerInput {
        sequences,
        seq_count,
        cigarx,
        index,
        offset,
        capacity,
    });
    (output.finish, output.alignments_count)
}

/// Finalize miner state after finishing all batches for the current guide.
pub fn post_mine() {
    ffi::post_mine();
}

/// Shutdown CUDA miner state for a given device.
pub fn shutdown(device: u32) {
    ffi::shutdown(device);
}
