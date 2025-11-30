use crate::{
    common::guide::Guide,
    memory::batch::{AlignmentRingBatch, SequenceRingBatch},
    common::thresholds::Thresholds
};

#[cxx::bridge(namespace = "cuda::miner")]
mod ffi {

    /// Output of the mining operation
    struct MinerOutput {
        alignments_count: usize,
        finish: bool,
    }

    unsafe extern "C++" {
        include!("crisprme-core/include/api.cuh");

        /// Invoked at the beginning of the program
        fn initialize(device: u32);

        /// Invoked before a new batch is mined
        unsafe fn pre_mine(guide: *const u8, glen: u32, slen: u32, ggap: u32, sgap: u32, mism: u32, strand: u8);

        /// Mines a sequence batch and generates a single alignment batch
        unsafe fn mine(
            batch: *const u8,
            batch_size: u32,
            alignments: *mut u8,
            capacity: u32,
        ) -> MinerOutput;

        /// Invoked after a batch has been mined
        fn post_mine();

        /// Invoked at the end of the program
        fn shutdown(device: u32);
    }
}

/// Invoked at the beginning of the program
pub fn initialize(device: u32) {
    ffi::initialize(device);
}

/// Invoked before a new batch is mined
pub fn pre_mine(guide: &Guide, seq_len: usize, thresholds: &Thresholds, strand: u8) {
    unsafe {
        ffi::pre_mine(
            guide.as_ptr() as *const u8,
            guide.len() as u32,
            seq_len as u32,
            thresholds.qgap,
            thresholds.tgap,
            thresholds.mism,
            strand
        );
    }
}

/// Mines a sequence batch and generates a single alignment batch
pub fn mine(batch: &SequenceRingBatch, alignments: &mut AlignmentRingBatch) -> bool {
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

/// Invoked after a batch has been mined
pub fn post_mine() {
    ffi::post_mine();
}

/// Invoked at the end of the program
pub fn shutdown(device: u32) {
    ffi::shutdown(device);
}
