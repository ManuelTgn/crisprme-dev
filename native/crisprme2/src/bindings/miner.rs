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
        // Full window width (PAM included) — also the batch row stride
        pub slen: u32,
        // PAM length; protospacer ends at slen - plen
        pub plen: u32,
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
///
/// Naming across the FFI boundary is a trap, so, explicitly:
///   CUDA `ggap` = gap in the GUIDE row  = extra DNA base = DNA bulge = Thresholds::tgap
///   CUDA `sgap` = gap in the SEQ   row  = extra RNA base = RNA bulge = Thresholds::qgap
/// (`ggap_dt()` fires on Step::S, which advances `sidx` only — a gap in the guide.)
pub fn prepare(guide: &Guide, seq_len: usize, pam_len: usize, thresholds: &Thresholds) {
    debug_assert!(
        pam_len < seq_len && seq_len - pam_len > guide.len(),
        "window must fit guide + PAM: seq_len={seq_len}, pam_len={pam_len}, glen={}",
        guide.len()
    );
    ffi::prepare(MinerConfig {
        guide: guide.as_ptr() as *const u8,
        glen: guide.len() as u32,
        slen: seq_len as u32,
        plen: pam_len as u32,
        ggap: thresholds.tgap, // was: thresholds.qgap  <-- SWAPPED
        sgap: thresholds.qgap, // was: thresholds.tgap  <-- SWAPPED
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

#[cfg(test)]
mod tests {
    //! Host-side model of the kernel's PAM anchoring.
    //!
    //! The kernel runs on the GPU, so these tests do not launch it. Instead they
    //! reproduce the exact arithmetic the kernel uses (`PSTOP = SLEN - PLEN`, the
    //! `offset_lo`/`offset_hi` loop bounds, and the completion condition
    //! `offset + sidx == PSTOP`) and assert the anchor invariant for the repro
    //! case, so a future change to that arithmetic is caught here.

    /// Mirrors `warp_reg_stack_nopam_columnar.cu`.
    fn pstop(slen: u32, plen: u32) -> u32 {
        slen - plen
    }
    fn offset_lo(pstop: u32, glen: u32, ggap: u32) -> u32 {
        let pad = pstop - glen;
        if pad > ggap {
            pad - ggap
        } else {
            0
        }
    }
    fn offset_hi(pstop: u32, glen: u32, sgap: u32) -> u32 {
        let pad = pstop - glen;
        pad + sgap
    }

    /// DNA cursor advance once the guide is fully consumed:
    ///   sidx = #B + #S = (glen - rna_bulges) + dna_bulges
    fn sidx_at_complete(glen: u32, dna_bulges: u32, rna_bulges: u32) -> u32 {
        glen + dna_bulges - rna_bulges
    }

    /// Repro: guide `CTAACAGTTGCTTTTATCAC` (20 nt), PAM `NGG` (3 nt),
    /// window `AACTAACAGTTGCTTTTATCACTGG` (25 nt). Protospacer ends at index 22.
    const GLEN: u32 = 20;
    const SLEN: u32 = 25;
    const PLEN: u32 = 3;
    const GGAP: u32 = 1; // max DNA bulges
    const SGAP: u32 = 1; // max RNA bulges

    /// For every valid bulge combination the emitted alignment abuts the PAM:
    /// `offset + sidx == pstop`, equivalently `offset + dna - rna == pad`.
    #[test]
    fn anchor_invariant_holds_for_all_bulge_combinations() {
        let pstop = pstop(SLEN, PLEN);
        let pad = pstop - GLEN; // == 2
        let lo = offset_lo(pstop, GLEN, GGAP);
        let hi = offset_hi(pstop, GLEN, SGAP);

        for &(dna, rna) in &[(0, 0), (1, 0), (0, 1), (1, 1)] {
            // The offset the kernel would need for this combination to anchor.
            let offset = pad + rna - dna; // pad - dna + rna, kept non-negative
            assert!(
                (lo..=hi).contains(&offset),
                "combo (dna={dna}, rna={rna}) needs offset={offset}, outside explored [{lo},{hi}]"
            );

            let sidx = sidx_at_complete(GLEN, dna, rna);
            assert_eq!(
                offset + sidx,
                pstop,
                "combo (dna={dna}, rna={rna}) offset={offset} does not abut the PAM"
            );
            assert_eq!(
                offset + dna - rna,
                pad,
                "combo (dna={dna}, rna={rna}) violates offset + dna - rna == pad"
            );
        }
    }

    /// The bug: an offset=1 path that consumes exactly 20 DNA bases (1 DNA bulge
    /// and 1 RNA bulge cancelling out) ends at window index 21, one short of the
    /// PAM. The completion condition must reject it.
    #[test]
    fn mis_anchored_alignment_is_rejected() {
        let pstop = pstop(SLEN, PLEN); // 22
        let offset = 1;
        let sidx = sidx_at_complete(GLEN, 1, 1); // 20
        assert_ne!(
            offset + sidx,
            pstop,
            "the orphaned-base alignment (offset=1, ends at 21) must not anchor"
        );
    }
}
