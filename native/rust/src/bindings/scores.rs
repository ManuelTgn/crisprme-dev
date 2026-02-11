//! GPU scoring kernels.
//!
//! Exposes CUDA scoring routines over IUPAC-encoded sequences.
//!
//! Layout assumptions:
//! - `strings` contains `n` strings, each of length `slen`, concatenated back-to-back.
//! - `result` has length at least `n` (one score per string), unless the CUDA API specifies otherwise.

use crate::sequence::iupac::Iupac;


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
    }
}


/// Compute scores for `n` strings against a query.
///
/// # Panics
/// - if buffer lengths are inconsistent
///
/// # Safety assumptions
/// `Iupac` must be a byte-sized POD type compatible with the CUDA kernel
/// (ideally `#[repr(u8)]`).
pub fn scores_into(
    query: &[Iupac],
    strings: &[Iupac],
    result: &mut [u8],
    slen: usize,
    n: usize,
) {
    assert!(query.len() <= slen);
    assert!(slen <= i32::MAX as usize);
    assert!(n <= i32::MAX as usize);

    let needed_strings = slen
        .checked_mul(n)
        .expect("scores_into: slen*n overflow");
    assert!(
        strings.len() >= needed_strings,
        "scores_into: strings too short: got {}, need {} (= slen*n)",
        strings.len(),
        needed_strings
    );
    assert!(
        result.len() >= n,
        "scores_into: result too short: got {}, need {} (= n)",
        result.len(),
        n
    );

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
