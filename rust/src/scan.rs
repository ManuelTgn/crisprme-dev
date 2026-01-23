use crate::pam::ParsedPAM;
use crate::iupac::{Iupac, matches_iupac};
use crate::target::Target;
use crate::threadpool;

use rayon::prelude::*;
use std::result::Result; 

// Type alias for the complex value data, just for cleaner code

/// Scans a sequence in parallel to find all candidate targets (k-mers) that match the PAM requirements.
///
/// The sequence is first converted to IUPAC bitmasks. The scanning is parallelized across
/// multiple threads using Rayon, with chunks extended to prevent boundary-crossing k-mers from being missed.
///
/// # Arguments
/// * `sequence` - The DNA/RNA sequence string to scan.
/// * `contig` - The identifier of the sequence (e.g., "chr1").
/// * `pam` - The parsed PAM structure containing the forward and reverse complement bitmasks.
/// * `size` - The length of the target k-mer (protospacer + pam + offset length).
/// * `right` - If `true`, the PAM is expected to be *right* of the k-mer in the scan window.
/// * `threads` - The number of threads to use for parallel processing.
///
/// # Returns
/// A `HashedTargets` object containing all unique targets and their occurrences, 
/// with grouping performed in Rust for maximum efficiency.
///
/// # Panics
/// Panics if the input sequence contains an unknown nucleotide character or if the Rayon 
/// thread pool cannot be built.
pub fn scan_targets(
    sequence: &str,
    contig: &str,
    pam: &ParsedPAM,
    size: usize,
    right: bool,
    threads: usize,
) -> Result<Vec<Target>, String> {
    // 1) Convert sequence to IUPAC bitmasks (return error instead of panic)
    let seq_bitmask: Vec<u8> = sequence
        .as_bytes()
        .iter()
        .enumerate()
        .map(|(i, &b)| {
            Iupac::from_ascii(b)
                .map(|iupac| iupac.0)
                .map_err(|e| format!("Invalid base at position {i}: {e}"))
        })
        .collect::<Result<_, _>>()?;

    let slen = seq_bitmask.len();
    if slen == 0 || size == 0 || size > slen {
        return Ok(Vec::new());
    }

    let pat = &pam.bytes;     // forward PAM bitmasks
    let rev = &pam.revcomp;   // reverse-complement PAM bitmasks
    let plen = pat.len();

    if plen == 0 || plen > size {
        return Err(format!("Invalid PAM length: plen={plen}, size={size}"));
    }

    // 2) Run the parallel scan inside a cached pool for this `threads` value
    threadpool::with_pool(threads, || {
        // chunk size (ceiling division)
        let chunk_size = (slen + threads - 1) / threads;

        (0..threads)
            .into_par_iter()
            .filter_map(|chunk_idx| {
                let orig_start = chunk_idx * chunk_size;
                if orig_start >= slen {
                    return None;
                }
                let orig_end = std::cmp::min(orig_start + chunk_size, slen);

                // extend by (size - 1) so k-mers starting within [orig_start, orig_end) are complete
                let extended_start = orig_start;
                let extended_end = std::cmp::min(orig_end + (size - 1), slen);

                let extended_chunk = &seq_bitmask[extended_start..extended_end];
                let chunk_len = extended_chunk.len();

                // Heuristic initial capacity; you can tune later
                let mut chunk_targets: Vec<Target> = Vec::new();

                if chunk_len >= size {
                    for i in 0..=(chunk_len - size) {
                        let global_pos = extended_start + i;

                        // avoid duplicates across chunks
                        if global_pos < orig_start || global_pos >= orig_end {
                            continue;
                        }

                        // get k-mer slice
                        let target_bitmask = unsafe { extended_chunk.get_unchecked(i..i + size) };

                        // skip if contains 'N' (any base)
                        if target_bitmask.iter().any(|&b| b == 0b1111) {
                            continue;
                        }

                        // NOTE: keep PAM slicing behavior
                        let (pam_slice_fwd, pam_slice_rev) = if right {
                            (&target_bitmask[0..plen], &target_bitmask[size - plen..size])
                        } else {
                            (&target_bitmask[size - plen..size], &target_bitmask[0..plen])
                        };

                        if matches_pattern(pam_slice_fwd, pat) {
                            chunk_targets.push(Target::new(
                                contig,
                                global_pos,
                                true,
                                target_bitmask.to_vec(),
                            ));
                        }

                        if matches_pattern(pam_slice_rev, rev) {
                            chunk_targets.push(Target::new(
                                contig,
                                global_pos,
                                false,
                                target_bitmask.to_vec(),
                            ));
                        }
                    }
                }

                Some(chunk_targets)
            })
            .flatten()
            .collect::<Vec<Target>>()
    })
}


/// Checks if a sequence bitmask slice matches a PAM pattern bitmask slice.
///
/// This uses the IUPAC matching logic (`matches_iupac`) to ensure every position 
/// in the sequence slice overlaps with the required pattern bits.
///
/// # Arguments
/// * `seq` - The sequence fragment bitmask slice (e.g., the PAM part of the k-mer window).
/// * `pam` - The PAM pattern bitmask slice (forward or reverse complement).
///
/// # Returns
/// * `true` if the sequence matches the pattern at all corresponding positions
#[inline(always)]
fn matches_pattern(seq: &[u8], pam: &[u8]) -> bool {
    seq.iter()
        .zip(pam.iter())
        .all(|(&a, &b)| matches_iupac(a, b))
}

// --------------------------------------------------------------------------------------------------

/// Extracts the PAM-defining slice from the k-mer window based on its relative position.
///
/// # Arguments
/// * `target` - The bitmask slice representing the k-mer *plus* the PAM window.
/// * `plen` - The length of the PAM sequence.
/// * `k` - The length of the target/protospacer (used to calculate the start position).
/// * `right` - If `true`, the PAM is on the left (at the start of the slice); 
///             if `false`, the PAM is on the right (at the end of the slice).
///
/// # Returns
/// A slice (`&[u8]`) containing only the bitmasks for the PAM region.
fn get_pam_slice(target: &[u8], plen: usize, k: usize, right: bool) -> &[u8] {
    if right {
        // if the PAM is specified to be *right* of the k-mer, then we assume the current
        // slice is structured as: [PAM_seq | K-mer]. The PAM is at the beginning.
        &target[0..plen]
    } else {
        // if the PAM is specified to be *left* of the k-mer, then the slice is 
        // structured as: [K-mer | PAM_seq]. The PAM is at the end.
        &target[k - plen..k]
    }
}