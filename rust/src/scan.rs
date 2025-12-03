use crate::hashing;
use crate::pam::ParsedPAM;
use crate::iupac::{Iupac, matches_iupac};
use crate::target::Target;

use std::sync::Arc;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::result::Result; 

// Type alias for the complex value data, just for cleaner code
type OccurrenceData = Vec<(String, usize, bool)>;

/// Scans a sequence in parallel to find all candidate targets (k-mers) that match the PAM requirements.
///
/// The sequence is first converted to IUPAC bitmasks. The scanning is parallelized across
/// multiple threads using Rayon, with chunks extended to prevent boundary-crossing k-mers from being missed.
///
/// # Arguments
/// * `sequence` - The DNA/RNA sequence string to scan.
/// * `contig` - The identifier of the sequence (e.g., "chr1").
/// * `pam` - The parsed PAM structure containing the forward and reverse complement bitmasks.
/// * `k` - The length of the target k-mer (protospacer length).
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
    k: usize, 
    right: bool, 
    is_first_chunk: bool,
    path: &str,
    threads: usize
// ) -> HashMap<Vec<u8>, OccurrenceData> {
) -> hashing::HashedTargets {
    // 1. Convert the entire sequence string into a single vector of IUPAC bitmasks
    let seq_bitmask: Result<Vec<u8>, String> = sequence.as_bytes()
        .iter()
        .map(|&b| Iupac::from_ascii(b).map(|iupac| iupac.0))
        .collect();

    // Unwrap the result, panicking if there's an error (e.g., unknown nucleotide), 
    // as the scanning cannot proceed with invalid input.
    let seq_bitmask: Vec<u8> = seq_bitmask
        .expect("Failed to process sequence: contains unknown nucleotide character.");

    // 2. Prepare data for parallel access.
    // let pat = &pam.bytes;  // forward PAM pattern (Arc allows shared access)
    let pat = Arc::new(pam.bytes.clone());
    // let rev = &pam.revcomp;  // reverse complement PAM pattern
    let rev = Arc::new(pam.revcomp.clone());
    let slen = seq_bitmask.len();  // total sequence length (in masks)
    let plen = pat.len();  // PAM pattern length (in masks)

    // // Calculate the maximum possible starting position for a full window
    // let max_start_pos = if slen >= k { slen - k + 1 } else { 0 };

    // // 3. Initialize the Rayon thread pool
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to create Rayon thread pool");

    // 4. Execute the parallel scanning logic
    let targets: Vec<Target> = pool.install(|| {
        // calculate chunk size (ceiling division)
        let chunk_size = (slen + threads - 1) / threads;

        // parallel iterate over chunk indices; produce Vec<Vec<(pos, String)>>
        (0..threads)
            .into_par_iter()
            .filter_map(|chunk_idx| {
                let orig_start = chunk_idx * chunk_size;
                // if the chunk start is beyond the sequence length, there's nothing left to scan
                if orig_start >= slen {
                    return None;
                }
                let orig_end = std::cmp::min(orig_start + chunk_size, slen);

                // define the extended window for the current chunk
                // the chunk must be extended by (k - 1) masks to ensure that the last
                // k-mer *starting* within the original range can be fully constructed.
                let extended_start = orig_start;
                let extended_end = std::cmp::min(orig_end + (k - 1), slen);

                // slice the bitmask data for the extended chunk
                let extended_chunk = &seq_bitmask[extended_start..extended_end];
                let chunk_len = extended_chunk.len();

                let mut chunk_targets = Vec::with_capacity(slen / 1000);
                if chunk_len >= k {
                    // iterate over all possible k-mer start positions within the extended chunk
                    for i in 0..=(chunk_len - k) {
                        let global_pos = extended_start + i;

                        // check if the k-mer's start position falls within the thread's 
                        // *original* boundaries to prevent duplicate results across threads.
                        if global_pos >= orig_start && global_pos < orig_end {

                            // retrieve the k-mer target sequence (bitmasks) from the extended chunk
                            let target_bitmask = unsafe { extended_chunk.get_unchecked(i..i + k) };

                            // skip the target if it contains the 'N' (Any base) bitmask
                            if target_bitmask.iter().any(|&b| b == 0b1111) {
                                continue;
                            }

                            let (pam_slice_fwd, pam_slice_rev) = if right {
                                (&target_bitmask[0..plen], &target_bitmask[k - plen..k])
                            } else {
                                (&target_bitmask[k - plen..k], &target_bitmask[0..plen])
                            };

                            // check for PAM match on the forward strand
                            if matches_pattern(pam_slice_fwd, &pat) {
                                chunk_targets.push(Target::new(contig, global_pos, true, target_bitmask.to_vec()));
                            }

                            // check for PAM match on the reverse strand
                            if matches_pattern(pam_slice_rev, &rev) {
                                chunk_targets.push(Target::new(contig, global_pos, false, target_bitmask.to_vec()));
                            }
                        }
                    }
                }
                Some(chunk_targets)  // return the vector of targets found in this chunk
            })
            .flatten()  // flatten the Vec<Vec<Target>> into a single Vec<Target>
            .collect()
    });

    // 5. Hash and group the results entirely in Rust for performance.
    // This is the step that replaces the slow Python loop
    let hashed_targets = hashing::hash_and_group_targets(targets);

    // --- NEW: Save the HashedTargets to indexed binary files ---
    hashed_targets.save_indexed_binary(path, is_first_chunk)
        .expect("FATAL: Failed to save targets to indexed binary files.");

    // Now return the structure to lib.rs (it's unused in lib.rs, but required for flow)
    // NOTE: Returning HashedTargets just to keep the current structure; 
    // lib.rs will simply discard it.
    hashed_targets
}

// --------------------------------------------------------------------------------------------------

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
    // for i in 0..seq.len() {
    //     if !matches_iupac(unsafe { *seq.get_unchecked(i) }, unsafe { *pam.get_unchecked(i) } ) {
    //         return false;
    //     }
    // }
    // true
    seq.iter()
        .zip(pam.iter())
        .all(|(a, b)| matches_iupac(*a, *b))
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