use rayon::prelude::*;
use std::result::Result;

use crate::crispr::pam::{build_sparse, PAM};
use crate::sequence::iupac::{matches_iupac, sequence_encoder};
use crate::utils::threadpool;

/// IUPAC bitmask for `N` (any base).
///
/// In this representation, `N` is `0b1111` (A|C|G|T). It is used both as a
/// degenerate PAM code and as an ambiguity marker in the input sequence.
/// In the scanner, k-mers containing `N` are skipped (see `scan_targets`).
const N_MASK: u8 = 0b1111;

// pub fn scan_targets(
//     sequence: &str,
//     pam: &PAM,
//     size: usize,
//     upstream: bool,
//     threads: usize,
// ) -> Result<Vec<usize>, String> {
//     // encode contig chunk sequence in bits
//     let seq_bitmask = sequence_encoder(sequence);

//     // extract target candidates from contig chunk
//     scan_targets_bitmask(&seq_bitmask, pam, size, upstream, threads)
// }

pub fn scan_targets_bitmask(
    seq_bitmask: &[u8],
    pam: &PAM,
    size: usize,
    upstream: bool,
    threads: usize,
) -> Result<Vec<usize>, String> {
    // get sequence length
    let slen = seq_bitmask.len();
    if slen == 0 || size == 0 || size > slen {
        return Ok((Vec::new()));
    }

    // get encoded pam
    let pat = if upstream { &pam.revcomp } else { &pam.bytes };
    let plen = pat.len();
    if plen == 0 || plen > size {
        return Err(format!("Invalid PAM length: plen={plen}, size={size}"));
    }

    // decide whether perform sparse pam matching (sparse pam example NNNRGG)
    let (idx, mask) = build_sparse(pat);
    let use_sparse = !pam.unconstrained && idx.len() < plen;

    // compute start positions for pam
    let pam_start = size - plen;

    threadpool::with_pool(threads, || {
        // define chunk sizes for threads spawning over sequences
        let chunk_size = (slen + threads - 1) / threads;

        let per_chunk: Vec<Vec<usize>> = (0..threads)
            .into_par_iter()
            .filter_map(|chunk_idx| {
                // compute start/stop positions for current chunk
                let orig_start = chunk_idx * chunk_size;
                if orig_start >= slen {
                    return None;
                }
                let orig_stop = std::cmp::min(orig_start + chunk_size, slen);

                // retrieve chunk
                let extended_start = orig_start;
                let extended_stop = std::cmp::min(orig_stop + (size - 1), slen);
                let extended_chunk = &seq_bitmask[extended_start..extended_stop];
                let chunk_len = extended_chunk.len();

                // define positions and strands arrays for current chunk
                let mut chunk_pos: Vec<usize> = Vec::new();

                if chunk_len >= size {
                    for i in 0..=(chunk_len - size) {
                        // global position within contig chunk
                        let global_pos = extended_start + i;

                        if global_pos < orig_start || global_pos >= orig_stop {
                            continue;
                        }

                        // retrieve target candidate
                        let target_bitmask = unsafe { extended_chunk.get_unchecked(i..i + size) };

                        if target_bitmask.iter().any(|&b| b == N_MASK) {
                            continue; // skip if candidate contains any N
                        }

                        if pam.unconstrained {
                            // degenerate PAM -> skip PAM matching
                            chunk_pos.push(global_pos);
                        } else {
                            // perform PAM matching to filter out candidates
                            let found_match = if use_sparse {
                                matches_pattern_sparse(target_bitmask, pam_start, &idx, &mask)
                            } else {
                                let pam_slice = &target_bitmask[pam_start..pam_start + plen];
                                matches_pattern(pam_slice, pat)
                            };

                            if found_match {
                                chunk_pos.push(global_pos);
                            }
                        }
                    }
                }

                Some((chunk_pos))
            })
            .collect();

        // collect results from targets extraction
        let total_hits: usize = per_chunk.iter().map(|(p)| p.len()).sum();
        let mut pos: Vec<usize> = Vec::with_capacity(total_hits);

        for (p) in per_chunk {
            pos.extend(p);
        }

        (pos)
    })
}

/// Checks if a sequence bitmask slice matches a PAM pattern bitmask slice (dense matching).
///
/// A dense match checks all `pam.len()` positions and requires that each
/// nucleotide bitmask overlaps with the corresponding PAM bitmask according
/// to IUPAC semantics (bitwise AND is non-zero).
///
/// This is typically optimal for low-degeneracy PAMs (few or no `N` positions).
///
/// # Arguments
/// * `seq` - Sequence fragment bitmask slice (PAM region within the scan window).
/// * `pam` - PAM bitmask pattern slice (forward or reverse complement).
///
/// # Returns
/// * `true` if all positions match under IUPAC semantics, otherwise `false`.
#[inline(always)]
fn matches_pattern(seq: &[u8], pam: &[u8]) -> bool {
    seq.iter()
        .zip(pam.iter())
        .all(|(&a, &b)| matches_iupac(a, b))
}

/// Checks whether a target sequence matches a PAM pattern using a sparse representation.
///
/// This function evaluates only the informative PAM positions (i.e., those not equal to `N`),
/// using IUPAC overlap matching (bitwise AND non-zero).
///
/// Compared to dense matching over the full PAM length, this approach can
/// significantly reduce the cost of PAM evaluation when the PAM contains many `N`s.
///
/// # Arguments
/// * `target_bitmask` - Full scan window encoded as IUPAC bitmasks (length `size`).
/// * `pam_start` - Start index of the PAM region within `target_bitmask`.
/// * `idx` - Indices of informative PAM positions (relative to `pam_start`).
/// * `mask` - IUPAC bitmasks corresponding to `idx`.
///
/// # Returns
/// * `true` if all informative PAM positions match, otherwise `false`.
///
/// # Safety
/// This function uses unchecked indexing for performance. Correctness relies on:
/// * `idx.len() == mask.len()`
/// * `pam_start + idx[i] < target_bitmask.len()` for all `i`
///
/// These invariants are guaranteed by construction in `scan_targets`.
#[inline(always)]
fn matches_pattern_sparse(
    target_bitmask: &[u8],
    pam_start: usize,
    idx: &[usize],
    mask: &[u8],
) -> bool {
    // idx and mask have the same length
    for t in 0..idx.len() {
        let seq_mask = unsafe { *target_bitmask.get_unchecked(pam_start + idx[t]) };
        let pam_mask = unsafe { *mask.get_unchecked(t) };
        if (seq_mask & pam_mask) == 0 {
            return false;
        }
    }

    true
}
