//! Parallel target candidate scanner.
//!
//! This module implements the hot-loop genome scanning logic used to extract
//! candidate target windows (k-mers) that satisfy a PAM definition.
//!
//! The scanner operates on a sequence converted into 4-bit IUPAC bitmasks,
//! enabling constant-time matching via bitwise operations. Scanning is
//! parallelized by splitting the sequence into disjoint chunks, each extended
//! by `(size - 1)` bases to ensure k-mers that cross chunk boundaries are not
//! missed. Duplicate reporting is avoided by only emitting k-mers whose start
//! position lies in the original (non-extended) chunk interval.
//!
//! Performance optimizations:
//! - Uses a cached Rayon thread pool (see `threadpool::with_pool`) to avoid pool
//!   construction overhead across repeated scans.
//! - Skips k-mers containing `N` (`0b1111`) to avoid ambiguous candidates.
//! - Implements fast-paths for fully unconstrained PAMs (`NNN...`) and partially
//!   unconstrained PAMs by using sparse matching over informative PAM positions.

use crate::pam::ParsedPAM;
use crate::iupac::{Iupac, matches_iupac};
use crate::threadpool;

use rayon::prelude::*;
use std::result::Result; 



/// IUPAC bitmask for `N` (any base).
///
/// In this representation, `N` is `0b1111` (A|C|G|T). It is used both as a
/// degenerate PAM code and as an ambiguity marker in the input sequence.
/// In the scanner, k-mers containing `N` are skipped (see `scan_targets`).
const N_MASK: u8 = 0b1111;


/// Scans a sequence in parallel and returns candidate target windows that satisfy PAM constraints.
///
/// The input `sequence` is first converted to IUPAC bitmasks using a fast lookup table
/// (`Iupac::from_ascii`). The scan then slides a fixed-size window of length `size`
/// across the sequence and evaluates PAM compatibility on both orientations:
/// - Forward orientation uses `pam.bytes`
/// - Reverse orientation uses `pam.revcomp`
///
/// Parallelization strategy:
/// - The sequence is split into `threads` chunks (ceiling division).
/// - Each chunk is *extended* by `(size - 1)` at the end to prevent missing k-mers
///   that cross chunk boundaries.
/// - Only k-mers whose start position lies within the original (non-extended) chunk
///   interval are emitted, avoiding duplicates between adjacent chunks.
///
/// PAM handling strategy:
/// - If `pam.unconstrained == true` (e.g., `NNN...`), PAM checks are skipped and
///   every k-mer (excluding those containing `N` in the sequence) is emitted.
/// - Otherwise, PAM matching is performed either:
///   - densely (check all PAM positions), or
///   - sparsely (check only informative PAM positions, skipping `N` entries in the PAM),
///     depending on the degeneracy heuristic.
/// 
/// # Arguments
/// * `sequence` - DNA/RNA sequence to scan (ASCII).
/// * `pam` - Parsed PAM containing forward and reverse-complement bitmasks plus
///           `unconstrained` flag.
/// * `size` - Scan window length (protospacer + PAM + optional bulge offset).
/// * `right` - Controls which end of the scan window is interpreted as the PAM.
///
///   **Important:** `right` reflects the scanner’s *window layout convention*.
///   With the current implementation:
///   - if `right == true`, the *forward PAM slice* starts at index `0`
///   - if `right == false`, the *forward PAM slice* starts at index `size - plen`
///
///   (and the reverse slice uses the opposite end). This preserves the behavior of
///   the existing Python-level pipeline.
/// * `threads` - Number of threads to use. Must be > 0.
///
/// # Returns
/// * `Ok(Vec<Target>)` containing all candidate targets found.
/// * `Err(String)` if the sequence contains invalid characters or PAM length is inconsistent.
///
/// # Errors
/// - Returns `Err` if `sequence` contains non-IUPAC characters (e.g., `*`, `-`, etc.).
/// - Returns `Err` if `pam.bytes.len() == 0` or `pam.bytes.len() > size`.
///
/// # Notes
/// - k-mers containing `N` in the sequence are skipped.
/// - Fully unconstrained PAMs (`NNN...`) preserve the current behavior of emitting
///   both orientations for every valid k-mer.
pub fn scan_targets(
    sequence: &str,
    pam: &ParsedPAM,
    size: usize,
    right: bool,
    threads: usize,
) -> Result<(Vec<usize>, Vec<u8>), String> {
    // Convert sequence to IUPAC bitmasks (return error instead of panic).
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
        return Ok((Vec::new(), Vec::new()));
    }

    let pat = &pam.bytes;   // forward PAM bitmasks
    let rev = &pam.revcomp; // reverse-complement PAM bitmasks
    let plen = pat.len();

    if plen == 0 || plen > size {
        return Err(format!("Invalid PAM length: plen={plen}, size={size}"));
    }

    // Build sparse representations for PAM matching (skip unconstrained N positions).
    let (idx_fwd, mask_fwd) = build_sparse(pat);
    let (idx_rev, mask_rev) = build_sparse(rev);

    // Heuristic: use sparse if it reduces the number of checked positions.
    let use_sparse_fwd = !pam.unconstrained && idx_fwd.len() < plen;
    let use_sparse_rev = !pam.unconstrained && idx_rev.len() < plen;

    // PAM start indices within the scan window (invariant across k-mers).
    let pam_start_fwd = if right { 0 } else { size - plen };
    let pam_start_rev = if right { size - plen } else { 0 };

    // Run the parallel scan inside a cached pool for this `threads` value.
    threadpool::with_pool(threads, || {
        // Chunk size (ceiling division).
        let chunk_size = (slen + threads - 1) / threads;

        // Each worker returns its own `(positions, strands)` buffers.
        let per_chunk: Vec<(Vec<usize>, Vec<u8>)> = (0..threads)
            .into_par_iter()
            .filter_map(|chunk_idx| {
                let orig_start = chunk_idx * chunk_size;
                if orig_start >= slen {
                    return None;
                }
                let orig_end = std::cmp::min(orig_start + chunk_size, slen);

                // Extend by (size - 1) so k-mers starting within [orig_start, orig_end) are complete.
                let extended_start = orig_start;
                let extended_end = std::cmp::min(orig_end + (size - 1), slen);

                let extended_chunk = &seq_bitmask[extended_start..extended_end];
                let chunk_len = extended_chunk.len();

                let mut chunk_pos: Vec<usize> = Vec::new();
                let mut chunk_strand: Vec<u8> = Vec::new(); // 1 = fwd, 0 = rev

                if chunk_len >= size {
                    for i in 0..=(chunk_len - size) {
                        let global_pos = extended_start + i;

                        // Avoid duplicates across chunks: only emit k-mers whose start
                        // lies within the non-extended chunk interval.
                        if global_pos < orig_start || global_pos >= orig_end {
                            continue;
                        }

                        // SAFETY: i..i+size is in-bounds because i <= chunk_len - size.
                        let target_bitmask = unsafe { extended_chunk.get_unchecked(i..i + size) };

                        // Skip k-mers containing an ambiguous base 'N' in the sequence.
                        if target_bitmask.iter().any(|&b| b == N_MASK) {
                            continue;
                        }

                        if pam.unconstrained {
                            // Fully unconstrained PAM (NNN...): emit both orientations.
                            chunk_pos.push(global_pos);
                            chunk_strand.push(1);
                            chunk_pos.push(global_pos);
                            chunk_strand.push(0);
                        } else {
                            // Forward orientation PAM match (sparse or dense).
                            let fwd_ok = if use_sparse_fwd {
                                matches_pattern_sparse(
                                    target_bitmask,
                                    pam_start_fwd,
                                    &idx_fwd,
                                    &mask_fwd,
                                )
                            } else {
                                let pam_slice_fwd =
                                    &target_bitmask[pam_start_fwd..pam_start_fwd + plen];
                                matches_pattern(pam_slice_fwd, pat)
                            };

                            // Reverse orientation PAM match (sparse or dense).
                            let rev_ok = if use_sparse_rev {
                                matches_pattern_sparse(
                                    target_bitmask,
                                    pam_start_rev,
                                    &idx_rev,
                                    &mask_rev,
                                )
                            } else {
                                let pam_slice_rev =
                                    &target_bitmask[pam_start_rev..pam_start_rev + plen];
                                matches_pattern(pam_slice_rev, rev)
                            };

                            if fwd_ok {
                                chunk_pos.push(global_pos);
                                chunk_strand.push(1);
                            }
                            if rev_ok {
                                chunk_pos.push(global_pos);
                                chunk_strand.push(0);
                            }
                        }
                    }
                }

                Some((chunk_pos, chunk_strand))
            })
            .collect();

        // Merge into final vectors (single-threaded, linear).
        let total_hits: usize = per_chunk.iter().map(|(p, _)| p.len()).sum();
        let mut pos: Vec<usize> = Vec::with_capacity(total_hits);
        let mut strand: Vec<u8> = Vec::with_capacity(total_hits);

        for (p, s) in per_chunk {
            // Invariant: positions and strands are parallel arrays.
            debug_assert_eq!(p.len(), s.len());
            pos.extend(p);
            strand.extend(s);
        }

        (pos, strand)
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


/// Builds a sparse representation of a PAM pattern by retaining only *informative* positions.
///
/// In IUPAC encoding, the mask `0b1111` (`N`) matches any nucleotide and therefore
/// does not constrain matching. This function filters out such positions and returns:
///   - the indices of PAM positions that impose constraints, and
///   - the corresponding IUPAC bitmasks.
///
/// This representation reduces per-k-mer matching work for partially-degenerate PAMs
/// (e.g., `NNGRRT`, `GGNRG`) by checking only informative positions.
///
/// # Arguments
/// * `pam` - Slice of IUPAC bitmasks representing the PAM sequence.
///
/// # Returns
/// A tuple `(idx, mask)` where:
/// * `idx[i]` is the position within the PAM of the `i`-th informative base.
/// * `mask[i]` is the corresponding IUPAC bitmask at that position.
///
/// # Notes
/// * If all PAM positions are unconstrained (`N`), both vectors will be empty.
/// * If no positions are unconstrained, `idx.len() == pam.len()`.
#[inline]
fn build_sparse(pam: &[u8]) -> (Vec<usize>, Vec<u8>) {
    // define vectors of indeexes and masks
    let mut idx: Vec<usize> = Vec::new();
    let mut mask: Vec<u8> = Vec::new();

    // iterate over pam nts
    for (i, &m) in pam.iter().enumerate() {
        if m != N_MASK {
            idx.push(i);
            mask.push(m);
        }
    }

    (idx, mask)

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
