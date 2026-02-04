use crate::pam::{ParsedPAM, build_sparse};
use crate::iupac::{matches_iupac, sequence_encoder};
use crate::threadpool;

use rayon::prelude::*;
use std::result::Result; 


/// IUPAC bitmask for `N` (any base).
///
/// In this representation, `N` is `0b1111` (A|C|G|T). It is used both as a
/// degenerate PAM code and as an ambiguity marker in the input sequence.
/// In the scanner, k-mers containing `N` are skipped (see `scan_targets`).
const N_MASK: u8 = 0b1111;

pub fn scan_targets(
    sequence: &str,
    pam: &ParsedPAM,
    size: usize,
    right: bool,
    threads: usize,
) -> Result<(Vec<usize>, Vec<u8>), String> {
    // encode contig chunk sequence in bits
    let seq_bitmask = sequence_encoder(sequence);

    // extract target candidates from contig chunk
    scan_targets_bitmask(&seq_bitmask, pam, size, right, threads)
}

pub fn scan_targets_bitmask(
    seq_bitmask: &[u8],
    pam: &ParsedPAM,
    size: usize,
    right: bool,
    threads: usize,
) -> Result<(Vec<usize>, Vec<u8>), String> {
    // get sequence length
    let slen = seq_bitmask.len();
    if slen == 0 || size == 0 || size > slen {
        return  Ok((Vec::new(), Vec::new()));
    }

    // get encoded pam 
    let pat = &pam.bytes;
    let rev = &pam.revcomp;
    let plen = pat.len();
    if plen == 0 || plen > size {
        return Err(format!("Invalid PAM length: plen={plen}, size={size}"));
    }

    // decide whether perform sparse pam matching (sparse pam example NNNRGG)
    let (idx_fwd, mask_fwd) = build_sparse(pat);
    let (idx_rev, mask_rev) = build_sparse(rev);

    let use_sparse_fwd = !pam.unconstrained && idx_fwd.len() < plen;
    let use_sparse_rev = !pam.unconstrained && idx_rev.len() < plen;

    // compute start positions for pam
    let pam_start_fwd = if right { 0 } else { size - plen };
    let pam_start_rev = if right { size - plen } else { 0 };

    threadpool::with_pool(threads, || {
        // define chunk sizes for threads spawning over sequences
        let chunk_size = (slen + threads - 1) / threads;

        let per_chunk: Vec<(Vec<usize>, Vec<u8>)> = (0..threads)
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
                let mut chunk_strand: Vec<u8> = Vec::new(); // 1 = +; 0 = -

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
                            continue;  // skip if candidate contains any N
                        }

                        if pam.unconstrained { // degenerate PAM -> skip PAM matching
                            chunk_pos.push(global_pos);
                            chunk_strand.push(1);
                            chunk_pos.push(global_pos);
                            chunk_strand.push(0);
                        } else {  // perform PAM matching to filter out candidates
                            let fwd_ok = if use_sparse_fwd {
                                matches_pattern_sparse(target_bitmask, pam_start_fwd, &idx_fwd, &mask_fwd)
                            } else {
                                let pam_slice_fwd = &target_bitmask[pam_start_fwd..pam_start_fwd + plen];
                                matches_pattern(pam_slice_fwd, pat)
                            };

                            let rev_ok = if use_sparse_rev {
                                matches_pattern_sparse(target_bitmask, pam_start_rev, &idx_rev, &mask_rev)
                            } else {
                                let pam_slice_rev = &target_bitmask[pam_start_rev..pam_start_rev + plen];
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

        // collect results from targets extraction
        let total_hits: usize = per_chunk.iter().map(|(p, _)| p.len()).sum();
        let mut pos: Vec<usize> = Vec::with_capacity(total_hits);
        let mut strand: Vec<u8> = Vec::with_capacity(total_hits);

        for (p, s) in per_chunk {
            debug_assert_eq!(p.len(), s.len());
            pos.extend(p);
            strand.extend(s);
        }

        (pos, strand)

    })
    
}


// pub fn scan_targets(
//     sequence: &str,
//     pam: &ParsedPAM,
//     size: usize,
//     right: bool,
//     threads: usize,
// ) -> Result<(Vec<usize>, Vec<u8>), String> {
//     // Convert sequence to IUPAC bitmasks (return error instead of panic).
//     let seq_bitmask: Vec<u8> = sequence
//         .as_bytes()
//         .iter()
//         .enumerate()
//         .map(|(i, &b)| {
//             Iupac::from_ascii(b)
//                 .map(|iupac| iupac.0)
//                 .map_err(|e| format!("Invalid base at position {i}: {e}"))
//         })
//         .collect::<Result<_, _>>()?;

//     let slen = seq_bitmask.len();
//     if slen == 0 || size == 0 || size > slen {
//         return Ok((Vec::new(), Vec::new()));
//     }

//     let pat = &pam.bytes;   // forward PAM bitmasks
//     let rev = &pam.revcomp; // reverse-complement PAM bitmasks
//     let plen = pat.len();

//     if plen == 0 || plen > size {
//         return Err(format!("Invalid PAM length: plen={plen}, size={size}"));
//     }

//     // Build sparse representations for PAM matching (skip unconstrained N positions).
//     let (idx_fwd, mask_fwd) = build_sparse(pat);
//     let (idx_rev, mask_rev) = build_sparse(rev);

//     // Heuristic: use sparse if it reduces the number of checked positions.
//     let use_sparse_fwd = !pam.unconstrained && idx_fwd.len() < plen;
//     let use_sparse_rev = !pam.unconstrained && idx_rev.len() < plen;

//     // PAM start indices within the scan window (invariant across k-mers).
//     let pam_start_fwd = if right { 0 } else { size - plen };
//     let pam_start_rev = if right { size - plen } else { 0 };

//     // Run the parallel scan inside a cached pool for this `threads` value.
//     threadpool::with_pool(threads, || {
//         // Chunk size (ceiling division).
//         let chunk_size = (slen + threads - 1) / threads;

//         // Each worker returns its own `(positions, strands)` buffers.
//         let per_chunk: Vec<(Vec<usize>, Vec<u8>)> = (0..threads)
//             .into_par_iter()
//             .filter_map(|chunk_idx| {
//                 let orig_start = chunk_idx * chunk_size;
//                 if orig_start >= slen {
//                     return None;
//                 }
//                 let orig_end = std::cmp::min(orig_start + chunk_size, slen);

//                 // Extend by (size - 1) so k-mers starting within [orig_start, orig_end) are complete.
//                 let extended_start = orig_start;
//                 let extended_end = std::cmp::min(orig_end + (size - 1), slen);

//                 let extended_chunk = &seq_bitmask[extended_start..extended_end];
//                 let chunk_len = extended_chunk.len();

//                 let mut chunk_pos: Vec<usize> = Vec::new();
//                 let mut chunk_strand: Vec<u8> = Vec::new(); // 1 = fwd, 0 = rev

//                 if chunk_len >= size {
//                     for i in 0..=(chunk_len - size) {
//                         let global_pos = extended_start + i;

//                         // Avoid duplicates across chunks: only emit k-mers whose start
//                         // lies within the non-extended chunk interval.
//                         if global_pos < orig_start || global_pos >= orig_end {
//                             continue;
//                         }

//                         // SAFETY: i..i+size is in-bounds because i <= chunk_len - size.
//                         let target_bitmask = unsafe { extended_chunk.get_unchecked(i..i + size) };

//                         // Skip k-mers containing an ambiguous base 'N' in the sequence.
//                         if target_bitmask.iter().any(|&b| b == N_MASK) {
//                             continue;
//                         }

//                         if pam.unconstrained {
//                             // Fully unconstrained PAM (NNN...): emit both orientations.
//                             chunk_pos.push(global_pos);
//                             chunk_strand.push(1);
//                             chunk_pos.push(global_pos);
//                             chunk_strand.push(0);
//                         } else {
//                             // Forward orientation PAM match (sparse or dense).
//                             let fwd_ok = if use_sparse_fwd {
//                                 matches_pattern_sparse(
//                                     target_bitmask,
//                                     pam_start_fwd,
//                                     &idx_fwd,
//                                     &mask_fwd,
//                                 )
//                             } else {
//                                 let pam_slice_fwd =
//                                     &target_bitmask[pam_start_fwd..pam_start_fwd + plen];
//                                 matches_pattern(pam_slice_fwd, pat)
//                             };

//                             // Reverse orientation PAM match (sparse or dense).
//                             let rev_ok = if use_sparse_rev {
//                                 matches_pattern_sparse(
//                                     target_bitmask,
//                                     pam_start_rev,
//                                     &idx_rev,
//                                     &mask_rev,
//                                 )
//                             } else {
//                                 let pam_slice_rev =
//                                     &target_bitmask[pam_start_rev..pam_start_rev + plen];
//                                 matches_pattern(pam_slice_rev, rev)
//                             };

//                             if fwd_ok {
//                                 chunk_pos.push(global_pos);
//                                 chunk_strand.push(1);
//                             }
//                             if rev_ok {
//                                 chunk_pos.push(global_pos);
//                                 chunk_strand.push(0);
//                             }
//                         }
//                     }
//                 }

//                 Some((chunk_pos, chunk_strand))
//             })
//             .collect();

//         // Merge into final vectors (single-threaded, linear).
//         let total_hits: usize = per_chunk.iter().map(|(p, _)| p.len()).sum();
//         let mut pos: Vec<usize> = Vec::with_capacity(total_hits);
//         let mut strand: Vec<u8> = Vec::with_capacity(total_hits);

//         for (p, s) in per_chunk {
//             // Invariant: positions and strands are parallel arrays.
//             debug_assert_eq!(p.len(), s.len());
//             pos.extend(p);
//             strand.extend(s);
//         }

//         (pos, strand)
//     })
// }


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


// pub fn scan_targets_streaming(
//     sequence: &str,
//     pam: &crate::pam::ParsedPAM,
//     size: usize,
//     right: bool,
//     threads: usize,
//     tx: Sender<HitBlock>,
// ) -> Result<(), String> {
//     // Convert to IUPAC bitmasks once.
//     let seq_bitmask: Vec<u8> = sequence_encoder(sequence);

//     let slen = seq_bitmask.len();
//     if slen == 0 || size == 0 || size > slen {
//         return Ok(());
//     }

//     let pat = &pam.bytes;
//     let rev = &pam.revcomp;
//     let plen = pat.len();
//     if plen == 0 || plen > size {
//         return Result::Err(format!("Invalid PAM length: plen={plen}, size={size}"));
//     }

//     let (idx_fwd, mask_fwd) = build_sparse(pat);
//     let (idx_rev, mask_rev) = build_sparse(rev);

//     let use_sparse_fwd = !pam.unconstrained && idx_fwd.len() < plen;
//     let use_sparse_rev = !pam.unconstrained && idx_rev.len() < plen;

//     let pam_start_fwd = if right { 0 } else { size - plen };
//     let pam_start_rev = if right { size - plen } else { 0 };

//     // Avoid division by zero / odd behavior.
//     let threads = threads.max(1);

//     threadpool::with_pool(threads, || {
//         let chunk_size = (slen + threads - 1) / threads;

//         (0..threads).into_par_iter().for_each(|chunk_idx| {
//             let orig_start = chunk_idx * chunk_size;
//             if orig_start >= slen {
//                 return;
//             }
//             let orig_end = std::cmp::min(orig_start + chunk_size, slen);

//             let extended_start = orig_start;
//             let extended_end = std::cmp::min(orig_end + (size - 1), slen);
//             let extended_chunk = &seq_bitmask[extended_start..extended_end];
//             let chunk_len = extended_chunk.len();

//             let mut chunk_pos: Vec<u32> = Vec::new();
//             let mut chunk_strand: Vec<u8> = Vec::new();

//             if chunk_len >= size {
//                 for i in 0..=(chunk_len - size) {
//                     let global_pos = extended_start + i;

//                     // avoid duplicates from overlap
//                     if global_pos < orig_start || global_pos >= orig_end {
//                         continue;
//                     }

//                     let target_bitmask = unsafe { extended_chunk.get_unchecked(i..i + size) };

//                     if target_bitmask.iter().any(|&b| b == N_MASK) {
//                         continue;
//                     }

//                     if pam.unconstrained {
//                         chunk_pos.push(global_pos as u32);
//                         chunk_strand.push(1);
//                         chunk_pos.push(global_pos as u32);
//                         chunk_strand.push(0);
//                     } else {
//                         let fwd_ok = if use_sparse_fwd {
//                             matches_pattern_sparse(
//                                 target_bitmask,
//                                 pam_start_fwd,
//                                 &idx_fwd,
//                                 &mask_fwd,
//                             )
//                         } else {
//                             let pam_slice_fwd =
//                                 &target_bitmask[pam_start_fwd..pam_start_fwd + plen];
//                             matches_pattern(pam_slice_fwd, pat)
//                         };

//                         let rev_ok = if use_sparse_rev {
//                             matches_pattern_sparse(
//                                 target_bitmask,
//                                 pam_start_rev,
//                                 &idx_rev,
//                                 &mask_rev,
//                             )
//                         } else {
//                             let pam_slice_rev =
//                                 &target_bitmask[pam_start_rev..pam_start_rev + plen];
//                             matches_pattern(pam_slice_rev, rev)
//                         };

//                         if fwd_ok {
//                             chunk_pos.push(global_pos as u32);
//                             chunk_strand.push(1);
//                         }
//                         if rev_ok {
//                             chunk_pos.push(global_pos as u32);
//                             chunk_strand.push(0);
//                         }
//                     }
//                 }
//             }

//             if !chunk_pos.is_empty() {
//                 // Backpressure is fine; ignore send error if receiver dropped.
//                 let _ = tx.send(HitBlock {
//                     pos: chunk_pos,
//                     strand: chunk_strand,
//                 });
//             }
//         });

//         Ok(())
//     });

//     Ok(())
// }
