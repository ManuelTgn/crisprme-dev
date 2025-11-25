use crate::pam::ParsedPAM;
use crate::iupac::matches_iupac;
use crate::target::Target;

use std::sync::Arc;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

pub fn scan_targets(
    sequence: &str, 
    contig: &str, 
    pam: &ParsedPAM, 
    k: usize, 
    right: bool, 
    threads: usize
) -> Vec<Target> {

    let seq = sequence.as_bytes();
    let pat = Arc::new(pam.bytes.clone());
    let rev = Arc::new(pam.revcomp.clone());
    let slen = seq.len();
    let plen = pat.len();

    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to create Rayon thread pool");

    let targets: Vec<Target> = pool.install(|| {
        // Calculate chunk size (ceiling division)
        let chunk_size = (slen + threads - 1) / threads;

        // Parallel iterate over chunk indices; produce Vec<Vec<(pos, String)>>
        (0..threads)
            .into_par_iter()
            .filter_map(|chunk_idx| {
                let orig_start = chunk_idx * chunk_size;
                if orig_start >= slen {
                    return None;
                }
                let orig_end = std::cmp::min(orig_start + chunk_size, slen);

                // Extend chunk by size-1 to capture k-mers crossing boundary
                let extended_start = orig_start;
                let extended_end = std::cmp::min(orig_end + (k - 1), slen);

                // Slice bytes (safe because DNA/RNA are ASCII)
                let chunk = &seq[extended_start..extended_end];
                let chunk_len = chunk.len();

                let mut chunk_targets = Vec::new();
                if chunk_len >= k {
                    for i in 0..=(chunk_len - k) {
                        let global_pos = extended_start + i;
                        // Emit only k-mers whose start is inside the original chunk (avoid duplicates)
                        if global_pos >= orig_start && global_pos < orig_end {
                            // retrieve target (bytes)
                            let target = &seq[i..i + k];
                            if target.contains(&b'N') {
                                continue;
                            }
                            // PAM match on forward strand
                            if matches_pattern(get_pam_slice(target, plen, k, right), &pat) {
                                chunk_targets.push(Target::new(contig, global_pos, true, std::str::from_utf8(target).unwrap()));
                            }
                            // PAM match on reverse strand
                            if matches_pattern(get_pam_slice(target, plen, k, !right), &rev) {
                                chunk_targets.push(Target::new(contig, global_pos, false, std::str::from_utf8(target).unwrap()));
                            }
                        }
                    }
                }
                Some(chunk_targets)
            })
            .flatten()
            .collect()
    });

    return targets;
}

fn matches_pattern(seq: &[u8], pam: &[u8]) -> bool {
    seq.iter()
        .zip(pam.iter())
        .all(|(a, b)| matches_iupac(*a, *b))
}


fn get_pam_slice(target: &[u8], plen: usize, k: usize, right: bool) -> &[u8] {
    if right {
        &target[0..plen]
    } else {
        &target[k - plen..k]
    }
}