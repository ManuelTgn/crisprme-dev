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
    let plen = pat.len();

    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to create Rayon thread pool");

    pool.install(|| {
        // parallel scan
        (0..=seq.len().saturating_sub(k))
            .into_par_iter()
            .filter_map(|i| {
                let target = &seq[i..i + k];

                //PAM at end or neginning of target depending on 'right'
                let pam_slice = if right {
                    &target[k - plen..k]
                } else {
                    &target[0..plen]
                };
                if matches_pattern(pam_slice, &pat) {
                    return Some(Target::new(contig, i, true, std::str::from_utf8(target).unwrap()));
                }

                //reverse complement
                if matches_pattern(pam_slice, &rev) {
                    return Some(Target::new(contig, i, false, std::str::from_utf8(target).unwrap()));
                }

                None

        })
        .collect()
    })
}

fn matches_pattern(seq: &[u8], pam: &[u8]) -> bool {
    seq.iter()
        .zip(pam.iter())
        .all(|(a, b)| matches_iupac(*a, *b))
}