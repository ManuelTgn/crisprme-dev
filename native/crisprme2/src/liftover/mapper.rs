//! Lift a single haplotype target position to the reference through a parsed 
//! [`ChainFile`], with a tolerance-gated ambiguity flag. 
//!
//! # Anchoring
//! Liftover is evaluated at **one** target coordinate per off-target — the
//! leftmost forward coordinate, the same anchor the hg38 clustering key uses —
//! so "does this site lift" is a single yes/no, not a fractional-overlap
//! judgement. Callers pass the site's leftmost forward target position.
//!
//! # Strand
//! A block from a reverse chain (`Chain::q_strand == Reverse`) flips the
//! reported reference strand relative to the haplotype strand. The mapper
//! returns the reference strand *relative to the haplotype hit's own strand*
//! (i.e. whether it was flipped).

use crate::liftover::chain::{Chain, ChainFile, Strand};

/// Outcome of lifting one target position, mirroring the report `MapStatus`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MapStatus {
    /// Uniquely mapped (or all covering chains agree within tolerance)
    Mapped,
    /// A lower-score covering chain disagrees by more than the tolerance
    Ambiguous,
    /// No chain covers this position (contig absent, or the base falls in a
    /// chain gap) — the site is assembly-specific and invisible to any
    /// reference-genome search
    Unmapped,
}

/// A successful lift: reference contig, 0-based position, strand relative to the
/// haplotype hit, and the quality status.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Lifted {
    pub ref_name: String,
    pub ref_pos: u64,
    /// Whether the reference strand is flipped vs. the haplotype hit's strand.
    pub flipped: bool,
    pub status: MapStatus,
}

/// Result of a lift attempt. `Unmapped`/`Novel` carry no coordinate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LiftResult {
    Ok(Lifted),
    NoMap(MapStatus),
}

/// Mapper over a parsed chain file.
pub struct LiftOver<'a> {
    chains: &'a ChainFile,
    tolerance: u64,
}

impl<'a> LiftOver<'a> {
    /// Build a mapper. `tolerance` (bp) is the ambiguity threshold
    pub fn new(chains: &'a ChainFile, tolerance: u64) -> Self {
        Self { chains, tolerance }
    }

    /// Lift one leftmost-forward target position on contig `t_name`
    pub fn lift(&self, t_name: &str, t_pos: u64) -> LiftResult {
        let mut primary: Option<(&Chain, u64, bool)> = None;
        let mut secondaries: Vec<u64> = Vec::new();
        for chain in self.chains.chains_for(t_name) {
            if let Some((q, flipped)) = map_in_chain(chain, t_pos) {
                match primary {
                    None => primary = Some((chain, q, flipped)),
                    Some(_) => secondaries.push(q),
                }
            }
        }
        // No covering chain (absent contig OR gap) -> assembly-specific
        let (primary_chain, q, flipped) = match primary {
            Some(p) => p,
            None => return LiftResult::NoMap(MapStatus::Unmapped),
        };
        let ambiguous = secondaries.iter().any(|&q2| abs_diff(q, q2) > self.tolerance);
        LiftResult::Ok(Lifted {
            ref_name: primary_chain.q_name.clone(),
            ref_pos: q,
            flipped,
            status: if ambiguous { MapStatus::Ambiguous } else { MapStatus::Mapped },
        })
    }
}

/// Map a target position within a single chain to its reference position.
///
/// Returns `(ref_pos, flipped)` if some block covers `t_pos`, else `None`
/// (the position lies in a chain gap). Binary-searches the block whose
/// `[t_start, t_end)` contains `t_pos`.
///
/// Within a block the target and reference run in lockstep (both plus-strand,
/// equal length): for a forward chain `q = q_start + (t_pos - t_start)`; for a
/// reverse chain the block's reference interval is reversed relative to target,
/// so `q = q_end - 1 - (t_pos - t_start)`.
fn map_in_chain(chain: &Chain, t_pos: u64) -> Option<(u64, bool)> {
    if t_pos < chain.t_start || t_pos >= chain.t_end {
        return None;
    }
    let blocks = &chain.blocks;
    // rightmost block with t_start <= t_pos
    let idx = match blocks.binary_search_by(|b| b.t_start.cmp(&t_pos)) {
        Ok(i) => i,
        Err(0) => return None, // before the first block
        Err(i) => i - 1,
    };
    let b = &blocks[idx];
    if t_pos >= b.t_end {
        return None; // in the gap between this block and the next
    }
    let off = t_pos - b.t_start;
    let (q, flipped) = match chain.q_strand {
        Strand::Forward => (b.q_start + off, false),
        Strand::Reverse => (b.q_end - 1 - off, true),
    };
    Some((q, flipped))
}

#[inline]
fn abs_diff(a: u64, b: u64) -> u64 {
    if a >= b { a - b } else { b - a }
}
