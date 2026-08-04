//! Strand-aware, single-linkage clustering over rows pre-sorted by
//! `(contig, strand, start)`, plus per-cluster primary selection.
//!
//! A cluster collapses the redundant representations of one on-strand alignment
//! (bulge/gap-shifted starts) into a primary and its alternatives. Clustering is
//! keyed on `(contig, strand, start)` — the strand leg guarantees a genuine
//! opposite-strand off-target at the same locus survives as its own primary
//! rather than being demoted by the tiebreak.

use std::ops::Range;

use super::{primary_cmp, AlignmentRow, PrimaryCriteria};

/// Do two *adjacent* rows (in `(contig, strand, start)` order) belong to the
/// same cluster? Strand-aware, single-linkage, inclusive of a gap of exactly `n`.
///
/// Contract: `cur` follows `prev` in a slice sorted by
/// `(contig, strand, start ascending)`, so `cur.start >= prev.start` within a
/// `(contig, strand)` group. `saturating_sub` keeps this defined if the contract
/// is ever violated (a violation can only *split* a run, never merge wrongly).
#[inline]
pub fn same_cluster(prev: &AlignmentRow, cur: &AlignmentRow, n: u32) -> bool {
    prev.contig == cur.contig
        && prev.strand == cur.strand
        && cur.start.saturating_sub(prev.start) <= n
}

/// Sweeps a sorted slice into maximal single-linkage cluster runs, yielding the
/// index [`Range`] of each run. Zero-allocation.
pub struct ClusterRuns<'a> {
    rows: &'a [AlignmentRow],
    n: u32,
    start: usize,
}

impl<'a> Iterator for ClusterRuns<'a> {
    type Item = Range<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start >= self.rows.len() {
            return None;
        }
        let begin = self.start;
        let mut end = begin + 1;
        while end < self.rows.len() && same_cluster(&self.rows[end - 1], &self.rows[end], self.n) {
            end += 1;
        }
        self.start = end;
        Some(begin..end)
    }
}

/// Iterate the cluster runs of a slice pre-sorted by `(contig, strand, start)`.
#[inline]
pub fn cluster_runs(rows: &[AlignmentRow], n: u32) -> ClusterRuns<'_> {
    ClusterRuns { rows, n, start: 0 }
}

/// Absolute index of the primary row within a (non-empty) cluster `run`, under
/// `criteria`. Every other row in the run is an alternative.
pub fn select_primary(
    rows: &[AlignmentRow],
    run: Range<usize>,
    criteria: &PrimaryCriteria,
) -> usize {
    run.min_by(|&i, &j| primary_cmp(&rows[i], &rows[j], criteria))
        .expect("cluster runs are non-empty by construction")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::occurence::Strand;
    use crate::partition::PrimaryCriteria;

    fn row(contig: &str, strand: Strand, start: u32, mm: u8, bdna: u8, brna: u8) -> AlignmentRow {
        AlignmentRow {
            contig: contig.into(),
            strand,
            start,
            mm,
            bdna,
            brna,
            scores: [f32::NAN; 4],
            target_aligned: "AAAA".into(),
        }
    }

    #[test]
    fn same_cluster_within_and_beyond_n() {
        let p = row("1", Strand::Forward, 100, 0, 0, 0);
        assert!(same_cluster(
            &p,
            &row("1", Strand::Forward, 103, 0, 0, 0),
            3
        ));
        assert!(!same_cluster(
            &p,
            &row("1", Strand::Forward, 104, 0, 0, 0),
            3
        ));
    }

    #[test]
    fn same_cluster_respects_strand_and_contig() {
        let p = row("1", Strand::Forward, 100, 0, 0, 0);
        assert!(!same_cluster(
            &p,
            &row("1", Strand::Reverse, 101, 0, 0, 0),
            3
        ));
        assert!(!same_cluster(
            &p,
            &row("2", Strand::Forward, 101, 0, 0, 0),
            3
        ));
    }

    #[test]
    fn runs_chain_by_single_linkage() {
        let rows = vec![
            row("1", Strand::Forward, 100, 0, 0, 0),
            row("1", Strand::Forward, 103, 0, 0, 0),
            row("1", Strand::Forward, 106, 0, 0, 0),
        ];
        assert_eq!(cluster_runs(&rows, 3).collect::<Vec<_>>(), vec![0..3]);

        let rows = vec![
            row("1", Strand::Forward, 100, 0, 0, 0),
            row("1", Strand::Forward, 104, 0, 0, 0),
        ];
        assert_eq!(cluster_runs(&rows, 3).collect::<Vec<_>>(), vec![0..1, 1..2]);
    }

    #[test]
    fn runs_split_on_strand_and_contig() {
        let rows = vec![
            row("1", Strand::Forward, 100, 0, 0, 0),
            row("1", Strand::Reverse, 101, 0, 0, 0),
        ];
        assert_eq!(cluster_runs(&rows, 3).collect::<Vec<_>>(), vec![0..1, 1..2]);

        let rows = vec![
            row("1", Strand::Forward, 100, 0, 0, 0),
            row("2", Strand::Forward, 101, 0, 0, 0),
        ];
        assert_eq!(cluster_runs(&rows, 3).collect::<Vec<_>>(), vec![0..1, 1..2]);
    }

    #[test]
    fn empty_and_singleton_sweeps() {
        assert!(cluster_runs(&[], 3).next().is_none());
        let rows = vec![row("1", Strand::Forward, 100, 0, 0, 0)];
        assert_eq!(cluster_runs(&rows, 3).collect::<Vec<_>>(), vec![0..1]);
    }

    #[test]
    fn select_primary_picks_min() {
        let c = PrimaryCriteria::default();
        let rows = vec![
            row("1", Strand::Forward, 100, 3, 0, 0),
            row("1", Strand::Forward, 101, 2, 0, 0), // primary (edit 2)
            row("1", Strand::Forward, 102, 4, 0, 0),
        ];
        assert_eq!(select_primary(&rows, 0..3, &c), 1);
    }

    #[test]
    fn select_primary_uses_tiebreak() {
        let c = PrimaryCriteria::default();
        let mut a = row("1", Strand::Forward, 100, 2, 0, 0);
        let mut b = row("1", Strand::Forward, 101, 2, 0, 0);
        a.target_aligned = "TTTT".into();
        b.target_aligned = "AAAA".into(); // smaller -> primary
        let rows = vec![a, b];
        assert_eq!(select_primary(&rows, 0..2, &c), 1);
    }
}
