//! Post-search finalization: split the single intermediate report into a
//! primary and an alternative report.
//!
//! This module is deliberately *outside* the streaming pipeline. Clustering is
//! a global, position-ordered operation and the parallel [`TsvWriterSink`] never
//! sees a whole cluster, so the work happens after the pipeline drains, over the
//! intermediate report on disk.
//!
//! Layout:
//!   - [`criteria`] — the user-configurable primary/alternative ranking.
//!   - [`cluster`]  — the strand-aware single-linkage sweep + primary selection.
//!
//! Step 3 (TSV parse -> sort -> route to two files, PyO3 + composition-root
//! wiring) builds on the pure functions re-exported here.
//!
//! [`TsvWriterSink`]: crate::pipeline::sink::writer::TsvWriterSink

pub mod cluster;
pub mod criteria;
pub mod partitioner;

use crate::error::crisprme_errors::PartitionError;
use crate::model::occurence::Strand;

pub use cluster::{cluster_runs, same_cluster, select_primary, ClusterRuns};
pub use criteria::{primary_cmp, PrimaryCriteria, PrimaryField, ScoreKind};
pub use partitioner::{partition_report, PartitionStats};

/// The decision-relevant fields of one intermediate-report row.
///
/// Parsed once from a report line; the verbatim line is carried separately by
/// the partitioner (step 3) so routing to primary/alternative preserves the
/// original formatting, scores, and annotation columns exactly.
///
/// `contig`, `strand`, and `start` together form the (non-negotiable) cluster
/// key; the score/bulge/mismatch fields feed the (configurable) comparator; and
/// `target_aligned` feeds the fixed tiebreak.
#[derive(Debug, Clone)]
pub struct AlignmentRow {
    /// Chromosome column, verbatim. Only equality + grouping matter here.
    pub contig: Box<str>,
    /// Reported strand (fixed part of the cluster key).
    pub strand: Strand,
    /// Reported left-most target position (fixed cluster key).
    pub start: u32,
    /// Mismatches (report `mismatches` column).
    pub mm: u8,
    /// DNA bulges (report `dna_bulges` column).
    pub bdna: u8,
    /// RNA bulges (report `rna_bulges` column).
    pub brna: u8,
    /// Specificity scores in report order (CFD, CRISTA, Elevation, aggregate).
    /// `NaN` encodes the report's `NA`.
    pub scores: [f32; 4],
    /// Aligned target column, used only by the fixed tiebreak.
    pub target_aligned: Box<str>,
}

impl AlignmentRow {
    /// Edit distance is fixed as `mm + bdna + brna` (widened to avoid overflow),
    /// the single source of truth for the [`PrimaryField::EditDistance`] term.
    #[inline(always)]
    pub fn edit_distance(&self) -> u16 {
        self.mm as u16 + self.bdna as u16 + self.brna as u16
    }
}
