//! Cross-sample global merge: merge every sample's rows into one report.
//!
//! Runs [`merge_sample`] per sample, concatenates the resulting `MergedRow`s,
//! then sweeps them again on the same sequence-aware key across samples —
//! `union`-ing the interned `sample_set`s (per-sample copy-wise OR via the
//! registry) so an off-target carried by several individuals collapses to one
//! row whose `samples` column lists them all. Representative selection reuses
//! the partitioner's `primary_cmp`, so "best" is defined once everywhere.

use crate::error::crisprme_errors::ReportError;
use crate::partition::criteria::{primary_cmp, PrimaryCriteria};
use crate::partition::AlignmentRow;
use crate::report::merge::{merge_sample, HaplotypeReport};
use crate::report::row::{ClusterKey, MergedRow};
use crate::report::samples::{SamplePresence, SampleSetRegistry, SampleTable};
use std::io::BufRead;

/// All haplotype reports for one sample.
pub struct SampleReports<'a> {
    pub sample: u32,
    pub haplotypes: Vec<HaplotypeReport<'a>>,
}

/// `AlignmentRow` view of a `MergedRow` for ranking (reference coords when
/// mapped, native otherwise — the space the row clusters in).
fn view(r: &MergedRow) -> AlignmentRow {
    let (contig, strand, start) = match &r.reference {
        Some(rf) => (rf.contig.as_str(), rf.strand, rf.start),
        None => (r.native.contig.as_str(), r.native.strand, r.native.start),
    };
    AlignmentRow {
        contig: contig.into(),
        strand,
        start,
        mm: r.mm,
        bdna: r.bdna,
        brna: r.brna,
        scores: [r.scores.cfd, r.scores.crista, r.scores.elevation, r.scores.aggregate],
        target_aligned: r.target_aligned.clone(),
    }
}

#[inline]
fn cluster_pos(r: &MergedRow) -> u32 {
    match &r.reference {
        Some(rf) => rf.start,
        None => r.native.start,
    }
}

/// Same cross-sample cluster iff identical key except coordinate, within
/// `merge_bp` of the running edge. Mapped never clusters with unmapped.
fn same_cluster(a: &MergedRow, b: &MergedRow, edge: u32, merge_bp: u32) -> bool {
    match (a.cluster_key(), b.cluster_key()) {
        (
            ClusterKey::Mapped { contig: ca, strand: sa, allele: la, .. },
            ClusterKey::Mapped { contig: cb, strand: sb, allele: lb, .. },
        )
        | (
            ClusterKey::Unmapped { contig: ca, strand: sa, allele: la, .. },
            ClusterKey::Unmapped { contig: cb, strand: sb, allele: lb, .. },
        ) => ca == cb && sa == sb && la == lb && b_cluster_pos_within(b, edge, merge_bp),
        _ => false,
    }
}
#[inline]
fn b_cluster_pos_within(b: &MergedRow, edge: u32, merge_bp: u32) -> bool {
    cluster_pos(b).saturating_sub(edge) <= merge_bp
}

/// Merge all samples into the final `MergedRow`s for one guide.
///
/// `merge_bp` is the single-linkage window (the search's `--merge`); `criteria`
/// selects representatives, threaded through from the CLI so it matches the
/// partitioner. `open` yields a reader per path (injected for testability).
pub fn merge_global(
    samples: &[SampleReports],
    table: &SampleTable,
    registry: &mut SampleSetRegistry,
    merge_bp: u32,
    criteria: &PrimaryCriteria,
    open: &dyn Fn(&str) -> std::io::Result<Box<dyn BufRead>>,
) -> Result<Vec<MergedRow>, ReportError> {
    // 1. within-sample merge (D2), already canonical per sample
    let mut rows: Vec<MergedRow> = Vec::new();
    for s in samples {
        let merged = merge_sample(&s.haplotypes, table, registry, merge_bp, criteria, open)?;
        rows.extend(merged);
    }

    // 2. global sort on the same sequence-aware key (mapped block, then
    //    unmapped; within each: contig, strand, start, allele)
    rows.sort_by(|a, b| a.cluster_key().cmp(&b.cluster_key()));

    // 3. cross-sample sweep: fold each cluster to one representative, union-ing
    //    the sample sets (per-sample copy-wise OR in the registry)
    let mut out: Vec<MergedRow> = Vec::with_capacity(rows.len());
    let mut i = 0;
    while i < rows.len() {
        let mut j = i + 1;
        let mut edge = cluster_pos(&rows[i]);
        while j < rows.len() && same_cluster(&rows[i], &rows[j], edge, merge_bp) {
            edge = cluster_pos(&rows[j]);
            j += 1;
        }
        out.push(fold_across_samples(&rows[i..j], registry, criteria));
        i = j;
    }
    Ok(out)
}

/// Fold a cross-sample cluster: representative by `primary_cmp`, annotation
/// OR-ed, and the `sample_set`s unioned (per-sample copy-wise OR via the
/// registry — the individual-carrier collapse).
fn fold_across_samples(
    cluster: &[MergedRow],
    registry: &mut SampleSetRegistry,
    criteria: &PrimaryCriteria,
) -> MergedRow {
    let rep = (0..cluster.len())
        .min_by(|&i, &j| primary_cmp(&view(&cluster[i]), &view(&cluster[j]), criteria))
        .expect("cluster is non-empty");
    let mut out = cluster[rep].clone();

    out.annotation = cluster.iter().fold(0u32, |acc, r| acc | r.annotation);

    // union all members' interned sets: gather their entries, re-intern once
    // (intern folds per-sample copies by OR, so HG1:1|0 + HG1:0|1 -> HG1:1|1)
    let mut pairs: Vec<SamplePresence> = Vec::new();
    for r in cluster {
        for &(sample, mask) in registry.entries_of(r.sample_set) {
            pairs.push(SamplePresence { sample, mask });
        }
    }
    out.sample_set = registry.intern(&pairs);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::merge::HaplotypeReport;
    use std::collections::HashMap;
    use std::io::Cursor;

    const HDR: &str = "chromosome\tstart\tstrand\tsgRNA_aligned\ttarget_aligned\tmismatches\tdna_bulges\trna_bulges\tCFD_score\tCRISTA_score\tElevation_score\taggregate_score\tannotation\thg38_chr\thg38_pos\thg38_flipped\tmap_status";

    fn r(contig: &str, target: &str) -> String {
        format!("{contig}\t100\t+\tG\t{target}\t0\t0\t0\t0.5\tNA\tNA\t0.5\t0\tchr2\t500\t0\tmapped")
    }
    fn table() -> SampleTable {
        // HG1, HG2 both diploid, PanSN haps {1,2}
        SampleTable::new(vec!["HG1".into(), "HG2".into()], vec![vec![1, 2], vec![1, 2]])
    }
    fn open(map: HashMap<String, String>) -> impl Fn(&str) -> std::io::Result<Box<dyn BufRead>> {
        move |p: &str| Ok(Box::new(Cursor::new(map.get(p).cloned().unwrap_or_default())) as Box<dyn BufRead>)
    }

    #[test]
    fn shared_site_lists_both_individuals() {
        // HG1 het (hap1 only), HG2 hom (both haps): same hg38 locus+allele
        let m: HashMap<String, String> = [
            ("g1h1".into(), format!("{HDR}\n{}", r("HG1#1#c", "ACGT"))),
            ("g1h2".into(), HDR.to_string()),
            ("g2h1".into(), format!("{HDR}\n{}", r("HG2#1#c", "ACGT"))),
            ("g2h2".into(), format!("{HDR}\n{}", r("HG2#2#c", "ACGT"))),
        ].into_iter().collect();
        let op = open(m);
        let samples = [
            SampleReports { sample: 0, haplotypes: vec![
                HaplotypeReport { path: "g1h1", sample: 0, expected_hap: 1 },
                HaplotypeReport { path: "g1h2", sample: 0, expected_hap: 2 },
            ]},
            SampleReports { sample: 1, haplotypes: vec![
                HaplotypeReport { path: "g2h1", sample: 1, expected_hap: 1 },
                HaplotypeReport { path: "g2h2", sample: 1, expected_hap: 2 },
            ]},
        ];
        let table = table();
        let mut reg = SampleSetRegistry::new();
        let out = merge_global(&samples, &table, &mut reg, 3, &PrimaryCriteria::default(), &op).unwrap();
        assert_eq!(out.len(), 1);
        let mut s = String::new();
        reg.render(out[0].sample_set, &table, &mut s);
        assert_eq!(s, "HG1:1|0,HG2:1|1"); // het HG1, hom HG2, one shared row
    }
}