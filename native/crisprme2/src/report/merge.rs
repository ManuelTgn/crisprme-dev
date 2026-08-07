//! Within-sample merge (D2): collapse one sample's lifted haplotype reports
//! into canonical [`MergedRow`]s.
//!
//! Reads every lifted report for a sample, stamps each row with its copy-
//! presence bit (PanSN contig -> hap_id -> layout bit), single-linkage-clusters
//! on the sequence-aware key, and folds each cluster to one representative:
//! metrics/scores from the primary-criteria winner, annotation OR-ed, and the
//! per-sample presence OR-ed (the homozygous collapse). Interning happens once
//! per cluster, so only folded sets enter the registry.

use std::io::BufRead;

use crate::error::crisprme_errors::ReportError;
use crate::liftover::mapper::MapStatus;
use crate::model::occurence::Strand;
use crate::partition::criteria::{primary_cmp, PrimaryCriteria};
use crate::partition::AlignmentRow;
use crate::report::row::{ClusterKey, MergedRow, NativeLoc, RefLoc, Scores};
use crate::report::samples::{Presence, SamplePresence, SampleSetRegistry, SampleTable};

/// Column names D2 resolves from the lifted report header (order-independent).
mod col {
    pub const CHROM: &str = "chromosome";
    pub const START: &str = "start";
    pub const STRAND: &str = "strand";
    pub const GUIDE: &str = "sgRNA_aligned";
    pub const TARGET: &str = "target_aligned";
    pub const MM: &str = "mismatches";
    pub const BDNA: &str = "dna_bulges";
    pub const BRNA: &str = "rna_bulges";
    pub const CFD: &str = "CFD_score";
    pub const CRISTA: &str = "CRISTA_score";
    pub const ELEVATION: &str = "Elevation_score";
    pub const AGGREGATE: &str = "aggregate_score";
    pub const ANNOTATION: &str = "annotation";
    pub const REF_CHR: &str = "hg38_chr";
    pub const REF_POS: &str = "hg38_pos";
    pub const REF_FLIP: &str = "hg38_flipped";
    pub const MAP_STATUS: &str = "map_status";
}

/// A parsed row plus its single-copy presence, before clustering/interning.
/// Carries everything `MergedRow` needs except the interned `sample_set`, which
/// is minted once per cluster in [`fold_cluster`].
struct PreRow {
    guide_aligned: Box<str>,
    target_aligned: Box<str>,
    mm: u8,
    bdna: u8,
    brna: u8,
    native: NativeLoc,
    reference: Option<RefLoc>,
    scores: Scores,
    annotation: u32,
    /// This row's sample + single-copy bit (`1 << layout_bit`).
    presence: SamplePresence,
}

impl PreRow {
    /// Same sequence-aware key as [`MergedRow::cluster_key`] — reference coords
    /// when mapped, native otherwise, sub-keyed by allele.
    fn cluster_key(&self) -> ClusterKey {
        match &self.reference {
            Some(r) => ClusterKey::Mapped {
                contig: r.contig.clone(),
                strand: r.strand.as_bit(),
                start: r.start,
                allele: self.target_aligned.clone(),
            },
            None => ClusterKey::Unmapped {
                contig: self.native.contig.clone(),
                strand: self.native.strand.as_bit(),
                start: self.native.start,
                allele: self.target_aligned.clone(),
            },
        }
    }

    #[inline]
    fn cluster_pos(&self) -> u32 {
        match &self.reference {
            Some(r) => r.start,
            None => self.native.start,
        }
    }

    /// Borrow-free `AlignmentRow` view for ranking via the partitioner's
    /// `primary_cmp`, so the merge and the primary/alternative split use one
    /// definition of "best". Scores are laid out CFD, CRISTA, Elevation,
    /// aggregate — the `ScoreKind` order.
    fn as_alignment_row(&self) -> AlignmentRow {
        AlignmentRow {
            contig: self.native.contig.as_str().into(),
            strand: self.native.strand,
            start: self.native.start,
            mm: self.mm,
            bdna: self.bdna,
            brna: self.brna,
            scores: [
                self.scores.cfd,
                self.scores.crista,
                self.scores.elevation,
                self.scores.aggregate,
            ],
            target_aligned: self.target_aligned.clone(),
        }
    }
}

/// One haplotype report of a sample to merge.
pub struct HaplotypeReport<'a> {
    pub path: &'a str,
    /// Sample vocabulary index (from the manifest / `SampleTable`).
    pub sample: u32,
    /// Expected PanSN hap_id for this file; every row's contig must agree.
    pub expected_hap: u32,
}

/// Merge one sample's lifted haplotype reports into canonical rows.
///
/// `merge_bp` is the single-linkage window (the same value as the underlying
/// search's `--merge`). `criteria` selects the cluster representative via the
/// partitioner's `primary_cmp`, so "best" is defined once across the merge and
/// the primary/alternative split. `open` yields a reader per path (injected so
/// tests need no filesystem).
pub fn merge_sample(
    reports: &[HaplotypeReport],
    table: &SampleTable,
    registry: &mut SampleSetRegistry,
    merge_bp: u32,
    criteria: &PrimaryCriteria,
    open: &dyn Fn(&str) -> std::io::Result<Box<dyn BufRead>>,
) -> Result<Vec<MergedRow>, ReportError> {
    // 1. parse every row of every haplotype into PreRows (presence stamped)
    let mut rows: Vec<PreRow> = Vec::new();
    for rep in reports {
        let reader = open(rep.path)
            .map_err(|e| ReportError::io(format!("opening lifted report {}", rep.path), e))?;
        parse_report(reader, rep, table, &mut rows)?;
    }

    // 2. sort by the sequence-aware key: mapped block then unmapped; within
    //    each, (contig, strand, start, allele)
    rows.sort_by(|a, b| a.cluster_key().cmp(&b.cluster_key()));

    // 3. single-linkage sweep within (space, contig, strand, allele); fold each
    //    cluster to one MergedRow, interning the OR-ed presence once
    let mut out: Vec<MergedRow> = Vec::with_capacity(rows.len());
    let mut i = 0;
    while i < rows.len() {
        let mut j = i + 1;
        let mut edge = rows[i].cluster_pos();
        while j < rows.len() && same_cluster(&rows[i], &rows[j], edge, merge_bp) {
            edge = rows[j].cluster_pos(); // chained-gap: advance the edge
            j += 1;
        }
        out.push(fold_cluster(&rows[i..j], registry, criteria));
        i = j;
    }
    Ok(out)
}

/// Same cluster iff identical key except coordinate, and `b`'s position is
/// within `merge_bp` of the running edge. Mapped never clusters with unmapped.
fn same_cluster(a: &PreRow, b: &PreRow, edge: u32, merge_bp: u32) -> bool {
    match (a.cluster_key(), b.cluster_key()) {
        (
            ClusterKey::Mapped {
                contig: ca,
                strand: sa,
                allele: la,
                ..
            },
            ClusterKey::Mapped {
                contig: cb,
                strand: sb,
                allele: lb,
                ..
            },
        )
        | (
            ClusterKey::Unmapped {
                contig: ca,
                strand: sa,
                allele: la,
                ..
            },
            ClusterKey::Unmapped {
                contig: cb,
                strand: sb,
                allele: lb,
                ..
            },
        ) => ca == cb && sa == sb && la == lb && b.cluster_pos().saturating_sub(edge) <= merge_bp,
        _ => false,
    }
}

/// Fold a cluster into one `MergedRow`: representative by `primary_cmp`
/// (`Ordering::Less` == more primary), annotation OR-ed, per-sample presence
/// OR-ed and interned once.
fn fold_cluster(
    cluster: &[PreRow],
    registry: &mut SampleSetRegistry,
    criteria: &PrimaryCriteria,
) -> MergedRow {
    let rep = (0..cluster.len())
        .min_by(|&i, &j| {
            primary_cmp(
                &cluster[i].as_alignment_row(),
                &cluster[j].as_alignment_row(),
                criteria,
            )
        })
        .expect("cluster is non-empty by construction");
    let r = &cluster[rep];

    let annotation = cluster.iter().fold(0u32, |acc, x| acc | x.annotation);
    let presence: Vec<SamplePresence> = cluster.iter().map(|x| x.presence).collect();
    let sample_set = registry.intern(&presence); // per-sample OR fold on intern

    MergedRow {
        guide_aligned: r.guide_aligned.clone(),
        target_aligned: r.target_aligned.clone(),
        mm: r.mm,
        bdna: r.bdna,
        brna: r.brna,
        native: r.native.clone(),
        reference: r.reference.clone(),
        scores: r.scores,
        annotation,
        sample_set,
    }
}

/// Parse one haplotype report into `PreRow`s. The presence bit comes from the
/// row's own PanSN contig, validated against `rep.expected_hap`.
fn parse_report(
    mut reader: Box<dyn BufRead>,
    rep: &HaplotypeReport,
    table: &SampleTable,
    out: &mut Vec<PreRow>,
) -> Result<(), ReportError> {
    let mut header = String::new();
    if reader
        .read_line(&mut header)
        .map_err(|e| ReportError::io("reading header", e))?
        == 0
    {
        return Err(ReportError::MissingHeader);
    }
    let sc = Schema::resolve(header.trim_end())?;

    let bit = table
        .bit_of(rep.sample, rep.expected_hap)
        .ok_or_else(|| ReportError::BadRow {
            line: 1,
            msg: format!(
                "hap_id {} not in sample {}'s declared layout",
                rep.expected_hap, rep.sample
            ),
        })?;
    let self_mask: Presence = 1 << bit;

    let mut line = String::new();
    let mut ln = 1usize;
    loop {
        line.clear();
        ln += 1;
        if reader
            .read_line(&mut line)
            .map_err(|e| ReportError::io("reading report", e))?
            == 0
        {
            break;
        }
        let t = line.trim_end();
        if t.is_empty() {
            continue;
        }
        let f: Vec<&str> = t.split('\t').collect();
        let get = |i: usize| -> Result<&str, ReportError> {
            f.get(i).copied().ok_or_else(|| ReportError::BadRow {
                line: ln,
                msg: format!("missing column index {i}"),
            })
        };
        let pu32 = |v: &str, what: &str| -> Result<u32, ReportError> {
            v.parse::<u32>().map_err(|_| ReportError::BadRow {
                line: ln,
                msg: format!("invalid {what} {v:?}"),
            })
        };
        let pu8 = |v: &str, what: &str| -> Result<u8, ReportError> {
            v.parse::<u8>().map_err(|_| ReportError::BadRow {
                line: ln,
                msg: format!("invalid {what} {v:?}"),
            })
        };
        let pf = |v: &str| -> f32 {
            if v == "NA" {
                f32::NAN
            } else {
                v.parse::<f32>().unwrap_or(f32::NAN)
            }
        };

        let native_strand = parse_strand(get(sc.strand)?, ln)?;
        let contig = get(sc.chrom)?.to_string();
        validate_hap(&contig, rep, ln)?;

        let status = get(sc.map_status)?;
        let reference = match status {
            "mapped" | "ambiguous" => {
                let flipped = get(sc.ref_flip)? == "1";
                Some(RefLoc {
                    contig: get(sc.ref_chr)?.to_string(),
                    start: pu32(get(sc.ref_pos)?, "hg38_pos")?,
                    strand: xor_strand(native_strand, flipped),
                    status: if status == "ambiguous" {
                        MapStatus::Ambiguous
                    } else {
                        MapStatus::Mapped
                    },
                })
            }
            "unmapped" => None,
            other => {
                return Err(ReportError::BadRow {
                    line: ln,
                    msg: format!("unknown map_status {other:?}"),
                })
            }
        };

        let annotation = match sc.annotation {
            Some(i) => pu32(get(i)?, "annotation").unwrap_or(0),
            None => 0,
        };

        out.push(PreRow {
            guide_aligned: get(sc.guide)?.into(),
            target_aligned: get(sc.target)?.into(),
            mm: pu8(get(sc.mm)?, "mismatches")?,
            bdna: pu8(get(sc.bdna)?, "dna_bulges")?,
            brna: pu8(get(sc.brna)?, "rna_bulges")?,
            native: NativeLoc {
                contig,
                start: pu32(get(sc.start)?, "start")?,
                strand: native_strand,
            },
            reference,
            scores: Scores {
                cfd: pf(get(sc.cfd)?),
                crista: pf(get(sc.crista)?),
                elevation: pf(get(sc.elevation)?),
                aggregate: pf(get(sc.aggregate)?),
            },
            annotation,
            presence: SamplePresence {
                sample: rep.sample,
                mask: self_mask,
            },
        });
    }
    Ok(())
}

/// Resolved column indices for one report header.
struct Schema {
    chrom: usize,
    start: usize,
    strand: usize,
    guide: usize,
    target: usize,
    mm: usize,
    bdna: usize,
    brna: usize,
    cfd: usize,
    crista: usize,
    elevation: usize,
    aggregate: usize,
    annotation: Option<usize>,
    ref_chr: usize,
    ref_pos: usize,
    ref_flip: usize,
    map_status: usize,
}

impl Schema {
    fn resolve(header: &str) -> Result<Self, ReportError> {
        let cols: Vec<&str> = header.split('\t').collect();
        let find = |name: &str| {
            cols.iter()
                .position(|&c| c == name)
                .ok_or_else(|| ReportError::MissingColumn(name.to_string()))
        };
        Ok(Self {
            chrom: find(col::CHROM)?,
            start: find(col::START)?,
            strand: find(col::STRAND)?,
            guide: find(col::GUIDE)?,
            target: find(col::TARGET)?,
            mm: find(col::MM)?,
            bdna: find(col::BDNA)?,
            brna: find(col::BRNA)?,
            cfd: find(col::CFD)?,
            crista: find(col::CRISTA)?,
            elevation: find(col::ELEVATION)?,
            aggregate: find(col::AGGREGATE)?,
            annotation: cols.iter().position(|&c| c == col::ANNOTATION),
            ref_chr: find(col::REF_CHR)?,
            ref_pos: find(col::REF_POS)?,
            ref_flip: find(col::REF_FLIP)?,
            map_status: find(col::MAP_STATUS)?,
        })
    }
}

fn parse_strand(s: &str, ln: usize) -> Result<Strand, ReportError> {
    match s {
        "+" => Ok(Strand::Forward),
        "-" => Ok(Strand::Reverse),
        other => Err(ReportError::BadRow {
            line: ln,
            msg: format!("invalid strand {other:?}"),
        }),
    }
}

#[inline]
fn xor_strand(s: Strand, flipped: bool) -> Strand {
    match (s, flipped) {
        (Strand::Forward, false) | (Strand::Reverse, true) => Strand::Forward,
        _ => Strand::Reverse,
    }
}

fn validate_hap(contig: &str, rep: &HaplotypeReport, ln: usize) -> Result<(), ReportError> {
    let mut it = contig.splitn(3, '#');
    let (_s, h) = (it.next(), it.next());
    let hap: u32 = h
        .and_then(|x| x.parse().ok())
        .ok_or_else(|| ReportError::BadRow {
            line: ln,
            msg: format!("contig {contig:?} is not PanSN sample#hap#contig"),
        })?;
    if hap != rep.expected_hap {
        return Err(ReportError::BadRow {
            line: ln,
            msg: format!("row hap {hap} != file's declared hap {}", rep.expected_hap),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const HDR: &str = "chromosome\tstart\tstrand\tsgRNA_aligned\ttarget_aligned\tmismatches\tdna_bulges\trna_bulges\tCFD_score\tCRISTA_score\tElevation_score\taggregate_score\tannotation\thg38_chr\thg38_pos\thg38_flipped\tmap_status";

    fn table() -> SampleTable {
        // one sample "HG1", diploid, PanSN haps {1, 2} -> bits {0, 1}
        SampleTable::new(vec!["HG1".into()], vec![vec![1, 2]])
    }

    fn opener(
        map: std::collections::HashMap<String, String>,
    ) -> impl Fn(&str) -> std::io::Result<Box<dyn BufRead>> {
        move |p: &str| {
            Ok(Box::new(Cursor::new(map.get(p).cloned().unwrap_or_default())) as Box<dyn BufRead>)
        }
    }

    fn row(
        contig: &str,
        start: u32,
        strand: &str,
        target: &str,
        mm: u8,
        refc: &str,
        refp: &str,
        flip: &str,
        status: &str,
    ) -> String {
        format!("{contig}\t{start}\t{strand}\tGUIDE\t{target}\t{mm}\t0\t0\t0.5\tNA\tNA\t0.5\t0\t{refc}\t{refp}\t{flip}\t{status}")
    }

    #[test]
    fn homozygous_site_collapses_to_1_1() {
        // same reference locus + allele on hap1 and hap2 -> one row, HG1:1|1
        let hap1 = format!(
            "{HDR}\n{}",
            row("HG1#1#c", 100, "+", "ACGT", 0, "chr2", "500", "0", "mapped")
        );
        let hap2 = format!(
            "{HDR}\n{}",
            row("HG1#2#c", 100, "+", "ACGT", 0, "chr2", "500", "0", "mapped")
        );
        let map = [("h1".to_string(), hap1), ("h2".to_string(), hap2)]
            .into_iter()
            .collect();
        let open = opener(map);
        let reports = [
            HaplotypeReport {
                path: "h1",
                sample: 0,
                expected_hap: 1,
            },
            HaplotypeReport {
                path: "h2",
                sample: 0,
                expected_hap: 2,
            },
        ];
        let table = table();
        let mut reg = SampleSetRegistry::new();
        let out = merge_sample(
            &reports,
            &table,
            &mut reg,
            3,
            &PrimaryCriteria::default(),
            &open,
        )
        .unwrap();
        assert_eq!(out.len(), 1);
        let mut s = String::new();
        reg.render(out[0].sample_set, &table, &mut s);
        assert_eq!(s, "HG1:1|1");
    }

    #[test]
    fn heterozygous_site_stays_one_copy() {
        // present only on hap1 -> HG1:1|0
        let hap1 = format!(
            "{HDR}\n{}",
            row("HG1#1#c", 100, "+", "ACGT", 0, "chr2", "500", "0", "mapped")
        );
        let hap2 = HDR.to_string(); // hap2 has no rows
        let map = [("h1".to_string(), hap1), ("h2".to_string(), hap2)]
            .into_iter()
            .collect();
        let open = opener(map);
        let reports = [
            HaplotypeReport {
                path: "h1",
                sample: 0,
                expected_hap: 1,
            },
            HaplotypeReport {
                path: "h2",
                sample: 0,
                expected_hap: 2,
            },
        ];
        let table = table();
        let mut reg = SampleSetRegistry::new();
        let out = merge_sample(
            &reports,
            &table,
            &mut reg,
            3,
            &PrimaryCriteria::default(),
            &open,
        )
        .unwrap();
        assert_eq!(out.len(), 1);
        let mut s = String::new();
        reg.render(out[0].sample_set, &table, &mut s);
        assert_eq!(s, "HG1:1|0");
    }

    #[test]
    fn different_allele_same_locus_stays_split() {
        // two alleles at one locus in one hap-1 file -> two rows
        let hap1 = format!(
            "{HDR}\n{}\n{}",
            row("HG1#1#c", 100, "+", "ACGT", 0, "chr2", "500", "0", "mapped"),
            row("HG1#1#c", 100, "+", "ACGA", 1, "chr2", "500", "0", "mapped")
        );
        let map = [("h1".to_string(), hap1)].into_iter().collect();
        let open = opener(map);
        let reports = [HaplotypeReport {
            path: "h1",
            sample: 0,
            expected_hap: 1,
        }];
        let table = table();
        let mut reg = SampleSetRegistry::new();
        let out = merge_sample(
            &reports,
            &table,
            &mut reg,
            3,
            &PrimaryCriteria::default(),
            &open,
        )
        .unwrap();
        assert_eq!(out.len(), 2); // two alleles at one locus -> two rows
    }

    #[test]
    fn unmapped_never_clusters_with_mapped() {
        let hap1 = format!(
            "{HDR}\n{}\n{}",
            row("HG1#1#c", 100, "+", "ACGT", 0, "chr2", "500", "0", "mapped"),
            row("HG1#1#c", 100, "+", "ACGT", 0, "NA", "NA", "NA", "unmapped")
        );
        let map = [("h1".to_string(), hap1)].into_iter().collect();
        let open = opener(map);
        let reports = [HaplotypeReport {
            path: "h1",
            sample: 0,
            expected_hap: 1,
        }];
        let table = table();
        let mut reg = SampleSetRegistry::new();
        let out = merge_sample(
            &reports,
            &table,
            &mut reg,
            3,
            &PrimaryCriteria::default(),
            &open,
        )
        .unwrap();
        assert_eq!(out.len(), 2); // one mapped, one assembly-specific
        assert_eq!(out.iter().filter(|r| r.is_mapped()).count(), 1);
    }
}
