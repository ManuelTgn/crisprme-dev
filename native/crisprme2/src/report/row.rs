//! Canonical merged-row schema.
//!
//! One `MergedRow` is a cluster representative in the assembly report. It owns
//! its aligned-target sequence bytes (the allele half of the sequence-aware
//! clustering key) and carries three coordinate/quality blocks:
//!
//! * **native** — authoritative assembly coordinates, never overwritten by
//!   liftover;
//! * **reference** — hg38 coordinates + `MapStatus`, present only when the row
//!   lifted (`status == Mapped | Ambiguous`);
//! * **scores / annotation** — carried through from the scan verbatim.
//!
//! plus the interned `SampleSetId` (see [`super::samples`]).
//!
//! # Clustering key
//! Global merge sorts and single-linkage-clusters mapped rows by
//! `(ref_contig, ref_strand, ref_start)` within `merge_bp`, sub-keyed by
//! `target` (the allele) so same-locus / different-protospacer hits stay
//! distinct. Unmapped rows carry no reference key and cluster only among
//! themselves by native coordinate — encoded here by `reference: Option<_>`.

use super::samples::SampleSetId;
use crate::liftover::mapper::MapStatus;
use crate::model::occurence::Strand;

/// Native (assembly) coordinates of a row — authoritative, never rewritten.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeLoc {
    /// Assembly contig name as written in the intermediate report (the PanSN
    /// `sample#hap#contig` string).
    pub contig: String,
    /// Leftmost forward coordinate (the anchor used for liftover).
    pub start: u32,
    pub strand: Strand,
}

/// Reference (hg38) coordinates + liftover quality. Absent when the row did not
/// lift (`MapStatus::Unmapped`); present for `Mapped` / `Ambiguous`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefLoc {
    pub contig: String,
    pub start: u32,
    /// Absolute hg38 strand = native strand XOR liftover flip. Harmonised here
    /// (the `flipped` flag is resolved away, per the C2/C3 agreement), so
    /// clustering and the final report read a single, absolute strand.
    pub strand: Strand,
    pub status: MapStatus, // Mapped | Ambiguous (Unmapped => reference is None)
}

/// Per-row score payload, carried through verbatim from the scan. `f32::NAN`
/// renders as `NA` and is always least-primary, as in the reference path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Scores {
    pub cfd: f32,
    pub crista: f32,
    pub elevation: f32,
    pub aggregate: f32,
}

impl Default for Scores {
    fn default() -> Self {
        Self { cfd: f32::NAN, crista: f32::NAN, elevation: f32::NAN, aggregate: f32::NAN }
    }
}

/// A cluster representative in the assembly report.
///
/// Owns its sequence bytes; therefore not `Copy`. The global sort moves rows
/// (or sorts an index of them); either is fine at report scale.
#[derive(Debug, Clone, PartialEq)]
pub struct MergedRow {
    // --- identity / allele (sequence-aware key) ---
    /// Guide as aligned, with bulge gaps (`-`) — carried for the report.
    pub guide_aligned: Box<str>,
    /// Aligned target sequence (protospacer+PAM as aligned, bulge gaps `-`).
    /// The allele half of the clustering key; owned by the row.
    pub target_aligned: Box<str>,

    // --- quality metrics (report columns + primary-selection inputs) ---
    pub mm: u8,
    pub bdna: u8,
    pub brna: u8,

    // --- coordinates ---
    pub native: NativeLoc,
    /// `Some` iff the row lifted; `None` is an assembly-specific site.
    pub reference: Option<RefLoc>,

    // --- payload carried through the merge ---
    pub scores: Scores,
    /// Functional-annotation bitmask (hg38-space; `0` when unannotated/unmapped).
    pub annotation: u32,
    /// Interned per-sample haplotype-presence set (see [`super::samples`]).
    pub sample_set: SampleSetId,
}

impl MergedRow {
    /// Edit distance — the single source of truth, `mm + bdna + brna`.
    #[inline]
    pub fn edit_distance(&self) -> u16 {
        self.mm as u16 + self.bdna as u16 + self.brna as u16
    }

    /// Categorical bulge tag, matching the reference sink (`X` / `DNA` / `RNA`
    /// / `DNA/RNA`).
    #[inline]
    pub fn bulge_type(&self) -> &'static str {
        match (self.bdna > 0, self.brna > 0) {
            (true, true) => "DNA/RNA",
            (true, false) => "DNA",
            (false, true) => "RNA",
            (false, false) => "X",
        }
    }

    /// `true` if the row lifted to the reference (Mapped or Ambiguous).
    #[inline]
    pub fn is_mapped(&self) -> bool {
        self.reference.is_some()
    }
}

/// The clustering/sort key. Mapped rows key on reference coordinates + allele;
/// unmapped rows have no reference key and cluster only among themselves by
/// native coordinate. Ordering `Unmapped > Mapped` on the outer enum keeps
/// unmapped rows in their own contiguous block after sorting.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ClusterKey {
    /// `(ref_contig, ref_strand_bit, ref_start, target_allele)`.
    Mapped {
        contig: String,
        strand: u8,
        start: u32,
        allele: Box<str>,
    },
    /// `(native_contig, native_strand_bit, native_start, target_allele)`.
    Unmapped {
        contig: String,
        strand: u8,
        start: u32,
        allele: Box<str>,
    },
}

impl MergedRow {
    /// Build the sort/cluster key for this row. The `merge_bp` single-linkage
    /// sweep (D3) runs *within* a `(contig, strand)` group over ascending
    /// `start`; the allele sub-key is compared only at equal coordinates so
    /// distinct protospacer alleles at one locus stay separate rows.
    pub fn cluster_key(&self) -> ClusterKey {
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::samples::SampleSetId;

    fn row(target: &str, mm: u8, bdna: u8, brna: u8, reference: Option<RefLoc>) -> MergedRow {
        MergedRow {
            guide_aligned: "GUIDE".into(),
            target_aligned: target.into(),
            mm, bdna, brna,
            native: NativeLoc { contig: "S#1#CM1".into(), start: 100, strand: Strand::Forward },
            reference,
            scores: Scores::default(),
            annotation: 0,
            sample_set: SampleSetId(0),
        }
    }

    #[test]
    fn edit_distance_and_bulge_type() {
        let r = row("ACGT", 2, 1, 0, None);
        assert_eq!(r.edit_distance(), 3);
        assert_eq!(r.bulge_type(), "DNA");
        assert!(!r.is_mapped());
    }

    #[test]
    fn mapped_key_uses_reference_coords() {
        let refloc = RefLoc { contig: "chr2".into(), start: 500, strand: Strand::Reverse, status: MapStatus::Mapped };
        let r = row("ACGT", 0, 0, 0, Some(refloc));
        match r.cluster_key() {
            ClusterKey::Mapped { contig, start, allele, .. } => {
                assert_eq!((contig.as_str(), start, allele.as_ref()), ("chr2", 500, "ACGT"));
            }
            _ => panic!("mapped row must produce a Mapped key"),
        }
    }

    #[test]
    fn unmapped_sorts_after_mapped() {
        let mapped = ClusterKey::Mapped { contig: "chr2".into(), strand: 1, start: 0, allele: "A".into() };
        let unmapped = ClusterKey::Unmapped { contig: "S#1#CM1".into(), strand: 1, start: 0, allele: "A".into() };
        assert!(mapped < unmapped); // Mapped variant orders first
    }

    #[test]
    fn same_locus_different_allele_keys_differ() {
        let rl = |s: &str| Some(RefLoc { contig: "chr2".into(), start: 500, strand: Strand::Forward, status: MapStatus::Mapped });
        assert_ne!(row("ACGT", 0, 0, 0, rl("x")).cluster_key(), row("ACGA", 0, 0, 0, rl("x")).cluster_key());
    }
}