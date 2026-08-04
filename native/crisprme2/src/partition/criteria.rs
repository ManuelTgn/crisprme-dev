//! The user-configurable primary/alternative ranking.
//!
//! A [`PrimaryCriteria`] is an ordered list of [`PrimaryField`]s. Direction is
//! *inherent* to each field — edit-family fields rank ascending ("lowest wins",
//! the original spec), score fields rank descending (higher score = more
//! primary, i.e. the representative alignment is the most concerning predicted
//! cut). [`primary_cmp`] walks the fields in order and, once they are exhausted
//! or all tie, applies the fixed, non-negotiable tiebreak
//! `(strand, target_aligned)`. The result is always a total, run-to-run
//! reproducible order — necessary because the intermediate report is written
//! unordered by parallel sink workers.

use std::cmp::Ordering;

use super::AlignmentRow;

/// One of the four report scores, indexing the [`AlignmentRow::scores`] slot.
#[repr(usize)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreKind {
    Cfd = 0,
    Crista = 1,
    Elevation = 2,
    Aggregate = 3,
}

/// "More primary" direction for a field. Internal detail: it is derived from the
/// field (see [`PrimaryField::direction`]), not supplied by the user.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Ascending,
    Descending,
}

/// A single field the user may rank on. Accepted spec tokens (case-insensitive):
/// `edit-dist`, `bdna`, `brna`, `mm`, `cfd`, `crista`, `elevation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimaryField {
    EditDistance,
    DnaBulges,
    RnaBulges,
    Mismatches,
    Score(ScoreKind),
}

impl PrimaryField {
    /// Numeric value of this field for `row`, as `f64`.
    /// Returns `NaN` for an `NA` (uncomputed) score.
    #[inline]
    fn value(self, row: &AlignmentRow) -> f64 {
        match self {
            PrimaryField::EditDistance => row.edit_distance() as f64,
            PrimaryField::DnaBulges => row.bdna as f64,
            PrimaryField::RnaBulges => row.brna as f64,
            PrimaryField::Mismatches => row.mm as f64,
            PrimaryField::Score(k) => row.scores[k as usize] as f64,
        }
    }

    /// Inherent direction: minimise the edit family, maximise scores.
    #[inline]
    fn direction(self) -> Direction {
        match self {
            PrimaryField::Score(_) => Direction::Descending,
            _ => Direction::Ascending,
        }
    }

    /// Compare `a` vs `b` on this single field. `Ordering::Less` means "a is
    /// *more* primary". An `NA` (NaN) value is always least primary, whatever
    /// the direction — a row with no computed score never beats one that has it.
    #[inline]
    fn cmp(self, a: &AlignmentRow, b: &AlignmentRow) -> Ordering {
        let (va, vb) = (self.value(a), self.value(b));
        match (va.is_nan(), vb.is_nan()) {
            (true, true) => Ordering::Equal,
            (true, false) => Ordering::Greater, // a is NA -> least primary
            (false, true) => Ordering::Less,    // b is NA -> a wins
            (false, false) => {
                let ord = va.partial_cmp(&vb).expect("finite values here");
                match self.direction() {
                    Direction::Ascending => ord,
                    Direction::Descending => ord.reverse(),
                }
            }
        }
    }

    /// Parse a criteria token (case-insensitive). The accepted vocabulary
    /// mirrors the Python `CRITERIA` list.
    fn from_token(s: &str) -> Result<Self, String> {
        Ok(match s.to_ascii_lowercase().as_str() {
            "edit-dist" => PrimaryField::EditDistance,
            "bdna" => PrimaryField::DnaBulges,
            "brna" => PrimaryField::RnaBulges,
            "mm" => PrimaryField::Mismatches,
            "cfd" => PrimaryField::Score(ScoreKind::Cfd),
            "crista" => PrimaryField::Score(ScoreKind::Crista),
            "elevation" => PrimaryField::Score(ScoreKind::Elevation),
            other => {
                return Err(format!(
                    "unknown criteria field {other:?} (expected one of: \
                     edit-dist, bdna, brna, mm, cfd, crista, elevation)"
                ))
            }
        })
    }
}

/// User-configurable ordering used to pick the primary within a cluster.
///
/// `order` is the ordered list of fields (user-settable). The fixed
/// `(strand, target_aligned)` tiebreak is *not* part of `order`; it is applied
/// by [`primary_cmp`].
#[derive(Debug, Clone)]
pub struct PrimaryCriteria {
    pub order: Vec<PrimaryField>,
}

impl Default for PrimaryCriteria {
    /// The original spec: lowest edit distance, then bdna, then brna, then mm.
    ///
    /// Because `edit = mm + bdna + brna`, the trailing `Mismatches` field is
    /// inert under this exact order (once edit, bdna, brna tie, mm is forced
    /// equal). It is kept for fidelity to the spec and becomes active as soon as
    /// the user reorders the fields.
    fn default() -> Self {
        use PrimaryField::*;
        Self {
            order: vec![EditDistance, DnaBulges, RnaBulges, Mismatches],
        }
    }
}

impl PrimaryCriteria {
    /// Build criteria from an ordered list of field tokens (e.g. handed in from
    /// the CLI across the FFI). The fixed `(strand, target_aligned)` tiebreak is
    /// applied by [`primary_cmp`] and is never part of the list. An empty list
    /// is rejected. `Err` carries a user-facing message that the FFI layer
    /// surfaces as a `ValueError`.
    pub fn from_spec(fields: &[String]) -> Result<Self, String> {
        if fields.is_empty() {
            return Err("criteria list is empty; provide at least one field".into());
        }
        let order = fields
            .iter()
            .map(|f| PrimaryField::from_token(f))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { order })
    }
}

/// Total order over rows: user `criteria` first, then the fixed tiebreak
/// `(strand, target_aligned)`. `Ordering::Less` means "a should be the primary
/// over b".
///
/// Strand is constant within a strand-aware cluster, so in practice the strand
/// leg of the tiebreak is inert and selection falls through to `target_aligned`.
pub fn primary_cmp(a: &AlignmentRow, b: &AlignmentRow, criteria: &PrimaryCriteria) -> Ordering {
    for field in &criteria.order {
        let ord = field.cmp(a, b);
        if ord != Ordering::Equal {
            return ord;
        }
    }
    a.strand
        .as_bit()
        .cmp(&b.strand.as_bit())
        .then_with(|| a.target_aligned.cmp(&b.target_aligned))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::occurence::Strand;

    fn row(strand: Strand, mm: u8, bdna: u8, brna: u8) -> AlignmentRow {
        AlignmentRow {
            contig: "1".into(),
            strand,
            start: 100,
            mm,
            bdna,
            brna,
            scores: [f32::NAN; 4],
            target_aligned: "AAAA".into(),
        }
    }
    fn spec(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn default_prefers_lowest_edit() {
        let c = PrimaryCriteria::default();
        assert_eq!(
            primary_cmp(
                &row(Strand::Forward, 2, 0, 0),
                &row(Strand::Forward, 3, 0, 0),
                &c
            ),
            Ordering::Less
        );
    }

    #[test]
    fn default_edit_tie_breaks_on_bdna() {
        let c = PrimaryCriteria::default();
        // both edit 3; lower bdna wins
        assert_eq!(
            primary_cmp(
                &row(Strand::Forward, 1, 2, 0),
                &row(Strand::Forward, 3, 0, 0),
                &c
            ),
            Ordering::Greater
        );
    }

    #[test]
    fn mm_inert_under_default_falls_to_tiebreak() {
        let c = PrimaryCriteria::default();
        let mut a = row(Strand::Forward, 2, 1, 1);
        let mut b = row(Strand::Forward, 2, 1, 1);
        a.target_aligned = "AAAA".into();
        b.target_aligned = "AAAT".into();
        assert_eq!(primary_cmp(&a, &b, &c), Ordering::Less);
    }

    #[test]
    fn score_field_is_descending() {
        let c = PrimaryCriteria::from_spec(&spec(&["cfd"])).unwrap();
        let mut a = row(Strand::Forward, 0, 0, 0);
        let mut b = row(Strand::Forward, 0, 0, 0);
        a.scores[ScoreKind::Cfd as usize] = 0.9;
        b.scores[ScoreKind::Cfd as usize] = 0.3;
        assert_eq!(primary_cmp(&a, &b, &c), Ordering::Less); // higher cfd = primary
    }

    #[test]
    fn na_score_least_primary() {
        let c = PrimaryCriteria::from_spec(&spec(&["cfd"])).unwrap();
        let mut a = row(Strand::Forward, 0, 0, 0);
        let b = row(Strand::Forward, 0, 0, 0); // b's cfd stays NaN
        a.scores[ScoreKind::Cfd as usize] = 0.5;
        assert_eq!(primary_cmp(&a, &b, &c), Ordering::Less);
    }

    #[test]
    fn reordering_makes_mm_active() {
        let c = PrimaryCriteria::from_spec(&spec(&["mm"])).unwrap();
        assert_eq!(
            primary_cmp(
                &row(Strand::Forward, 1, 5, 0),
                &row(Strand::Forward, 4, 0, 0),
                &c
            ),
            Ordering::Less
        );
    }

    #[test]
    fn tiebreak_strand_then_target() {
        let c = PrimaryCriteria::default();
        assert_eq!(
            primary_cmp(
                &row(Strand::Reverse, 1, 0, 0),
                &row(Strand::Forward, 1, 0, 0),
                &c
            ),
            Ordering::Less
        );
        let mut x = row(Strand::Forward, 1, 0, 0);
        let mut y = row(Strand::Forward, 1, 0, 0);
        x.target_aligned = "AAAA".into();
        y.target_aligned = "TTTT".into();
        assert_eq!(primary_cmp(&x, &y, &c), Ordering::Less);
    }

    #[test]
    fn from_spec_case_insensitive() {
        let c = PrimaryCriteria::from_spec(&spec(&["EDIT-DIST", "BDNA", "brna", "Mm"])).unwrap();
        assert!(matches!(c.order[0], PrimaryField::EditDistance));
        assert!(matches!(c.order[3], PrimaryField::Mismatches));
    }

    #[test]
    fn from_spec_rejects_bad_and_empty() {
        assert!(PrimaryCriteria::from_spec(&spec(&["edit_distance"])).is_err()); // old token
        assert!(PrimaryCriteria::from_spec(&spec(&["aggregate"])).is_err()); // not in vocabulary
        assert!(PrimaryCriteria::from_spec(&[]).is_err());
    }

    #[test]
    fn from_spec_reproduces_default() {
        let via = PrimaryCriteria::from_spec(&spec(&["edit-dist", "bdna", "brna", "mm"])).unwrap();
        let a = row(Strand::Forward, 1, 2, 0);
        let b = row(Strand::Forward, 3, 0, 0);
        assert_eq!(
            primary_cmp(&a, &b, &via),
            primary_cmp(&a, &b, &PrimaryCriteria::default())
        );
    }
}
