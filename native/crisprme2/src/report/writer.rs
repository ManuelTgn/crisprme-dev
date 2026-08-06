//! Final assembly report writer (D5).
//!
//! Emits merged rows in the reference-report column shape: `chromosome`,
//! `start`, `strand` carry the **hg38** coordinates for mapped rows (native
//! authoritative internally, hg38 surfaced), and the **native** assembly
//! coordinates for `unmapped` (assembly-specific) rows — `map_status`, the
//! trailing column, tells the two apart. Two columns are appended after the
//! reference `REPORT_HEADER`: `samples` (per-sample copy-presence, e.g.
//! `HG1:1|1,HG2:1|0`) then `map_status` (`mapped` | `ambiguous` | `unmapped`).

use std::io::Write;

use crate::error::crisprme_errors::ReportError;
use crate::liftover::mapper::MapStatus;
use crate::report::row::MergedRow;
use crate::report::samples::{SampleSetRegistry, SampleTable};

/// Trailing columns appended to the reference `REPORT_HEADER`, in order.
pub const ASSEMBLY_TRAILING_HEADER: [&str; 2] = ["samples", "map_status"];

/// Render one merged row as a TSV line into `out` (no trailing newline).
///
/// Column order matches the reference report; `chromosome`/`start`/`strand`
/// hold hg38 coordinates when the row mapped, native coordinates otherwise.
/// `NA` renders for uncomputed scores. The `samples` and `map_status` columns
/// trail.
fn render_row(
    row: &MergedRow,
    registry: &SampleSetRegistry,
    table: &SampleTable,
    out: &mut String,
) {
    // coordinate block: hg38 when mapped, native (assembly-specific) otherwise
    let (contig, start, strand, status): (&str, u32, &str, &str) = match &row.reference {
        Some(r) => (
            r.contig.as_str(),
            r.start,
            r.strand.as_str(),
            match r.status {
                MapStatus::Ambiguous => "ambiguous",
                _ => "mapped",
            },
        ),
        None => (
            row.native.contig.as_str(),
            row.native.start,
            row.native.strand.as_str(),
            "unmapped",
        ),
    };

    // --- reference REPORT_HEADER column order ---
    out.push_str(contig);
    out.push('\t');
    push_u32(out, start);
    out.push('\t');
    out.push_str(strand);
    out.push('\t');
    out.push_str(&row.guide_aligned);
    out.push('\t');
    out.push_str(&row.target_aligned);
    out.push('\t');
    push_u8(out, row.mm);
    out.push('\t');
    push_u8(out, row.bdna);
    out.push('\t');
    push_u8(out, row.brna);
    out.push('\t');
    push_u16(out, row.edit_distance());
    out.push('\t');
    out.push_str(row.bulge_type());
    out.push('\t');
    push_score(out, row.scores.cfd);
    out.push('\t');
    push_score(out, row.scores.crista);
    out.push('\t');
    push_score(out, row.scores.elevation);
    out.push('\t');
    push_score(out, row.scores.aggregate);
    // annotation column(s) intentionally omitted for now (D4 deferred)

    // --- trailing assembly columns ---
    out.push('\t');
    registry.render(row.sample_set, table, out); // samples: HG1:1|1,HG2:1|0
    out.push('\t');
    out.push_str(status); // map_status
}

/// Write the full report: header + one line per merged row.
pub fn write_report<W: Write>(
    mut writer: W,
    header: &str,
    rows: &[MergedRow],
    registry: &SampleSetRegistry,
    table: &SampleTable,
) -> Result<(), ReportError> {
    writeln!(writer, "{header}\t{}", ASSEMBLY_TRAILING_HEADER.join("\t"))
        .map_err(|e| ReportError::io("writing report header", e))?;
    let mut line = String::with_capacity(256);
    for row in rows {
        line.clear();
        render_row(row, registry, table, &mut line);
        writeln!(writer, "{line}").map_err(|e| ReportError::io("writing report row", e))?;
    }
    Ok(())
}

#[inline]
fn push_u8(out: &mut String, v: u8) {
    let mut b = itoa::Buffer::new();
    out.push_str(b.format(v));
}
#[inline]
fn push_u16(out: &mut String, v: u16) {
    let mut b = itoa::Buffer::new();
    out.push_str(b.format(v));
}
#[inline]
fn push_u32(out: &mut String, v: u32) {
    let mut b = itoa::Buffer::new();
    out.push_str(b.format(v));
}
/// `NaN` -> `NA`; finite -> its decimal form (matches the reference sink).
#[inline]
fn push_score(out: &mut String, v: f32) {
    if v.is_nan() {
        out.push_str("NA");
    } else {
        let mut b = ryu::Buffer::new();
        out.push_str(b.format(v));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::occurence::Strand;
    use crate::report::row::{NativeLoc, RefLoc, Scores};
    use crate::report::samples::{SamplePresence, SampleSetRegistry, SampleTable};

    fn table() -> SampleTable {
        SampleTable::new(vec!["HG1".into(), "HG2".into()], vec![vec![1, 2], vec![1, 2]])
    }

    fn mapped_row(reg: &mut SampleSetRegistry) -> MergedRow {
        let ss = reg.intern(&[
            SamplePresence { sample: 0, mask: 0b11 }, // HG1 hom
            SamplePresence { sample: 1, mask: 0b01 }, // HG2 het
        ]);
        MergedRow {
            guide_aligned: "GUIDE".into(),
            target_aligned: "ACGTACGT".into(),
            mm: 1, bdna: 0, brna: 0,
            native: NativeLoc { contig: "HG1#1#c".into(), start: 100, strand: Strand::Forward },
            reference: Some(RefLoc {
                contig: "chr2".into(), start: 500, end: 523,
                strand: Strand::Reverse, status: MapStatus::Mapped,
            }),
            scores: Scores { cfd: 0.5, crista: f32::NAN, elevation: f32::NAN, aggregate: 0.5 },
            annotation: 0,
            sample_set: ss,
        }
    }

    #[test]
    fn mapped_row_uses_hg38_coords_and_trailing_cols() {
        let mut reg = SampleSetRegistry::new();
        let row = mapped_row(&mut reg);
        let mut s = String::new();
        render_row(&row, &reg, &table(), &mut s);
        let f: Vec<&str> = s.split('\t').collect();
        assert_eq!(f[0], "chr2");            // hg38 contig replaces chromosome
        assert_eq!(f[1], "500");             // hg38 start replaces start
        assert_eq!(f[2], "-");               // resolved hg38 strand
        assert_eq!(f[f.len() - 2], "HG1:1|1,HG2:1|0"); // samples
        assert_eq!(f[f.len() - 1], "mapped");          // map_status last
    }

    #[test]
    fn unmapped_row_uses_native_coords_and_flag() {
        let mut reg = SampleSetRegistry::new();
        let ss = reg.intern(&[SamplePresence { sample: 0, mask: 0b01 }]);
        let mut row = mapped_row(&mut reg);
        row.reference = None;
        row.sample_set = ss;
        let mut s = String::new();
        render_row(&row, &reg, &table(), &mut s);
        let f: Vec<&str> = s.split('\t').collect();
        assert_eq!(f[0], "HG1#1#c");         // native contig surfaces
        assert_eq!(f[1], "100");
        assert_eq!(f[2], "+");
        assert_eq!(f[f.len() - 1], "unmapped");
    }
}