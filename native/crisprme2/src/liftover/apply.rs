//! Apply liftover to a per-haplotype intermediate report.
//!
//! Reads a native-coordinate TSV report, lifts each row's leftmost-forward
//! anchor (the `start` column) through a parsed chain, and appends the hg38
//! columns additively — native coordinates stay authoritative and every row is
//! preserved (unmapped rows kept with `NA`/`unmapped`). Std-only; the PyO3
//! wrapper and typed Python error live in the bindings module.

use std::io::{BufRead, Write};

use crate::error::crisprme_errors::LiftoverError;
use crate::liftover::mapper::{LiftOver, LiftResult, MapStatus};

/// Columns appended to the report by liftover, in order.
pub const LIFTOVER_HEADER: [&str; 4] = ["hg38_chr", "hg38_pos", "hg38_flipped", "map_status"];

#[derive(Debug, Default, Clone, Copy)]
pub struct LiftReportStats {
    pub mapped: u64,
    pub ambiguous: u64,
    pub unmapped: u64,
}

/// PanSN `sample#hap#contig` -> bare contig (the chain `t_name`). Defensive:
/// returns the input unchanged when there is no `#`.
#[inline]
fn normalize_contig(pansn: &str) -> &str {
    pansn.rsplit('#').next().unwrap_or(pansn)
}

/// Lift every row of a report from `reader` to `writer`, appending the four
/// [`LIFTOVER_HEADER`] columns. `contig_col`/`start_col` are looked up by name
/// in the header, so column order is not assumed.
pub fn lift_report_core<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    lifter: &LiftOver,
    contig_col: &str,
    start_col: &str,
) -> Result<LiftReportStats, LiftoverError> {
    let mut lines = reader.lines();
    let header = match lines.next() {
        Some(r) => r.map_err(|e| LiftoverError::io("reading report header", e))?,
        None => return Err(LiftoverError::MissingHeader),
    };
    let cols: Vec<&str> = header.split('\t').collect();
    let ci = cols.iter().position(|&c| c == contig_col)
        .ok_or_else(|| LiftoverError::MissingColumn(contig_col.into()))?;
    let si = cols.iter().position(|&c| c == start_col)
        .ok_or_else(|| LiftoverError::MissingColumn(start_col.into()))?;

    writeln!(writer, "{}\t{}", header, LIFTOVER_HEADER.join("\t"))
        .map_err(|e| LiftoverError::io("writing lifted report", e))?;

    let mut stats = LiftReportStats::default();
    for (idx, line_res) in lines.enumerate() {
        let line = line_res.map_err(|e| LiftoverError::io("reading report", e))?;
        if line.is_empty() {
            continue;
        }
        let line_no = idx + 2; // 1-based, +1 for the header
        let fields: Vec<&str> = line.split('\t').collect();
        let contig = *fields.get(ci)
            .ok_or_else(|| LiftoverError::bad_row(line_no, "missing contig field"))?;
        let start_s = *fields.get(si)
            .ok_or_else(|| LiftoverError::bad_row(line_no, "missing start field"))?;
        let start: u64 = start_s.parse()
            .map_err(|_| LiftoverError::bad_row(line_no, format!("invalid start {start_s:?}")))?;

        let (chr, pos, flip, status) = match lifter.lift(normalize_contig(contig), start) {
            LiftResult::Ok(l) => {
                let status = match l.status {
                    MapStatus::Ambiguous => { stats.ambiguous += 1; "ambiguous" }
                    _ => { stats.mapped += 1; "mapped" }
                };
                (l.ref_name, l.ref_pos.to_string(), if l.flipped { "1" } else { "0" }, status)
            }
            LiftResult::NoMap(_) => {
                stats.unmapped += 1;
                ("NA".to_string(), "NA".to_string(), "NA", "unmapped")
            }
        };
        writeln!(writer, "{line}\t{chr}\t{pos}\t{flip}\t{status}")
            .map_err(|e| LiftoverError::io("writing lifted report", e))?;
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::liftover::chain::ChainFile;
    use crate::liftover::mapper::LiftOver;
    use std::io::Cursor;

    #[test]
    fn normalizes_pansn() {
        assert_eq!(normalize_contig("HG06807#2#CM101429.1"), "CM101429.1");
        assert_eq!(normalize_contig("CM101429.1"), "CM101429.1");
    }

    #[test]
    fn lifts_rows_and_keeps_unmapped() {
        // chain: target CM1 [0,5) -> ref chr2 [100,105)
        let cf = ChainFile::parse(Cursor::new(
            "chain 100 CM1 20 + 0 5 chr2 300 + 100 105 1\n5\n",
        )).unwrap();
        let lifter = LiftOver::new(&cf, 10);
        let report = "chromosome\tstart\tstrand\nS#1#CM1\t2\t+\nS#1#CM1\t9\t-\n";
        let mut out = Vec::new();
        let stats = lift_report_core(Cursor::new(report), &mut out, &lifter, "chromosome", "start").unwrap();
        let out = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "chromosome\tstart\tstrand\thg38_chr\thg38_pos\thg38_flipped\tmap_status");
        assert_eq!(lines[1], "S#1#CM1\t2\t+\tchr2\t102\t0\tmapped");   // 2 in [0,5) -> 102
        assert_eq!(lines[2], "S#1#CM1\t9\t-\tNA\tNA\tNA\tunmapped");   // 9 outside chain
        assert_eq!((stats.mapped, stats.unmapped), (1, 1));
    }
}