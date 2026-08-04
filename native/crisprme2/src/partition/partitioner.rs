//! Split the intermediate report into `primary` and `alternative`.
//!
//! In-memory first cut: reads the whole intermediate into RAM, parses each
//! line's decision fields, sorts by the fixed cluster key, sweeps into clusters,
//! and routes each **original line verbatim** (byte-identical, annotation tail
//! included) to one of the two outputs. The sort makes clustering correct across
//! the chunk/batch boundaries the parallel sink can't respect.
//!
//! Scale follow-ups (deferred): an external / spilling sort and zero-copy
//! borrowed rows for whole-genome runs whose intermediate exceeds RAM.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::Path;

use super::{primary_cmp, same_cluster, AlignmentRow, PrimaryCriteria};

use crate::error::crisprme_errors::{PartitionError, FIXED_COLS, MAX_ANNOTATIONS};
use crate::model::occurence::Strand;

// Fixed-schema column indices into every report line.
const COL_CHROM: usize = 0;
const COL_START: usize = 1;
const COL_STRAND: usize = 2;
const COL_TARGET: usize = 4;
const COL_MM: usize = 5;
const COL_BDNA: usize = 6;
const COL_BRNA: usize = 7;
const COL_CFD: usize = 10; // CFD, CRISTA, Elevation, aggregate are contiguous

/// Counts returned for logging after a successful partition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PartitionStats {
    pub clusters: usize,
    pub primary: usize,
    pub alternative: usize,
}

/// A byte range into the file buffer, *including* the trailing `\n`.
#[derive(Clone, Copy)]
struct Span {
    start: usize,
    end: usize,
}

impl Span {
    #[inline]
    fn bytes<'a>(&self, buf: &'a [u8]) -> &'a [u8] {
        &buf[self.start..self.end]
    }
}

/// A parsed row paired with the span of its verbatim source line.
struct ParsedRow {
    row: AlignmentRow,
    span: Span,
}

/// Split a buffer into line spans, each including its trailing `\n`
/// (the final line is included even if unterminated).
fn line_spans(buf: &[u8]) -> Vec<Span> {
    let mut spans = Vec::new();
    let mut start = 0;
    for (i, &b) in buf.iter().enumerate() {
        if b == b'\n' {
            spans.push(Span { start, end: i + 1 });
            start = i + 1;
        }
    }
    if start < buf.len() {
        spans.push(Span {
            start,
            end: buf.len(),
        });
    }
    spans
}

/// Tab-separated fields of a line span, with the terminator trimmed.
fn fields<'a>(buf: &'a [u8], span: Span, line_no: usize) -> Result<Vec<&'a str>, PartitionError> {
    let s = std::str::from_utf8(span.bytes(buf))
        .map_err(|_| PartitionError::NotUtf8 { line: line_no })?;
    Ok(s.trim_end_matches(['\n', '\r']).split('\t').collect())
}

fn parse_u8(v: &str, line: usize, col: usize) -> Result<u8, PartitionError> {
    v.parse().map_err(|_| PartitionError::BadField {
        line,
        col,
        value: v.into(),
        what: "integer",
    })
}

fn parse_u32(v: &str, line: usize, col: usize) -> Result<u32, PartitionError> {
    v.parse().map_err(|_| PartitionError::BadField {
        line,
        col,
        value: v.into(),
        what: "position",
    })
}

fn parse_strand(v: &str, line: usize, col: usize) -> Result<Strand, PartitionError> {
    match v {
        "+" => Ok(Strand::Forward),
        "-" => Ok(Strand::Reverse),
        _ => Err(PartitionError::BadField {
            line,
            col,
            value: v.into(),
            what: "strand",
        }),
    }
}

/// Report scores render as `{:.2}` or `NA`; the latter becomes `NaN`.
fn parse_score(v: &str, line: usize, col: usize) -> Result<f32, PartitionError> {
    if v == "NA" {
        return Ok(f32::NAN);
    }
    v.parse().map_err(|_| PartitionError::BadField {
        line,
        col,
        value: v.into(),
        what: "score",
    })
}

fn parse_row(
    buf: &[u8],
    span: Span,
    line_no: usize,
    ncols: usize,
) -> Result<AlignmentRow, PartitionError> {
    let f = fields(buf, span, line_no)?;
    if f.len() != ncols {
        return Err(PartitionError::FieldCount {
            line: line_no,
            expected: ncols,
            found: f.len(),
        });
    }
    let mut scores = [f32::NAN; 4];
    for (i, slot) in scores.iter_mut().enumerate() {
        *slot = parse_score(f[COL_CFD + i], line_no, COL_CFD + i)?;
    }
    Ok(AlignmentRow {
        contig: f[COL_CHROM].into(),
        start: parse_u32(f[COL_START], line_no, COL_START)?,
        strand: parse_strand(f[COL_STRAND], line_no, COL_STRAND)?,
        mm: parse_u8(f[COL_MM], line_no, COL_MM)?,
        bdna: parse_u8(f[COL_BDNA], line_no, COL_BDNA)?,
        brna: parse_u8(f[COL_BRNA], line_no, COL_BRNA)?,
        scores,
        target_aligned: f[COL_TARGET].into(),
    })
}

/// Sweep span-paired rows (pre-sorted by the cluster key) into cluster runs.
///
/// Reuses the tested [`same_cluster`] predicate; the `&[AlignmentRow]` wrappers
/// [`super::cluster_runs`] / [`super::select_primary`] are the equivalent forms
/// for callers that hold a bare row slice rather than span-paired rows.
fn for_each_cluster(rows: &[ParsedRow], n: u32, mut f: impl FnMut(Range<usize>)) {
    let mut start = 0;
    while start < rows.len() {
        let mut end = start + 1;
        while end < rows.len() && same_cluster(&rows[end - 1].row, &rows[end].row, n) {
            end += 1;
        }
        f(start..end);
        start = end;
    }
}

fn primary_index(rows: &[ParsedRow], run: Range<usize>, criteria: &PrimaryCriteria) -> usize {
    run.min_by(|&i, &j| primary_cmp(&rows[i].row, &rows[j].row, criteria))
        .expect("cluster runs are non-empty by construction")
}

/// Partition `intermediate` into `primary_out` and `alternative_out`.
///
/// `n` is the clustering window; `criteria` selects the primary within each
/// cluster (its fixed `(strand, target_aligned)` tiebreak is applied last).
/// Output lines are byte-identical to the intermediate; each file carries the
/// original header, and rows appear in `(contig, strand, start)` order.
pub fn partition_report(
    intermediate: &Path,
    primary_out: &Path,
    alternative_out: &Path,
    n: u32,
    criteria: &PrimaryCriteria,
) -> Result<PartitionStats, PartitionError> {
    let buf = fs::read(intermediate)
        .map_err(|e| PartitionError::io(format!("reading {}", intermediate.display()), e))?;
    let spans = line_spans(&buf);

    let header = *spans.first().ok_or(PartitionError::MissingHeader)?;
    let ncols = fields(&buf, header, 1)?.len();
    if !(FIXED_COLS..=FIXED_COLS + MAX_ANNOTATIONS).contains(&ncols) {
        return Err(PartitionError::BadHeaderShape { cols: ncols });
    }

    // Parse data rows (header is line 1; data starts at line 2).
    let mut rows = Vec::with_capacity(spans.len().saturating_sub(1));
    for (i, &span) in spans.iter().enumerate().skip(1) {
        let row = parse_row(&buf, span, i + 1, ncols)?;
        rows.push(ParsedRow { row, span });
    }

    // Sort by the fixed cluster key, then by full line bytes so that rows which
    // tie on every field the comparator reads still resolve to a deterministic,
    // run-to-run reproducible primary.
    rows.sort_unstable_by(|a, b| {
        a.row
            .contig
            .cmp(&b.row.contig)
            .then_with(|| a.row.strand.as_bit().cmp(&b.row.strand.as_bit()))
            .then_with(|| a.row.start.cmp(&b.row.start))
            .then_with(|| a.span.bytes(&buf).cmp(b.span.bytes(&buf)))
    });

    let mut prim = BufWriter::new(
        File::create(primary_out)
            .map_err(|e| PartitionError::io(format!("creating {}", primary_out.display()), e))?,
    );
    let mut alt =
        BufWriter::new(File::create(alternative_out).map_err(|e| {
            PartitionError::io(format!("creating {}", alternative_out.display()), e)
        })?);

    let hbytes = header.bytes(&buf);
    prim.write_all(hbytes)
        .map_err(|e| PartitionError::io("writing primary header", e))?;
    alt.write_all(hbytes)
        .map_err(|e| PartitionError::io("writing alternative header", e))?;

    let mut stats = PartitionStats {
        clusters: 0,
        primary: 0,
        alternative: 0,
    };
    let mut io_err: Option<PartitionError> = None;

    for_each_cluster(&rows, n, |run| {
        if io_err.is_some() {
            return;
        }
        stats.clusters += 1;
        let p = primary_index(&rows, run.clone(), criteria);
        for i in run {
            let is_primary = i == p;
            let target = if is_primary { &mut prim } else { &mut alt };
            if let Err(e) = target.write_all(rows[i].span.bytes(&buf)) {
                io_err = Some(PartitionError::io("writing report row", e));
                return;
            }
            if is_primary {
                stats.primary += 1;
            } else {
                stats.alternative += 1;
            }
        }
    });

    if let Some(e) = io_err {
        return Err(e);
    }
    prim.flush()
        .map_err(|e| PartitionError::io("flushing primary report", e))?;
    alt.flush()
        .map_err(|e| PartitionError::io("flushing alternative report", e))?;

    tracing::info!(
        "partitioned report: {} clusters -> {} primary, {} alternative",
        stats.clusters,
        stats.primary,
        stats.alternative
    );
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "chromosome\tstart\tstrand\tsgRNA_aligned\ttarget_aligned\tmismatches\tdna_bulges\trna_bulges\tbulge_type\tedit_distance\tCFD_score\tCRISTA_score\tElevation_score\taggregate_score\tgene\tregulatory";

    // Crafted for: single-linkage chaining, strand split, contig split, tie-on-bdna.
    const R1: &str = "chr1\t100\t+\tGUIDE\tTARGETA\t2\t0\t0\tX\t2\tNA\tNA\tNA\tNA\tEX1\tNA";
    const R2: &str = "chr1\t102\t+\tGUIDE\tTARGETB\t1\t1\t0\tDNA\t2\tNA\tNA\tNA\tNA\tNA\tNA";
    const R3: &str = "chr1\t101\t-\tGUIDE\tTARGETC\t0\t0\t0\tX\t0\tNA\tNA\tNA\tNA\tNA\tREG1";
    const R4: &str = "chr1\t200\t+\tGUIDE\tTARGETD\t3\t0\t0\tX\t3\tNA\tNA\tNA\tNA\tNA\tNA";
    const R5: &str = "chr2\t100\t+\tGUIDE\tTARGETE\t1\t0\t0\tX\t1\tNA\tNA\tNA\tNA\tNA\tNA";

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("crisprme_part_{}_{}", std::process::id(), name))
    }
    fn write(path: &Path, lines: &[&str]) {
        let mut s = String::new();
        for l in lines {
            s.push_str(l);
            s.push('\n');
        }
        fs::write(path, s).unwrap();
    }

    #[test]
    fn end_to_end_split() {
        let inp = tmp("in.tsv");
        let p = tmp("prim.tsv");
        let a = tmp("alt.tsv");
        write(&inp, &[HEADER, R4, R2, R5, R1, R3]); // scrambled

        let stats = partition_report(&inp, &p, &a, 3, &PrimaryCriteria::default()).unwrap();
        assert_eq!(
            stats,
            PartitionStats {
                clusters: 4,
                primary: 4,
                alternative: 1
            }
        );

        assert_eq!(
            fs::read_to_string(&p).unwrap(),
            format!("{HEADER}\n{R3}\n{R1}\n{R4}\n{R5}\n")
        );
        assert_eq!(fs::read_to_string(&a).unwrap(), format!("{HEADER}\n{R2}\n"));
    }

    #[test]
    fn byte_identical_lines_preserved() {
        let inp = tmp("in2.tsv");
        let p = tmp("prim2.tsv");
        let a = tmp("alt2.tsv");
        write(&inp, &[HEADER, R1, R2]);
        partition_report(&inp, &p, &a, 3, &PrimaryCriteria::default()).unwrap();
        assert!(fs::read_to_string(&p).unwrap().contains(&format!("{R1}\n"))); // annotation tail intact
        assert!(fs::read_to_string(&a).unwrap().contains(&format!("{R2}\n")));
    }

    #[test]
    fn header_only_is_ok() {
        let inp = tmp("in3.tsv");
        let p = tmp("prim3.tsv");
        let a = tmp("alt3.tsv");
        write(&inp, &[HEADER]);
        let stats = partition_report(&inp, &p, &a, 3, &PrimaryCriteria::default()).unwrap();
        assert_eq!(
            stats,
            PartitionStats {
                clusters: 0,
                primary: 0,
                alternative: 0
            }
        );
        assert_eq!(fs::read_to_string(&p).unwrap(), format!("{HEADER}\n"));
    }

    #[test]
    fn bad_field_count_errors() {
        let inp = tmp("in4.tsv");
        let p = tmp("prim4.tsv");
        let a = tmp("alt4.tsv");
        write(&inp, &[HEADER, "chr1\t100\t+\tSHORT"]);
        let e = partition_report(&inp, &p, &a, 3, &PrimaryCriteria::default()).unwrap_err();
        assert!(matches!(
            e,
            PartitionError::FieldCount {
                line: 2,
                expected: 16,
                found: 4
            }
        ));
    }
}
