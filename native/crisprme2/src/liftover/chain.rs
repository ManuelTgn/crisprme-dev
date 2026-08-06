//! Standalone UCSC chain-file parser for haplotype -> reference liftover.
//!
//! Parses `.chain` files (as produced by `axtChain`) into an eagerly
//! materialised set of aligned [`Block`]s with absolute, 0-based, half-open,
//! **plus-strand** coordinates on both the target (haplotype) and query
//! (reference) sides.
//!
//! # Orientation
//! In UCSC chain format the *target* `t` is the first sequence (here the
//! haplotype) and the *query* `q` is the second (here GRCh38); blocks map
//! `t -> q`, which is exactly the liftover direction. The target strand is
//! always `+`. When the query strand is `-`, the file's query coordinates are
//! measured on the minus strand (from `qSize`); this parser converts them to
//! plus-strand at parse time, so every [`Block`] stores `q_start < q_end` on
//! the forward reference. The chain's [`Chain::q_strand`] records whether the
//! within-block mapping is orientation-reversed (target `t_start` aligns to
//! query `q_end - 1`), which the C2 mapper uses to place an individual base.

use std::collections::HashMap;
use std::io::BufRead;

use crate::error::crisprme_errors::LiftoverError;

/// Query (reference) strand of a chain. The target (haplotype) side is always
/// forward in UCSC chain format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    fn parse(tok: &str, line: usize) -> Result<Self, LiftoverError> {
        match tok {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            other => Err(LiftoverError::malformed_chain(
                line,
                format!("invalid strand {other:?}"),
            )),
        }
    }
}

/// One ungapped aligned block, absolute and plus-strand on both sides.
///
/// Target `[t_start, t_end)` maps to query `[q_start, q_end)`; both are
/// 0-based, half-open, plus-strand, and equal length. For a reverse chain the
/// mapping *within* the block runs the other way (see [`Chain::q_strand`])
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block {
    pub t_start: u64,
    pub t_end: u64,
    pub q_start: u64,
    pub q_end: u64,
}

impl Block {
    #[inline]
    pub fn t_len(&self) -> u64 {
        self.t_end - self.t_start
    }
}

/// A single chain: header metadata plus its eagerly materialised blocks
/// (sorted ascending by `t_start`)
#[derive(Debug, Clone)]
pub struct Chain {
    pub id: u64,
    pub score: i64,
    pub t_name: String,
    pub t_size: u64,
    pub t_start: u64,
    pub t_end: u64,
    pub q_name: String,
    pub q_size: u64,
    /// Query strand; target strand is always forward.
    pub q_strand: Strand,
    /// Query span in plus-strand coordinates (`q_start < q_end`).
    pub q_start: u64,
    pub q_end: u64,
    pub blocks: Vec<Block>,
}

/// A parsed chain file: all chains plus an index from target name to its
/// chains, ordered by descending score (primary chain first)
#[derive(Debug, Clone)]
pub struct ChainFile {
    pub chains: Vec<Chain>,
    by_tname: HashMap<String, Vec<usize>>,
}

impl ChainFile {
    /// Parse a chain file from any buffered reader.
    pub fn parse<R: BufRead>(reader: R) -> Result<Self, LiftoverError> {
        let mut chains: Vec<Chain> = Vec::new();
        // (header, data lines as (size, Option<(dt, dq)>), header line number)
        let mut pending: Option<(Header, Vec<(u64, Option<(u64, u64)>)>, usize)> = None;

        for (idx, line_res) in reader.lines().enumerate() {
            let line_no = idx + 1;
            let raw = line_res.map_err(|e| LiftoverError::io("reading chain file", e))?;
            let line = raw.trim();

            if line.is_empty() {
                if let Some((h, data, fl)) = pending.take() {
                    chains.push(finalize(h, &data, fl)?);
                }
                continue;
            }
            if line.starts_with('#') {
                continue; // ##matrix / ##gapPenalties and any other comment
            }
            let first = line.split_whitespace().next().unwrap_or("");
            if first == "chain" {
                if let Some((h, data, fl)) = pending.take() {
                    chains.push(finalize(h, &data, fl)?);
                }
                pending = Some((parse_header(line, line_no)?, Vec::new(), line_no));
            } else {
                match pending.as_mut() {
                    Some((_, data, _)) => data.push(parse_data_line(line, line_no)?),
                    None => {
                        return Err(LiftoverError::malformed_chain(
                            line_no,
                            "alignment data line outside of any chain",
                        ))
                    }
                }
            }
        }
        if let Some((h, data, fl)) = pending.take() {
            chains.push(finalize(h, &data, fl)?);
        }

        // index by target name, primary (highest score) chain first
        let mut by_tname: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, c) in chains.iter().enumerate() {
            by_tname.entry(c.t_name.clone()).or_default().push(i);
        }
        for idxs in by_tname.values_mut() {
            idxs.sort_by(|&a, &b| chains[b].score.cmp(&chains[a].score));
        }
        Ok(ChainFile { chains, by_tname })
    }

    /// Open and parse a chain file from a path (plain text).
    pub fn from_path(path: impl AsRef<std::path::Path>) -> Result<Self, LiftoverError> {
        let file = std::fs::File::open(&path).map_err(|e| {
            LiftoverError::io(format!("opening chain file {}", path.as_ref().display()), e)
        })?;
        Self::parse(std::io::BufReader::new(file))
    }

    /// Chains covering target `t_name`, ordered by descending score. Empty if
    /// the target is absent.
    pub fn chains_for<'a>(&'a self, t_name: &str) -> impl Iterator<Item = &'a Chain> + 'a {
        self.by_tname
            .get(t_name)
            .into_iter()
            .flatten()
            .map(move |&i| &self.chains[i])
    }

    /// Distinct target (haplotype) contig names present in the file.
    pub fn target_names(&self) -> impl Iterator<Item = &str> {
        self.by_tname.keys().map(String::as_str)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.chains.len()
    }
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.chains.is_empty()
    }
}

struct Header {
    id: u64,
    score: i64,
    t_name: String,
    t_size: u64,
    t_start: u64,
    t_end: u64,
    q_name: String,
    q_size: u64,
    q_strand: Strand,
    q_start_file: u64, // file frame: minus-strand-relative when q_strand == Reverse
    q_end_file: u64,
}

fn parse_header(line: &str, ln: usize) -> Result<Header, LiftoverError> {
    let f: Vec<&str> = line.split_whitespace().collect();
    if f.len() != 12 && f.len() != 13 {
        return Err(LiftoverError::malformed_chain(
            ln,
            format!("chain header expects 12 or 13 fields, got {}", f.len()),
        ));
    }
    let score = parse_i64(f[1], ln, "score")?;
    let t_name = f[2].to_string();
    let t_size = parse_u64(f[3], ln, "tSize")?;
    if f[4] != "+" {
        return Err(LiftoverError::malformed_chain(
            ln,
            format!("target strand must be '+', got {:?}", f[4]),
        ));
    }
    let t_start = parse_u64(f[5], ln, "tStart")?;
    let t_end = parse_u64(f[6], ln, "tEnd")?;
    let q_name = f[7].to_string();
    let q_size = parse_u64(f[8], ln, "qSize")?;
    let q_strand = Strand::parse(f[9], ln)?;
    let q_start_file = parse_u64(f[10], ln, "qStart")?;
    let q_end_file = parse_u64(f[11], ln, "qEnd")?;
    let id = if f.len() == 13 {
        parse_u64(f[12], ln, "id")?
    } else {
        0
    };

    if t_end < t_start || q_end_file < q_start_file {
        return Err(LiftoverError::malformed_chain(ln, "end < start in header"));
    }
    if t_end > t_size || q_end_file > q_size {
        return Err(LiftoverError::malformed_chain(ln, "end > size in header"));
    }
    Ok(Header {
        id,
        score,
        t_name,
        t_size,
        t_start,
        t_end,
        q_name,
        q_size,
        q_strand,
        q_start_file,
        q_end_file,
    })
}

fn parse_data_line(line: &str, ln: usize) -> Result<(u64, Option<(u64, u64)>), LiftoverError> {
    let f: Vec<&str> = line.split_whitespace().collect();
    match f.len() {
        1 => Ok((parse_u64(f[0], ln, "size")?, None)),
        3 => Ok((
            parse_u64(f[0], ln, "size")?,
            Some((parse_u64(f[1], ln, "dt")?, parse_u64(f[2], ln, "dq")?)),
        )),
        n => Err(LiftoverError::malformed_chain(
            ln,
            format!("alignment data line expects 1 or 3 fields, got {n}"),
        )),
    }
}

fn finalize(
    h: Header,
    data: &[(u64, Option<(u64, u64)>)],
    ln: usize,
) -> Result<Chain, LiftoverError> {
    if data.is_empty() {
        return Err(LiftoverError::malformed_chain(
            ln,
            "chain has no aligned blocks",
        ));
    }
    let last = data.len() - 1;
    let mut blocks = Vec::with_capacity(data.len());
    let mut t = h.t_start;
    let mut q = h.q_start_file;

    for (i, &(size, gap)) in data.iter().enumerate() {
        if size == 0 {
            return Err(LiftoverError::malformed_chain(ln, "block size must be > 0"));
        }
        if gap.is_none() && i != last {
            return Err(LiftoverError::malformed_chain(
                ln,
                "single-field block appears before end of chain",
            ));
        }
        let (q_start, q_end) = match h.q_strand {
            Strand::Forward => (q, q + size),
            // minus-strand file interval [q, q+size) -> plus-strand
            Strand::Reverse => (h.q_size - (q + size), h.q_size - q),
        };
        blocks.push(Block {
            t_start: t,
            t_end: t + size,
            q_start,
            q_end,
        });

        let (dt, dq) = gap.unwrap_or((0, 0));
        t += size + dt;
        q += size + dq;
    }

    // block advances must reproduce the header spans exactly (file frame)
    if t != h.t_end {
        return Err(LiftoverError::malformed_chain(
            ln,
            format!("target blocks span to {t} but header tEnd={}", h.t_end),
        ));
    }
    if q != h.q_end_file {
        return Err(LiftoverError::malformed_chain(
            ln,
            format!("query blocks span to {q} but header qEnd={}", h.q_end_file),
        ));
    }

    let (q_start, q_end) = match h.q_strand {
        Strand::Forward => (h.q_start_file, h.q_end_file),
        Strand::Reverse => (h.q_size - h.q_end_file, h.q_size - h.q_start_file),
    };
    debug_assert!(blocks.windows(2).all(|w| w[0].t_start < w[1].t_start));

    Ok(Chain {
        id: h.id,
        score: h.score,
        t_name: h.t_name,
        t_size: h.t_size,
        t_start: h.t_start,
        t_end: h.t_end,
        q_name: h.q_name,
        q_size: h.q_size,
        q_strand: h.q_strand,
        q_start,
        q_end,
        blocks,
    })
}

fn parse_u64(tok: &str, ln: usize, field: &str) -> Result<u64, LiftoverError> {
    tok.parse::<u64>()
        .map_err(|_| LiftoverError::malformed_chain(ln, format!("invalid {field} value {tok:?}")))
}

fn parse_i64(tok: &str, ln: usize, field: &str) -> Result<i64, LiftoverError> {
    tok.parse::<i64>()
        .map_err(|_| LiftoverError::malformed_chain(ln, format!("invalid {field} value {tok:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn parse_str(s: &str) -> Result<ChainFile, LiftoverError> {
        ChainFile::parse(Cursor::new(s))
    }

    #[test]
    fn forward_chain_blocks_and_spans() {
        // t: [0,5) then +2 gap -> [7,11); q: [0,5) then +3 gap -> [8,12)
        let cf = parse_str("chain 100 T 20 + 0 11 Q 30 + 0 12 1\n5 2 3\n4\n").unwrap();
        let c = &cf.chains[0];
        assert_eq!(c.q_strand, Strand::Forward);
        assert_eq!(
            c.blocks,
            vec![
                Block {
                    t_start: 0,
                    t_end: 5,
                    q_start: 0,
                    q_end: 5
                },
                Block {
                    t_start: 7,
                    t_end: 11,
                    q_start: 8,
                    q_end: 12
                },
            ]
        );
        assert_eq!((c.t_end, c.q_end), (11, 12));
    }

    #[test]
    fn reverse_chain_maps_to_plus_strand() {
        // qSize=30: file [0,5)->plus[25,30); file [8,12)->plus[18,22)
        let cf = parse_str("chain 100 T 20 + 0 11 Q 30 - 0 12 1\n5 2 3\n4\n").unwrap();
        let c = &cf.chains[0];
        assert_eq!(c.q_strand, Strand::Reverse);
        assert_eq!(
            c.blocks,
            vec![
                Block {
                    t_start: 0,
                    t_end: 5,
                    q_start: 25,
                    q_end: 30
                },
                Block {
                    t_start: 7,
                    t_end: 11,
                    q_start: 18,
                    q_end: 22
                },
            ]
        );
        assert_eq!((c.q_start, c.q_end), (18, 30));
    }

    #[test]
    fn parses_real_header_fields_and_comments() {
        let text = "\
##matrix=axtChain 16 91,-114
##gapPenalties=axtChain O=400 E=30
chain 22419049526 CM101396.1 244346724 + 0 5 chr2 242193529 + 0 5 1
5
";
        let c = &parse_str(text).unwrap().chains[0];
        assert_eq!(c.t_name, "CM101396.1");
        assert_eq!(c.q_name, "chr2");
        assert_eq!((c.t_size, c.q_size), (244346724, 242193529));
        assert_eq!(c.score, 22419049526);
        assert_eq!(c.q_strand, Strand::Forward);
        assert_eq!(c.id, 1);
    }

    #[test]
    fn cumulative_mismatch_is_error() {
        let e = parse_str("chain 1 T 20 + 0 99 Q 20 + 0 5 1\n5\n").unwrap_err();
        assert!(matches!(e, LiftoverError::Malformed { .. }));
    }

    #[test]
    fn data_outside_chain_is_error() {
        assert!(matches!(
            parse_str("5 2 3\n").unwrap_err(),
            LiftoverError::Malformed { .. }
        ));
    }

    #[test]
    fn single_field_before_end_is_error() {
        let e = parse_str("chain 1 T 20 + 0 11 Q 20 + 0 11 1\n5\n4 0 0\n").unwrap_err();
        assert!(matches!(e, LiftoverError::Malformed { .. }));
    }

    #[test]
    fn multiple_chains_sorted_by_score_desc() {
        let text = "\
chain 10 T 20 + 0 5 Q 20 + 0 5 1
5

chain 999 T 20 + 5 10 Q 20 + 5 10 2
5
";
        let cf = parse_str(text).unwrap();
        let order: Vec<i64> = cf.chains_for("T").map(|c| c.score).collect();
        assert_eq!(order, vec![999, 10]);
    }
}
