//! Feature registry: the name <-> bit mapping for one annotation BED file.
//!
//! # Role in the annotation pipeline
//!
//! Functional annotation attaches, to every off-target, the set of annotation
//! *terms* it overlaps. Each input BED file becomes exactly one annotation
//! column in the report, and the terms live in the **4th column** of that BED
//! (the BED `name` field). This registry is the single source of truth that
//! maps those terms to bits and back:
//!
//! - **annotate time** — the transform fetches the BED records overlapping a
//!   target, and for each record turns its term into a bit with [`mask_of`] and
//!   ORs it into the alignment's `u32` annotation slot
//!   (`Alignment::features[slot]`).
//! - **sink time** — the writer takes the accumulated `u32` and expands it back
//!   into a comma-separated list of term names with [`decode`].
//!
//! Because there is exactly one bit per term and the slot is a `u32`, a single
//! BED may define at most [`MAX_FEATURES`] (= 32) distinct terms. This is by
//! design: functional-annotation vocabularies are small and bounded (chromatin
//! states, cCRE classes, regulatory-build types, ...). Gene-level annotation,
//! which is *not* bounded, is handled by a separate path and SHOULD NEVER flow
//! through this registry. Exceeding the cap is a hard error, not a silent
//! truncation — see [`AnnotationError::TooManyFeatures`].
//!
//! # Determinism
//!
//! Bits are assigned in **first-seen order** while scanning the file, so a
//! given BED always produces the same term->bit mapping, and [`decode`] always
//! lists a target's terms in a stable order (ascending bit index). The mapping
//! is per file: two BEDs sharing a term name assign it independently.
//!
//! # Concurrency
//!
//! A `FeatureRegistry` is immutable after construction (no interior
//! mutability), so it is `Send + Sync` and is meant to be built once and shared
//! read-only across worker threads behind an `Arc`.
//!
//! # Input format
//!
//! BED files are accepted either plain or bgzip-compressed (`.bed.gz`, the
//! tabix-indexed form pysam expects). Compression is detected by gzip magic
//! bytes, so the extension is irrelevant. BGZF is a series of concatenated
//! gzip members terminated by an empty member, which is why decompression uses
//! [`MultiGzDecoder`] rather than a single-member decoder.

use crate::error::crisprme_errors::AnnotationError;

use flate2::read::MultiGzDecoder;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Maximum number of distinct terms a single BED may define.
///
/// One bit per term in the `u32` annotation slot, hence 32. Widening the slot
/// to `u64` in `model::alignment` would let this grow to 64, but bounded
/// functional-annotation vocabularies sit comfortably under 32.
pub const MAX_FEATURES: usize = 32;

/// Immutable term <-> bit mapping for one annotation BED file.
///
/// See the [module docs](self) for the surrounding pipeline. The registry owns
/// two parallel views of the same small set of terms:
///
/// - `names[bit]` -> term (used by [`decode`](Self::decode) / [`term`](Self::term))
/// - `index[term]` -> bit (used by [`mask_of`](Self::mask_of) / [`bit_of`](Self::bit_of))
#[derive(Debug)]
pub struct FeatureRegistry {
    /// Term name for each bit index; `names.len() <= MAX_FEATURES`.
    names: Box<[Box<str>]>,
    /// Reverse map: term name -> bit index in `0..names.len()`.
    index: HashMap<Box<str>, u32>,
}

impl FeatureRegistry {
    /// Build a registry from a BED file (plain or bgzip-compressed).
    ///
    /// The 4th column of every non-comment line is taken as an annotation term.
    /// Distinct terms are assigned bits in first-seen order. Blank lines and
    /// lines beginning with `#` (tabix meta lines) are skipped.
    ///
    /// # Errors
    /// - [`AnnotationError::IoError`] — the file cannot be opened or read (also
    ///   covers a corrupt gzip stream).
    /// - [`AnnotationError::MalformedBed`] — a data line has fewer than 4
    ///   tab-separated fields, or its 4th field is empty.
    /// - [`AnnotationError::TooManyFeatures`] — more than [`MAX_FEATURES`]
    ///   distinct terms are present.
    pub fn from_bed<P: AsRef<Path>>(path: P) -> Result<Self, AnnotationError> {
        let path = path.as_ref();
        let reader = open_bed(path)?;
        Self::from_reader(reader, &path.display().to_string())
    }

    /// Core builder over any `BufRead`. `source` is only used to label errors.
    ///
    /// Split out from [`from_bed`](Self::from_bed) so the parsing logic can be
    /// unit-tested against in-memory readers without touching the filesystem.
    fn from_reader<R: BufRead>(reader: R, source: &str) -> Result<Self, AnnotationError> {
        let mut names: Vec<Box<str>> = Vec::new();
        let mut index: HashMap<Box<str>, u32> = HashMap::new();

        for (line_no, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| AnnotationError::IoError(e.to_string()))?;
            // `lines()` strips '\n' but not a trailing '\r' on CRLF files, which
            // would otherwise contaminate the term of a BED4 record.
            let line = line.trim_end_matches('\r');

            if line.is_empty() || line.starts_with('#') {
                continue; // blank line or tabix meta/comment line
            }

            // BED layout: chrom, start, end, name(=term), [score, strand, ...].
            // We only need the first four fields; any extra columns are ignored.
            let mut fields = line.split('\t');
            let term = match (fields.next(), fields.next(), fields.next(), fields.next()) {
                (Some(_chrom), Some(_start), Some(_end), Some(term)) => term,
                _ => return Err(AnnotationError::MalformedBed { line: line_no + 1 }),
            };
            if term.is_empty() {
                return Err(AnnotationError::MalformedBed { line: line_no + 1 });
            }

            if index.contains_key(term) {
                continue; // already registered; bitset naturally dedups terms
            }

            let bit = names.len();
            if bit >= MAX_FEATURES {
                return Err(AnnotationError::TooManyFeatures {
                    path: source.to_owned(),
                    count: bit + 1, // this term is the (bit+1)-th distinct one
                    max: MAX_FEATURES,
                });
            }

            let boxed: Box<str> = Box::from(term);
            names.push(boxed.clone());
            index.insert(boxed, bit as u32);
        }

        Ok(Self {
            names: names.into_boxed_slice(),
            index,
        })
    }

    /// Number of distinct terms (equivalently, bits in use).
    #[inline(always)]
    pub fn num_features(&self) -> usize {
        self.names.len()
    }

    /// `true` when the BED defined no terms at all.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Bit **mask** (`1 << bit`) for a term, ready to OR into an annotation
    /// slot. `None` if the term is not in this registry.
    ///
    /// This is the hot-path accessor used by the annotation transform.
    #[inline]
    pub fn mask_of(&self, term: &str) -> Option<u32> {
        self.index.get(term).map(|&bit| 1u32 << bit)
    }

    /// Bit **index** (`0..num_features`) for a term, or `None` if absent.
    #[inline]
    pub fn bit_of(&self, term: &str) -> Option<u32> {
        self.index.get(term).copied()
    }

    /// Term name for a bit index, or `None` if the index is unused.
    #[inline]
    pub fn term(&self, bit: u32) -> Option<&str> {
        self.names.get(bit as usize).map(|s| &**s)
    }

    /// Iterate all registered terms in bit-index order (e.g. for a header).
    pub fn terms(&self) -> impl Iterator<Item = &str> + '_ {
        self.names.iter().map(|s| &**s)
    }

    /// Expand an accumulated annotation `u32` into the names of its set bits,
    /// in ascending bit order — the sink joins these with commas.
    ///
    /// Bits with no registered term are ignored rather than panicking. In this
    /// pipeline the slot is only ever filled with masks from *this* registry,
    /// so stray high bits should not occur; ignoring them is purely defensive,
    /// since a panic in a sink worker would silently truncate the report.
    pub fn decode(&self, bits: u32) -> impl Iterator<Item = &str> + '_ {
        let mut bits = bits & self.valid_mask();
        std::iter::from_fn(move || {
            if bits == 0 {
                return None;
            }
            let bit = bits.trailing_zeros();
            bits &= bits - 1; // clear lowest set bit
            Some(&*self.names[bit as usize])
        })
    }

    /// Mask covering exactly the bits that have a registered term.
    #[inline]
    fn valid_mask(&self) -> u32 {
        let n = self.names.len() as u32;
        // n is in 0..=32; `1 << 32` would overflow a u32, so special-case it.
        if n >= 32 {
            u32::MAX
        } else {
            (1u32 << n) - 1
        }
    }
}

/// Open a BED path as a line reader, transparently decompressing bgzip/gzip.
///
/// Detection is by the two-byte gzip magic (`1f 8b`) peeked without consuming,
/// so a mislabeled or extension-less file still reads correctly.
fn open_bed(path: &Path) -> Result<Box<dyn BufRead>, AnnotationError> {
    let file = File::open(path)
        .map_err(|e| AnnotationError::IoError(format!("{}: {e}", path.display())))?;
    let mut reader = BufReader::new(file);

    let gzipped = {
        // `fill_buf` peeks; it does not advance the cursor, so the magic bytes
        // are still delivered to whichever reader we hand `reader` to below.
        let head = reader
            .fill_buf()
            .map_err(|e| AnnotationError::IoError(e.to_string()))?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };

    if gzipped {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn build(lines: &[&str]) -> Result<FeatureRegistry, AnnotationError> {
        let bytes = lines.join("\n").into_bytes();
        FeatureRegistry::from_reader(Cursor::new(bytes), "<test>")
    }

    #[test]
    fn dedup_and_first_seen_bit_order() {
        let reg = build(&[
            "chr1\t10\t20\tCTCF",
            "chr1\t30\t40\tpromoter",
            "chr1\t50\t60\tCTCF", // duplicate -> no new bit
            "chr1\t70\t80\tenhancer",
        ])
        .unwrap();

        assert_eq!(reg.num_features(), 3);
        assert_eq!(reg.mask_of("CTCF"), Some(0b001));
        assert_eq!(reg.mask_of("promoter"), Some(0b010));
        assert_eq!(reg.mask_of("enhancer"), Some(0b100));
        assert_eq!(reg.bit_of("enhancer"), Some(2));
        assert_eq!(reg.mask_of("missing"), None);
        assert_eq!(reg.term(1), Some("promoter"));
        assert_eq!(
            reg.terms().collect::<Vec<_>>(),
            ["CTCF", "promoter", "enhancer"]
        );
    }

    #[test]
    fn decode_round_trip_low_to_high() {
        let reg = build(&["c\t0\t1\tCTCF", "c\t0\t1\tpromoter", "c\t0\t1\tenhancer"]).unwrap();

        let bits = reg.mask_of("enhancer").unwrap() | reg.mask_of("CTCF").unwrap();
        assert_eq!(bits, 0b101);
        assert_eq!(reg.decode(bits).collect::<Vec<_>>(), ["CTCF", "enhancer"]);
        assert_eq!(reg.decode(0).collect::<Vec<_>>(), Vec::<&str>::new());
        // stray high bit with no registered term is ignored, not a panic
        assert_eq!(
            reg.decode(bits | (1 << 20)).collect::<Vec<_>>(),
            ["CTCF", "enhancer"]
        );
    }

    #[test]
    fn exactly_max_features_ok() {
        let lines: Vec<String> = (0..MAX_FEATURES)
            .map(|k| format!("c\t0\t1\tt{k}"))
            .collect();
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let reg = build(&refs).unwrap();
        assert_eq!(reg.num_features(), MAX_FEATURES);
        // all 32 bits valid -> a full-slot decode lists every term
        assert_eq!(reg.decode(u32::MAX).count(), MAX_FEATURES);
    }

    #[test]
    fn one_past_max_features_fails() {
        let lines: Vec<String> = (0..=MAX_FEATURES)
            .map(|k| format!("c\t0\t1\tt{k}"))
            .collect();
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        match build(&refs) {
            Err(AnnotationError::TooManyFeatures { count, max, .. }) => {
                assert_eq!(max, MAX_FEATURES);
                assert_eq!(count, MAX_FEATURES + 1);
            }
            other => panic!("expected TooManyFeatures, got {other:?}"),
        }
    }

    #[test]
    fn skips_comments_blanks_and_strips_crlf() {
        let reg = build(&[
            "# tabix meta line",
            "",
            "chr2\t1\t2\tintron\r", // CRLF: trailing '\r' must not join the term
        ])
        .unwrap();
        assert_eq!(reg.num_features(), 1);
        assert_eq!(reg.bit_of("intron"), Some(0));
        assert_eq!(reg.bit_of("intron\r"), None);
    }

    #[test]
    fn malformed_lines_rejected() {
        assert!(matches!(
            build(&["chrX\t1\t2"]), // only 3 fields
            Err(AnnotationError::MalformedBed { line: 1 })
        ));
        assert!(matches!(
            build(&["chrX\t1\t2\t"]), // empty term
            Err(AnnotationError::MalformedBed { line: 1 })
        ));
    }

    #[test]
    fn reads_gzip_compressed_bed() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let payload = b"chrM\t0\t1\tCpG_island\nchrM\t5\t6\tCpG_island\n";
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(payload).unwrap();
        let gz = enc.finish().unwrap();
        assert_eq!(&gz[..2], &[0x1f, 0x8b]); // gzip magic

        let path = std::env::temp_dir().join(format!("crisprme_reg_{}.bed.gz", std::process::id()));
        std::fs::write(&path, &gz).unwrap();

        let reg = FeatureRegistry::from_bed(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(reg.num_features(), 1);
        assert_eq!(reg.bit_of("CpG_island"), Some(0));
    }
}
