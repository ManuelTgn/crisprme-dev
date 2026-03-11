//! Alignment core types.
//!
//! This module defines:
//! - [`AlignmentState`]: a compact, bit-packed state used during alignment exploration
//!   (CPU side search / traversal of edit operations).
//! - [`AlignmentOp`]: the operation cursor used when expanding a state.
//! - [`Alignment`]: the final alignment record written by the aligner,
//!   represented primarily by an encoded CIGARX.
//! - `visualize`: a debugging helper to render an alignment against guide/target.
//!
//! # Design goals
//!
//! 1. **Throughput**: `AlignmentState` is a single `u32` to minimize memory traffic.
//! 2. **Predictable limits**: field widths cap max query/target length and max error counts.
//! 3. **FFI friendliness**: [`Alignment`] is `#[repr(C)]` and can be transferred as raw bytes.
//!
//! # Important invariants / limits
//!
//! AlignmentState bit layout (LSB -> MSB):
//! ```text
//! [oper:2][mism:4][qgap:4][tgap:4][qidx:5][tidx:5][reserved:8]
//! ```
//!
//! - `oper` is a cursor for enumerating transitions (`Match`, `Deletion`, `Insertion`, `Exhausted`).
//! - `mism`, `qgap`, `tgap` are capped at 15 (4 bits).
//! - `qidx`, `tidx` are capped at 31 (5 bits), effectively limiting sequences to length <= 32
//!   in this state representation.
//!

use super::cigarx::{Cigarx, CigarxOp};
use super::thresholds::Thresholds;

/// The next operation to try when expanding an [`AlignmentState`].
///
/// This enum is used as a small "cursor" rather than a semantic description of
/// an already-applied operation. The `travel()` method advances the internal
/// cursor (`oper` field) to enumerate possible next transitions.
///
/// The values are intentionally small (2-bit) to fit into [`AlignmentState`].
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentOp {
    // Try consuming one base from both guide and target
    Match = 0,
    // Try consuming one base from guide only (gap in target)
    Insertion = 1,
    // Try consuming one base from target only (gap in guide)
    Deletion = 2,
    // No more operations remain for this state (cursor exhausted)
    Exhausted = 3,
}

/// Packed state for alignment exploration.
///
/// The state encodes indices, error counts, and an operation cursor.
///
/// Bit layout (LSB -> MSB):
/// ```text
/// [oper:2][mism:4][qgap:4][tgap:4][qidx:5][tidx:5][reserved:8]
/// ```
///
/// Field semantics:
/// - `qidx`: current query index (0..31)
/// - `tidx`: current target index (0..31)
/// - `qgap`: number of inserted bases (gap in target) so far
/// - `tgap`: number of deleted bases (gap in query) so far
/// - `mism`: number of mismatches so far
/// - `oper`: operation cursor used by [`AlignmentState::travel`]
///
/// # Limitations
/// This state representation is valid only for sequences with max length 32 and
/// max mismatch/gap counts 15.
#[derive(Clone, Copy)]
pub struct AlignmentState(u32);

macro_rules! bitfield_get {
    ($name:ident, $offset:ident, $bits:ident) => {
        #[inline(always)]
        pub fn $name(self) -> u32 {
            (self.0 >> Self::$offset) & (1u32.wrapping_shl(Self::$bits) - 1)
        }
    };
}

macro_rules! bitfield_set {
    ($name:ident, $offset:ident, $bits:ident) => {
        #[inline(always)]
        pub fn $name(mut self, value: u32) -> Self {
            let mask = 1u32.wrapping_shl(Self::$bits) - 1;
            self.0 &= !(mask << Self::$offset);
            self.0 |= (value & mask) << Self::$offset;
            self
        }
    };
}

impl AlignmentState {
    // ---- Bit widths (8 bits are reserved at the top) ----
    const BITS_OPER: u32 = 2; // Maximum operations: 4
    const BITS_MISM: u32 = 4; // Maximum number of mismatches: 16
    const BITS_QGAP: u32 = 4; // Maximum number of query gaps: 16
    const BITS_TGAP: u32 = 4; // Maximum number of target gaps: 16
    const BITS_QIDX: u32 = 5; // Maximum len of query: 32
    const BITS_TIDX: u32 = 5; // Maximum len of target: 32

    // ---- Bit offsets ----
    const OFFSET_OPER: u32 = 0;
    const OFFSET_MISM: u32 = Self::OFFSET_OPER + Self::BITS_OPER;
    const OFFSET_QGAP: u32 = Self::OFFSET_MISM + Self::BITS_MISM;
    const OFFSET_TGAP: u32 = Self::OFFSET_QGAP + Self::BITS_QGAP;
    const OFFSET_QIDX: u32 = Self::OFFSET_TGAP + Self::BITS_TGAP;
    const OFFSET_TIDX: u32 = Self::OFFSET_QIDX + Self::BITS_QIDX;

    // ---- Getters ----
    bitfield_get!(oper, OFFSET_OPER, BITS_OPER);
    bitfield_get!(mism, OFFSET_MISM, BITS_MISM);
    bitfield_get!(qgap, OFFSET_QGAP, BITS_QGAP);
    bitfield_get!(tgap, OFFSET_TGAP, BITS_TGAP);
    bitfield_get!(qidx, OFFSET_QIDX, BITS_QIDX);
    bitfield_get!(tidx, OFFSET_TIDX, BITS_TIDX);

    // ---- Setters ----
    bitfield_set!(with_oper, OFFSET_OPER, BITS_OPER);
    bitfield_set!(with_mism, OFFSET_MISM, BITS_MISM);
    bitfield_set!(with_qgap, OFFSET_QGAP, BITS_QGAP);
    bitfield_set!(with_tgap, OFFSET_TGAP, BITS_TGAP);
    bitfield_set!(with_qidx, OFFSET_QIDX, BITS_QIDX);
    bitfield_set!(with_tidx, OFFSET_TIDX, BITS_TIDX);

    /// Create an initial state positioned at a given target index.
    ///
    /// Common usage: start exploring alignments anchored at a candidate
    /// starting location in the target sequence.
    ///
    /// The query index, error counters, and operation cursor are initialized to 0.
    #[inline]
    pub fn initial(tidx: u32) -> Self {
        // Defensive debug check: tidx must fit in the packed field
        debug_assert!(tidx < (1 << Self::BITS_TIDX));
        AlignmentState(0).with_tidx(tidx)
    }

    /// Returns `true` if the state is still within the given thresholds.
    ///
    /// This is a fast pruning predicate: any state exceeding mismatch or gap
    /// limits cannot lead to a valid alignment.
    #[inline(always)]
    pub fn below_thresholds(self, thresholds: &Thresholds) -> bool {
        self.mism() <= thresholds.mism
            && self.tgap() <= thresholds.tgap
            && self.qgap() <= thresholds.qgap
    }

    /// Returns `true` if the operation cursor has no remaining transitions.
    ///
    /// This does **not** mean the alignment is complete; it only means that for
    /// this `(qidx, tidx, errors)` state we have already enumerated all allowed
    /// "next operations".
    #[inline(always)]
    pub fn exhausted(self) -> bool {
        self.oper() == AlignmentOp::Exhausted as u32
    }

    /// Enumerate the next operation to try from this state.
    ///
    /// This method returns:
    /// - the updated state with advanced cursor (`oper += 1`)
    /// - the operation corresponding to the previous cursor value
    ///
    /// Operation order (as currently implemented):
    /// 1. Match
    /// 2. Deletion
    /// 3. Insertion
    /// 4. Exhausted
    ///
    /// The returned state always has its `oper` cursor advanced.
    #[inline(always)]
    pub fn travel(self) -> (Self, AlignmentOp) {
        let oper = self.oper();
        let result = self.with_oper(oper + 1);
        match oper {
            0 => (result, AlignmentOp::Match),
            1 => (result, AlignmentOp::Deletion),
            2 => (result, AlignmentOp::Insertion),
            _ => (result, AlignmentOp::Exhausted),
        }
    }

    /// Apply an insertion transition:
    /// - consumes one query base (`qidx += 1`)
    /// - increases query gaps (`qgap += 1`)
    /// - resets operation cursor (`oper = 0`)
    ///
    /// Interpretation: gap in target (query has extra base).
    #[inline(always)]
    pub fn op_insert(self) -> Self {
        self.with_qgap(self.qgap() + 1).with_qidx(self.qidx() + 1).with_oper(0)
    }

    /// Apply a deletion transition:
    /// - consumes one target base (`tidx += 1`)
    /// - increases target gaps (`tgap += 1`)
    /// - resets operation cursor (`oper = 0`)
    ///
    /// Interpretation: gap in query (target has extra base).
    #[inline(always)]
    pub fn op_delete(self) -> Self {
        self.with_tgap(self.tgap() + 1).with_tidx(self.tidx() + 1).with_oper(0)
    }

    /// Apply a match transition:
    /// - consumes one base from both query and target
    /// - resets operation cursor (`oper = 0`)
    #[inline(always)]
    pub fn op_match(self) -> Self {
        self.with_qidx(self.qidx() + 1).with_tidx(self.tidx() + 1).with_oper(0)
    }

    /// Apply a mismatch transition:
    /// - increases mismatch count
    /// - then applies a match-like advance (consume one base from both)
    ///
    /// Equivalent to: `mism += 1; qidx += 1; tidx += 1; oper = 0`.
    #[inline(always)]
    pub fn op_mismatch(self) -> Self {
        self.with_mism(self.mism() + 1).op_match().with_oper(0)
    }

    /// Returns `true` if the state violates thresholds.
    #[inline(always)]
    pub fn invalid(self, config: &Thresholds) -> bool {
        self.qgap() > config.qgap || self.tgap() > config.tgap || self.mism() > config.mism
    }

    /// Edit distance for this state (Levenshtein under this model).
    ///
    /// Defined as: `mism + qgap + tgap`.
    #[inline(always)]
    pub fn ed(self) -> u32 {
        self.qgap() + self.tgap() + self.mism()
    }
}

/// A complete alignment record.
///
/// This is the unit that alignment pipeline outputs (e.g., in an
/// `AlignmentRingBatch`), typically later copied GPU->CPU and serialized.
///
/// Layout is `#[repr(C)]` to keep field order stable for FFI and binary IO.
///
/// Fields:
/// - `cigarx`: compact operation encoding (supports mismatch/ins/del)
/// - `id`: target identifier (e.g., index into a target list or global ID)
/// - `offset`: alignment start offset relative to the target window
/// - `strand`: strand indicator (you likely use 0/1 or '+'/'-')
#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct Alignment {
    pub cigarx: Cigarx<u64>,
    pub id: u32,
    pub offset: u8,
    pub strand: u8,
}

impl Alignment {
    /// Construct an alignment record.
    ///
    /// `cigarx` may be provided as a concrete `Cigarx<u64>` or any type that
    /// converts into it.
    pub fn new<C: Into<Cigarx<u64>>>(id: u32, cigarx: C, offset: u8, strand: u8) -> Self {
        Self {
            id,
            cigarx: cigarx.into(),
            offset,
            strand
        }
    }

    /// Compute the edit distance implied by the CIGARX operations.
    ///
    /// Costs:
    /// - Match: 0
    /// - Mismatch / Insertion / Deletion: 1
    ///
    /// Note: this iterates all operations.
    pub fn ed(&self) -> u32 {
        self.cigarx
            .operations()
            .map(|op| match op {
                CigarxOp::Insertion | CigarxOp::Deletion | CigarxOp::Mismatch => 1,
                CigarxOp::Match => 0,
            })
            .sum()
    }
}

/// Build a human-readable alignment visualization.
///
/// This function is intended for debugging, unit tests, and sanity checks.
/// It returns a formatted multi-line string rather than printing directly.
/// Callers can `println!` if desired.
///
/// # Inputs
/// - `query`: ASCII bytes of the guide/query (e.g., b"ACGT...")
/// - `target`: ASCII bytes of the target sequence
/// - `cigar`: ASCII CIGAR operations (M/= /X/D/I)
/// - `start_pos`: alignment starting position in target
///
/// # Notes / limitations
/// - This is not optimized and should not be used in hot paths.
/// - Unknown CIGAR op bytes are rendered as '?' rather than panicking.
pub fn visualize(query: &[u8], target: &[u8], cigar: &[u8], start_pos: usize) -> String {
    let mut qline = String::new();
    let mut mline = String::new();
    let mut cline = String::new();
    let mut tline = String::new();

    let mut qidx: usize = 0;
    let mut tidx: usize = 0;

    // Add target prefix (unaligned)
    let prefix_len = start_pos.min(target.len());
    for i in 0..prefix_len {
        tline.push(target[i] as char);
        qline.push(' ');
        cline.push(' ');
        mline.push(' ');
        tidx += 1;
    }

    // Alignment
    for &op in cigar {
        cline.push(op as char);
        match op {
            b'M' | b'=' | b'X' => {
                let qc = query.get(qidx).copied().unwrap_or(b'?') as char;
                let tc = target.get(tidx).copied().unwrap_or(b'?') as char;

                qline.push(qc);
                tline.push(tc);
                mline.push(if qc == tc { '|' } else { ' ' });

                qidx += 1;
                tidx += 1;
            }

            b'D' => {
                let qc = query.get(qidx).copied().unwrap_or(b'?') as char;
                qline.push(qc);
                tline.push('-');
                mline.push(' ');
                qidx += 1;
            }
            b'I' => {
                let tc = target.get(tidx).copied().unwrap_or(b'?') as char;
                qline.push('-');
                tline.push(tc);
                mline.push(' ');
                tidx += 1;
            }
            _ => {
                // Unknwon operation
                qline.push('?');
                tline.push('?');
                mline.push(' ');
            }
        }
    }

    // Add target suffix (unaligned)
    while tidx < target.len() {
        tline.push(target[tidx] as char);
        qline.push(' ');
        mline.push(' ');
        tidx += 1;
    }

    // println!("target: {}", tline);
    // println!("cigarx: {}", cline);
    // println!(" guide: {}", qline);

    format!(
        "target: {tline}\n\
         cigarx: {cline}\n\
          guide: {qline}\n\
                {mline}\n"
    )
}

// =============================================================================
// STD implementations

impl std::fmt::Debug for AlignmentState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "AlignmentState(oper: {}, qidx: {}, tidx: {}, tgap: {}, qgap: {}, mism: {})",
            self.oper(),
            self.qidx(),
            self.tidx(),
            self.tgap(),
            self.qgap(),
            self.mism()
        )
    }
}

impl std::fmt::Debug for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Alignment(id: {}, offset: {}, strand: {}, cigarx: {:?})",
            self.id, self.offset, self.strand as char, self.cigarx
        )
    }
}

impl std::fmt::Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.cigarx)
    }
}
