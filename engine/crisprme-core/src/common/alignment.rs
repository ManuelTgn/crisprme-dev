use bytemuck::{Pod, Zeroable};

use super::{
    cigarx::{Cigarx, CigarxOp},
    thresholds::Thresholds,
};

/// The current alignment operation on a state
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentOp {
    Match = 0,
    Insertion = 1,
    Deletion = 2,
    Exhausted = 3,
}

/// Packed state for alignment exploration
/// [tidx:5][qidx:5][tgap:4][qgap:4][mism:4][oper:2][reserved:8]
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
    // Bits size of all the sections (8 bits are reserved)
    const BITS_OPER: u32 = 2; // Maximum operations: 4
    const BITS_MISM: u32 = 4; // Maximum number of mismatches: 16
    const BITS_QGAP: u32 = 4; // Maximum number of query gaps: 16
    const BITS_TGAP: u32 = 4; // Maximum number of target gaps: 16
    const BITS_QIDX: u32 = 5; // Maximum len of query: 32
    const BITS_TIDX: u32 = 5; // Maximum len of target: 32

    // Bits offset of all the sections
    const OFFSET_OPER: u32 = 0;
    const OFFSET_MISM: u32 = Self::OFFSET_OPER + Self::BITS_OPER;
    const OFFSET_QGAP: u32 = Self::OFFSET_MISM + Self::BITS_MISM;
    const OFFSET_TGAP: u32 = Self::OFFSET_QGAP + Self::BITS_QGAP;
    const OFFSET_QIDX: u32 = Self::OFFSET_TGAP + Self::BITS_TGAP;
    const OFFSET_TIDX: u32 = Self::OFFSET_QIDX + Self::BITS_QIDX;

    // Get the value of a section
    bitfield_get!(oper, OFFSET_OPER, BITS_OPER);
    bitfield_get!(mism, OFFSET_MISM, BITS_MISM);
    bitfield_get!(qgap, OFFSET_QGAP, BITS_QGAP);
    bitfield_get!(tgap, OFFSET_TGAP, BITS_TGAP);
    bitfield_get!(qidx, OFFSET_QIDX, BITS_QIDX);
    bitfield_get!(tidx, OFFSET_TIDX, BITS_TIDX);

    // Set the value of a section
    bitfield_set!(with_oper, OFFSET_OPER, BITS_OPER);
    bitfield_set!(with_mism, OFFSET_MISM, BITS_MISM);
    bitfield_set!(with_qgap, OFFSET_QGAP, BITS_QGAP);
    bitfield_set!(with_tgap, OFFSET_TGAP, BITS_TGAP);
    bitfield_set!(with_qidx, OFFSET_QIDX, BITS_QIDX);
    bitfield_set!(with_tidx, OFFSET_TIDX, BITS_TIDX);

    /// Create an initial state from a starting target sequence index
    pub fn initial(tidx: u32) -> Self {
        AlignmentState(0).with_tidx(tidx)
    }

    /// Check state against a set of thresholds
    #[inline(always)]
    pub fn below_thresholds(self, thresholds: &Thresholds) -> bool {
        self.mism() <= thresholds.mism
            && self.tgap() <= thresholds.tgap
            && self.qgap() <= thresholds.qgap
    }

    /// We explored all possible operations
    #[inline(always)]
    pub fn exhausted(self) -> bool {
        self.oper() == AlignmentOp::Exhausted as u32
    }

    /// Get next travel operation from state, until it is exausted
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

    /// Create a new state with a insertion
    #[inline(always)]
    pub fn op_insert(self) -> Self {
        self.with_qgap(self.qgap() + 1).with_qidx(self.qidx() + 1).with_oper(0)
    }

    /// Create a new state with a deletion
    #[inline(always)]
    pub fn op_delete(self) -> Self {
        self.with_tgap(self.tgap() + 1).with_tidx(self.tidx() + 1).with_oper(0)
    }

    /// Create a new state with a match
    #[inline(always)]
    pub fn op_match(self) -> Self {
        self.with_qidx(self.qidx() + 1).with_tidx(self.tidx() + 1).with_oper(0)
    }

    /// Create a new state with a mismatch
    #[inline(always)]
    pub fn op_mismatch(self) -> Self {
        self.with_mism(self.mism() + 1).op_match().with_oper(0)
    }

    #[inline(always)]
    pub fn invalid(self, config: &Thresholds) -> bool {
        self.qgap() > config.qgap || self.tgap() > config.tgap || self.mism() > config.mism
    }

    #[inline(always)]
    pub fn ed(self) -> u32 {
        self.qgap() + self.tgap() + self.mism()
    }
}

/// An alignment as a encoded CIGARX and an inner offset
#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct Alignment {
    pub cigarx: Cigarx<u64>,
    pub id: u32,
    pub offset: u8,
    pub strand: u8,
}

impl Alignment {
    /// Create alignment from partial
    pub fn new<C: Into<Cigarx<u64>>>(id: u32, cigarx: C, offset: u8, strand: u8) -> Self {
        Self {
            id,
            cigarx: cigarx.into(),
            offset,
            strand
        }
    }

    /// Returns the edit-distance of this alignment
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

/// Visualize the alignment relative to a target sequence and a guide
pub fn visualize(query: &[u8], target: &[u8], cigar: &[u8], start_pos: usize) {
        let mut qline = String::new();
        let mut mline = String::new();
        let mut cline = String::new();
        let mut tline = String::new();

        let mut qidx: usize = 0;
        let mut tidx: usize = 0;

        // Step 1: Add target prefix (unaligned)
        for i in 0..start_pos {
            tline.push(target[i] as char);
            qline.push(' ');
            cline.push(' ');
            mline.push(' ');
            tidx += 1;
        }

        // Step 2: Alignment
        for op in cigar {
            cline.push(*op as char);
            match op {
                b'M' | b'=' | b'X' => {
                    let qc = if qidx < query.len() { query[qidx] as char } else { '?' };
                    let tc = if tidx < target.len() { target[tidx] as char } else { '?' };
                    qline.push(qc);
                    tline.push(tc);
                    mline.push(if qc == tc { '|' } else { ' ' });
                    qidx += 1;
                    tidx += 1;
                }
                b'D' => {
                    let qc = if qidx < query.len() { query[qidx] as char } else { '?' };
                    qline.push(qc);
                    tline.push('-');
                    mline.push(' ');
                    qidx += 1;
                }
                b'I' => {
                    let tc = if tidx < target.len() { target[tidx] as char } else { '?' };
                    qline.push('-');
                    tline.push(tc);
                    mline.push(' ');
                    tidx += 1;
                }
                _ => unimplemented!(),
            }
        }

        // Step 3: Add target suffix (unaligned)
        while tidx < target.len() {
            tline.push(target[tidx] as char);
            qline.push(' ');
            mline.push(' ');
            tidx += 1;
        }

        println!("target: {}", tline);
        println!("cigarx: {}", cline);
        println!(" guide: {}", qline);
    }



// ===================================================================================================
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
