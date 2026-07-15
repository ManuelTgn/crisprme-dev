use std::fmt;

use bytemuck::{Pod, Zeroable};

/// Genomic strand of an occurrence.
///
/// # Bit encoding
/// The encoding is fixed by `sequence::scanner`, which pushes `1` for a
/// forward-strand PAM hit and `0` for a reverse-strand hit. It is **not**
/// the intuitive `0 = forward`. Do not open-code the comparison; go through
/// [`Strand::from_bit`].
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Reverse = 0,
    Forward = 1,
}

impl Strand {
    #[inline(always)]
    pub const fn from_bit(bit: u8) -> Self {
        if bit & 1 == 1 {
            Self::Forward
        } else {
            Self::Reverse
        }
    }

    #[inline(always)]
    pub const fn as_bit(self) -> u8 {
        self as u8
    }

    /// Report representation: `"+"` forward, `"-"` reverse (BED / GFF convention).
    #[inline(always)]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "+",
            Self::Reverse => "-",
        }
    }

    /// Was a window on this strand found by scanning the **reverse-complemented**
    /// chunk?
    ///
    /// The miner requires the PAM at the window's right edge, so Python feeds the
    /// scanner the chunk orientation that puts it there. Which physical chunk that
    /// is depends on *both* the reported strand and where the PAM sits relative to
    /// the protospacer:
    ///
    /// | PAM placement          | `+` strand   | `−` strand   |
    /// |------------------------|--------------|--------------|
    /// | downstream (Cas9, NGG) | forward chunk| **RC chunk** |
    /// | upstream (Cas12, TTTV) | **RC chunk** | forward chunk|
    ///
    /// which collapses to `is_forward == upstream`.
    #[inline(always)]
    pub const fn scanned_on_revcomp(self, upstream: bool) -> bool {
        matches!(self, Self::Forward) == upstream
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

#[repr(transparent)]
#[derive(Debug, Clone, Copy)]
pub struct Occurence(pub u64);

impl Occurence {
    #[inline(always)]
    pub fn new(contig: u16, pam: u16, position: u32, strand: Strand) -> Self {
        Self(
            (((contig & 0x7FFF) as u64) << 49)
                | ((pam as u64) << 33)
                | ((position as u64) << 1)
                | (strand.as_bit() as u64),
        )
    }

    #[inline(always)]
    pub fn contig(&self) -> u16 {
        ((self.0 >> 49) & 0x7FFF) as u16
    }

    #[inline(always)]
    pub fn pam(&self) -> u16 {
        ((self.0 >> 33) & 0xFFFF) as u16
    }

    /// Contig-local genomic position.
    ///
    /// `position` occupies bits 1..=32, so the mask is 32 bits wide.
    #[inline(always)]
    pub fn position(&self) -> u32 {
        ((self.0 >> 1) & 0x_FFFF_FFFF) as u32
    }

    #[inline(always)]
    pub const fn strand(&self) -> Strand {
        Strand::from_bit((self.0 & 1) as u8)
    }
}

unsafe impl Zeroable for Occurence {}
unsafe impl Pod for Occurence {}
