//! Compact CIGARX encoding.
//!
//! This module defines a compact encoding for alignment operation strings
//! (CIGAR-like, but explicitly distinguishing match and mismatch):
//!
//! - `=` : match
//! - `X` : mismatch
//! - `D` : deletion
//! - `I` : insertion
//!
//! The primary purpose is to store per-alignment edit operations in a small,
//! fixed-size integer (`u64` by default) for:
//! - GPU -> CPU transfer (fixed-size records are easy to DMA and serialize)
//! - minimal memory footprint in ring buffers
//! - fast edit-distance computation by iterating ops
//!
//! # Encoding strategy (current implementation)
//!
//! The current implementation uses a **fixed 2-bit encoding** per operation.
//!
//! - Each operation is represented by 2 bits (values 0..3).
//! - The encoded sequence is stored in an integer `T` (`u32`, `u64`, ...).
//! - `bits` stores how many bits are currently used (always a multiple of 2).
//!
//! ## Order convention
//!
//! `push(op)` appends one operation at the *end* of the sequence.
//! Internally, we shift left by 2 and OR the new op bits.
//!
//! Therefore, the *oldest* operation resides at the highest used bits,
//! and the *newest* operation resides in the lowest used bits.
//!
//! - `operations()` iterates from oldest -> newest (left to right).
//! - `pop()` removes from the end (newest) -> oldest.
//!
//! This is analogous to a stack where `push` appends and `pop` removes the last op.
//!
//! # Capacity
//!
//! For a `T` of `N` bits (e.g., `u64` -> 64 bits), maximum operations is `N / 2`.
//!
//! If you need more than `N/2` operations:
//! - switch to a larger `T` (e.g., `u128`), or
//! - store multiple words, or
//! - revive the Huffman encoding approach (variable-length symbols).
//!

use std::hash::Hash;
use std::mem::size_of;

use num_traits::{PrimInt, Unsigned, Zero};

/// CIGARX operations encoded in 2 bits each.
///
/// Encoding (2-bit):
/// - `00` = Match (`=`),
/// - `01` = Mismatch (`X`),
/// - `10` = Deletion (`D`),
/// - `11` = Insertion (`I`).
///
/// `#[repr(u8)]` makes the numeric values stable and FFI-friendly.
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarxOp {
    // `=` (sequence match)
    Match = 0b00,
    // `X` (sequence mismatch)
    Mismatch = 0b01,
    // `D` (deletion from reference)
    Deletion = 0b10,
    // `I` (insertion to reference)
    Insertion = 0b11,
}

impl CigarxOp {
    /// Convert an operation to its single-character CIGARX representation.
    #[inline(always)]
    pub fn to_utf8(self) -> char {
        match self {
            CigarxOp::Match => '=',
            CigarxOp::Mismatch => 'X',
            CigarxOp::Deletion => 'D',
            CigarxOp::Insertion => 'I',
        }
    }

    /// Fallible conversion from a single-character representation.
    ///
    /// Used in parsing paths that must not panic.
    #[inline(always)]
    pub fn try_from_utf8(c: char) -> Option<CigarxOp> {
        match c {
            '=' => Some(CigarxOp::Match),
            'X' => Some(CigarxOp::Mismatch),
            'I' => Some(CigarxOp::Insertion),
            'D' => Some(CigarxOp::Deletion),
            _ => None,
        }
    }

    /// Infallible conversion from a single-character representation.
    ///
    /// # Panics
    /// Panics if `c` is not one of `= X D I`.
    #[inline(always)]
    pub fn from_utf8(c: char) -> CigarxOp {
        Self::try_from_utf8(c).expect("invalid CIGARX character")
    }

    /// Get the 2-bit integer value representing this operation.
    #[inline(always)]
    pub fn to_integer(self) -> u8 {
        self as u8
    }

    /// Create an operation from its low 2 bits.
    ///
    /// This always succeeds for 2-bit values 0..3.
    #[inline(always)]
    pub fn from_bits(bits: u8) -> Option<Self> {
        match bits & 0b11 {
            0b00 => Some(CigarxOp::Match),
            0b01 => Some(CigarxOp::Mismatch),
            0b10 => Some(CigarxOp::Deletion),
            0b11 => Some(CigarxOp::Insertion),
            _ => None,
        }
    }
}

/// Compact CIGARX stored in an integer backbone.
///
/// `storage` holds the packed bits.
/// `bits` indicates how many bits are in use (always multiple of 2 in fixed mode).
///
/// # Type parameter `T`
/// - Must be an unsigned integer type implementing `PrimInt + Unsigned + Zero`.
/// - Typical choices: `u32` (max 16 ops), `u64` (max 32 ops), `u128` (max 64 ops).
///
/// # Semantics
/// - `push()` appends an operation to the end.
/// - `pop()` removes the last appended operation.
/// - `operations()` iterates from the first pushed to the last pushed.
#[repr(C)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Cigarx<T: Hash> {
    storage: T,
    bits: u8,
}

impl<T> Cigarx<T>
where
    T: PrimInt + Zero + Unsigned + Hash,
{
    /// Maximum number of operations storable in this `Cigarx<T>` using 2-bit encoding.
    #[inline(always)]
    pub fn capacity_ops() -> u8 {
        (size_of::<T>() as u8 * 8) / 2
    }

    /// Number of operations currently stored.
    #[inline(always)]
    pub fn len_ops(&self) -> u8 {
        self.bits / 2
    }

    /// Returns `true` if no operations are stored.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.bits == 0
    }

    /// Append one operation to the end of the packed sequence.
    ///
    /// Complexity: O(1).
    ///
    /// # Panics
    /// Panics if there is not enough capacity in `T` to store another op.
    pub fn push(&mut self, op: CigarxOp) {
        let cap_bits = size_of::<T>() as u8 * 8;
        // Need room for 2 more bits
        assert!(
            self.bits + 2 <= cap_bits,
            "Cigarx overflow: bits={}, capacity_bits={}",
            self.bits,
            cap_bits,
        );

        self.storage = (self.storage << 2) | T::from(op.to_integer()).unwrap();
        self.bits += 2;
    }

    /// Remove and return the last operation (LIFO).
    ///
    /// This returns operations in reverse order of `push()`.
    ///
    /// Complexity: O(1).
    #[inline(always)]
    pub fn pop(&mut self) -> Option<CigarxOp> {
        if self.bits == 0 {
            return None;
        }

        // The last pushed op lives in the lowest 2 bits
        let val: T = self.storage & T::from(0b11).unwrap();
        self.storage = self.storage >> 2;
        self.bits -= 2;

        CigarxOp::from_bits(val.to_u8().unwrap())
    }

    /// Iterate all operations from oldest to newest (push order).
    ///
    /// Complexity: O(k) where k is number of operations.
    #[inline]
    pub fn operations(&self) -> CigarxIter<T> {
        CigarxIter {
            remaining_bits: self.bits,
            storage: self.storage,
        }
    }
}

/// Iterator over packed operations (oldest -> newest).
///
/// This iterator does not allocate and is `Copy`-cheap.
pub struct CigarxIter<T: Hash> {
    remaining_bits: u8,
    storage: T,
}

impl<T> Iterator for CigarxIter<T>
where
    T: PrimInt + Unsigned + Zero + Hash,
{
    type Item = CigarxOp;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_bits == 0 {
            return None;
        }

        // Oldest op lives at the top of the used bits
        let shift = self.remaining_bits as usize - 2;
        let val = (self.storage >> shift) & T::from(0b11).unwrap();
        self.remaining_bits -= 2;

        CigarxOp::from_bits(val.to_u8().unwrap())
    }
}

// =============================================================================
// STD implementations

impl<T> Default for Cigarx<T>
where
    T: Zero + Hash,
{
    /// Create an empty Cigarx (no operations).
    #[inline]
    fn default() -> Self {
        Self {
            storage: T::zero(),
            bits: 0,
        }
    }
}

impl<T> From<&str> for Cigarx<T>
where
    T: PrimInt + Zero + Unsigned + Hash,
{
    /// Parse an ASCII CIGARX string (e.g. "==XDI") into a packed representation.
    ///
    /// # Panics
    /// Panics if the string contains invalid characters.
    fn from(value: &str) -> Self {
        let mut result = Self::default();
        for c in value.chars() {
            result.push(CigarxOp::from_utf8(c));
        }
        result
    }
}

impl<T> std::fmt::Display for Cigarx<T>
where
    T: PrimInt + Unsigned + Zero + Hash,
{
    // Pretty-print as a CIGARX string.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for op in self.operations() {
            write!(f, "{}", op.to_utf8())?;
        }
        Ok(())
    }
}

impl<T> std::fmt::Debug for Cigarx<T>
where
    T: PrimInt + Unsigned + Zero + Hash,
{
    // Debug print includes length and rendered operations
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let cap_ops = (size_of::<T>() * 8) / 2;
        write!(f, "Cigarx(len: {}/{}, ops: \"", self.bits / 2, cap_ops)?;
        for op in self.operations() {
            write!(f, "{}", op.to_utf8())?;
        }
        write!(f, "\")")
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn encode_then_decode() {
        const OPERATIONS: [CigarxOp; 6] = [
            CigarxOp::Match,
            CigarxOp::Mismatch,
            CigarxOp::Mismatch,
            CigarxOp::Deletion,
            CigarxOp::Insertion,
            CigarxOp::Deletion,
        ];

        let mut cigarx = Cigarx::<u64>::default();
        for op in &OPERATIONS {
            cigarx.push(*op);
        }

        let result: Vec<CigarxOp> = cigarx.operations().collect();
        assert_eq!(result, OPERATIONS);
    }

    #[test]
    fn push_then_pop_lifo() {
        let mut cigarx = Cigarx::<u64>::default();
        cigarx.push(CigarxOp::Match);
        cigarx.push(CigarxOp::Deletion);
        cigarx.push(CigarxOp::Insertion);

        assert_eq!(cigarx.pop(), Some(CigarxOp::Insertion));
        assert_eq!(cigarx.pop(), Some(CigarxOp::Deletion));
        assert_eq!(cigarx.pop(), Some(CigarxOp::Match));
        assert_eq!(cigarx.pop(), None);
    }
}
