use std::hash::Hash;

use num_traits::{PrimInt, Unsigned, Zero};

/// Huffman codes for a cigarx 
///
/// We expect the majority of characters to be '=', then 'X' and then 'D' or 'I'
const CIGARX_HUFFMAN_CODE: [u64; 4] = [ 0b0, 0b10, 0b110, 0b111 ];

/// Huffman codes' len for a cigarx
const CIGARX_HUFFMAN_CLEN: [u8; 4] = [ 1, 2, 3, 3 ];

/// Cigarx operations encoded in 2 bits each
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarxOp {
    Match     = 0b00, // = (sequence match)
    Mismatch  = 0b01, // X (sequence mismatch)  
    Deletion  = 0b10, // D (deletion from reference)
    Insertion = 0b11, // I (insertion to reference)
}

impl CigarxOp {

    /// Get single character representation
    pub fn to_utf8(self) -> char {
        match self {
            CigarxOp::Match      => '=',
            CigarxOp::Mismatch   => 'X',
            CigarxOp::Deletion   => 'D',
            CigarxOp::Insertion  => 'I',
        }
    }

    /// From single character representation
    pub fn from_utf8(c: char) -> CigarxOp {
        match c {
            '=' => CigarxOp::Match,
            'X' => CigarxOp::Mismatch,
            'D' => CigarxOp::Deletion,
            'I' => CigarxOp::Insertion,
            _ => unimplemented!()
        }
    }


    /// Get 2-bit representation
    #[inline(always)]
    pub fn to_integer(self) -> u8 {
        self as u8
    }

    /// Create CigarOp from 2-bit value
    pub fn from_bits(bits: u8) -> Option<Self> {
        match bits & 0b11 {
            0b00 => Some(CigarxOp::Match),
            0b01 => Some(CigarxOp::Mismatch),
            0b10 => Some(CigarxOp::Deletion),
            0b11 => Some(CigarxOp::Insertion),
            _ => None,
        }
    }

    /// Convert a cigarx symbol to a huffman index, code and len
    pub fn huffman(&self) -> (u64, u8) {
        let index = self.to_integer() as usize;
        (
            CIGARX_HUFFMAN_CODE[index],
            CIGARX_HUFFMAN_CLEN[index]
        )
    }
}

// Huffman encoded CIGARX for processing with T integer backbone
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Cigarx<T: Hash> {
    storage: T,
    bits: u8
}

/*
impl<T> Cigarx<T>
where
    T: PrimInt + Zero + Unsigned + Hash
{
    /// Push an operation to the storage bits
    pub fn push(&mut self, op: CigarxOp) {
        let (hcode, hbits) = op.huffman();
        debug_assert!(self.bits + hbits <= 64, "{:064b}", self.bits);
        self.storage = (self.storage << hbits as usize) | T::from(hcode).unwrap();
        self.bits += hbits;
    }

    /// Pop an operation from the storage bits
    pub fn pop(&mut self) -> Option<CigarxOp> {
        
        // Keep only the bits inside the current len
        let k = T::from(1u64.wrapping_shl(self.bits as u32) - 1).unwrap();
        let masked = self.storage & k;

        // Huffman codes converted to T
        let he = T::from(CIGARX_HUFFMAN_CODE[0]).unwrap();
        let hx = T::from(CIGARX_HUFFMAN_CODE[1]).unwrap();
        let hd = T::from(CIGARX_HUFFMAN_CODE[2]).unwrap();
        let hi = T::from(CIGARX_HUFFMAN_CODE[3]).unwrap();

        // No more bits
        if self.bits != 0 {
            let a = masked >> (self.bits - 1) as usize;
            
            if a == he {
                self.bits -= 1;
                return Some(CigarxOp::Match);
            }
        }

        if self.bits >= 2 {
            let b = masked >> (self.bits - 2) as usize;
            
            if b == hx {
                self.bits -= 2;
                return Some(CigarxOp::Mismatch);
            }
        }

        if self.bits >= 3 {
            let c = masked >> (self.bits - 3) as usize;

            if c == hd {
                self.bits -= 3;
                return Some(CigarxOp::Deletion);
            }

            if c == hi {
                self.bits -= 3;
                return Some(CigarxOp::Insertion);
            }
        }



        None
    }

    /// Returns all operations present
    pub fn operations(&self) -> impl Iterator<Item = CigarxOp> {
        CigarxIter::<T> { cigarx: *self }
    }
}
*/

impl<T> Cigarx<T>
where
    T: PrimInt + Zero + Unsigned + Hash
{
    /// Push an operation to the storage bits
    pub fn push(&mut self, op: CigarxOp) {
        assert!(self.bits < size_of::<T>() as u8 * 8, "{:064b}", self.bits);
        self.storage = (self.storage << 2) | T::from(op.to_integer()).unwrap();
        self.bits += 2;
    }

    /// Pop an operation from the storage bits
    pub fn pop(&mut self) -> Option<CigarxOp> {
        if self.bits == 0 { return None; }

        let val: T = self.storage & T::from(0b11).unwrap(); // last 2 bits (LSB)
        self.storage = self.storage >> 2;
        self.bits -= 2;

        CigarxOp::from_bits(val.to_u8().unwrap())
    }

    /// Returns all operations present
    pub fn operations(&self) -> impl Iterator<Item = CigarxOp> {
        CigarxIter::<T> {
            remaining_bits: self.bits,
            storage: self.storage
        }
    }

}

/// Iterator for a encoded CIGARX
pub struct CigarxIter<T: Hash> {
    remaining_bits: u8,
    storage: T,
}

impl<T> Iterator for CigarxIter<T> 
where
    T: PrimInt + Unsigned + Zero + Hash
{
    type Item = CigarxOp;
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_bits == 0 {
            return None;
        }

        let shift = self.remaining_bits as usize - 2;
        let val = (self.storage >> shift) & T::from(0b11).unwrap();
        self.remaining_bits -= 2;

        CigarxOp::from_bits(val.to_u8().unwrap())
    }
}

/*
/// High bits count to low bits count
impl From<Cigarx<u64>> for Cigarx<u32> {
    fn from(value: Cigarx<u64>) -> Self {
        assert!(value.bits <= 32);
        Self {
            storage: value.storage as u32,
            bits: value.bits
        }
    }
}
*/

// ===================================================================================================
// STD implementations

impl<T> Default for Cigarx<T>
where
    T: Zero + Hash
{
    fn default() -> Self {
        Self {
            storage: T::zero(),
            bits: 0
        }
    }
}

impl<T> From<&str> for Cigarx<T> 
where 
    T: PrimInt + Zero + Unsigned + Hash 
{
    fn from(value: &str) -> Self {
        let mut result = Self::default();
        for c in value.chars() {
            result.push(CigarxOp::from_utf8(c));
        }
        result
    }
}

// Pretty print as an alignment string
impl<T> std::fmt::Display for Cigarx<T> 
where
    T: PrimInt + Unsigned + Zero + Hash
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for op in self.operations() {
            write!(f, "{}", op.to_utf8())?;
        }
        Ok(())
    }
}

// Debug show as an aligment string
impl<T> std::fmt::Debug for Cigarx<T> 
where
    T: PrimInt + Unsigned + Zero + Hash
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Cigarx(len: {}/{}, storage: \"", self.bits / 2, std::mem::size_of::<T>() * 8 / 2)?;
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
}
