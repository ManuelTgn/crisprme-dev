use std::{
    ffi::os_str::Display,
    fmt::{write, Debug},
};

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarxOp {
    Match = 0b00,     // `=` (sequence match)
    Mismatch = 0b01,  // `X` (sequence mismatch)
    Deletion = 0b10,  // `D` (deletion from reference)
    Insertion = 0b11, // `I` (insertion to reference)
}

/// All common operations supported by a cigarx
pub trait Cigarx {
    /// Number of operations that can be stored
    const CAPACITY: usize;

    /// Number of operations present
    fn len(&self) -> usize;

    fn push(&mut self, op: CigarxOp);
    fn pop(&mut self) -> Option<CigarxOp>;

    /// Returns the operations from left to right
    fn iter(&self) -> impl Iterator<Item = CigarxOp>;
}

/// Cigarx implemented as u64 with sentinel bit,
/// it supports up to 31 cigarx operations
#[repr(transparent)]
#[derive(Clone, Copy)]
pub struct Cigarx64(u64);

unsafe impl bytemuck::Pod for Cigarx64 {}
unsafe impl bytemuck::Zeroable for Cigarx64 {
    fn zeroed() -> Self {
        // This in an invalid cigarx, but a valid storage
        Self(0)
    }
}

/// A valid cigarx starts with the sentinel at the first bit
impl Default for Cigarx64 {
    fn default() -> Self {
        Self(1)
    }
}

impl Cigarx for Cigarx64 {
    // One bit is reserved for the sentinel
    const CAPACITY: usize = (64 - 1) / 2;

    fn len(&self) -> usize {
        if self.0 == 0 {
            return 0;
        }
        // zeroed state for bytemuck compatibility
        else {
            let available_bits = u64::BITS - self.0.leading_zeros() - 1;
            (available_bits / 2) as usize
        }
    }

    fn push(&mut self, op: CigarxOp) {
        assert!(self.len() < Self::CAPACITY, "Cigarx64 overflow");
        if self.0 == 0 {
            self.0 = 1
        }; // default state for bytemuck compatibility
        self.0 = (self.0 << 2) | (op as u64);
    }

    fn pop(&mut self) -> Option<CigarxOp> {
        if self.len() == 0 {
            None
        } else {
            let bits = (self.0 & 0b11) as u8;
            Some(match bits {
                0 => CigarxOp::Match,
                1 => CigarxOp::Mismatch,
                2 => CigarxOp::Deletion,
                _ => CigarxOp::Insertion,
            })
        }
    }

    fn iter(&self) -> impl Iterator<Item = CigarxOp> {
        Cigarx64Iter {
            remaining: self.len(),
            storage: self.0,
        }
    }
}

pub struct Cigarx64Iter {
    remaining: usize,
    storage: u64,
}

impl Iterator for Cigarx64Iter {
    type Item = CigarxOp;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            None
        } else {
            let shift = (self.remaining - 1) * 2;
            let bits = ((self.storage >> shift) & 0b11) as u8;
            self.remaining -= 1;
            Some(match bits {
                0 => CigarxOp::Match,
                1 => CigarxOp::Mismatch,
                2 => CigarxOp::Deletion,
                _ => CigarxOp::Insertion,
            })
        }
    }
}

impl Debug for Cigarx64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for elem in self.iter() {
            write!(
                f,
                "{}",
                match elem {
                    CigarxOp::Match => '=',
                    CigarxOp::Mismatch => 'X',
                    CigarxOp::Deletion => 'D',
                    CigarxOp::Insertion => 'I',
                }
            )?
        }
        Ok(())
    }
}

/// Cigarx implemented as u128 with sentinel bit,
/// it supports up to 63 cigarx operations
#[repr(transparent)]
#[derive(Clone, Copy)]
pub struct Cigarx128(u128);

unsafe impl bytemuck::Pod for Cigarx128 {}
unsafe impl bytemuck::Zeroable for Cigarx128 {
    fn zeroed() -> Self {
        // This in an invalid cigarx, but a valid storage
        Self(0)
    }
}

/// A valid cigarx starts with the sentinel at the first bit
impl Default for Cigarx128 {
    fn default() -> Self {
        Self(1)
    }
}

impl Cigarx for Cigarx128 {
    // One bit is reserved for the sentinel
    const CAPACITY: usize = (128 - 1) / 2;

    fn len(&self) -> usize {
        if self.0 == 0 {
            return 0;
        }
        // zeroed state for bytemuck compatibility
        else {
            let available_bits = u128::BITS - self.0.leading_zeros() - 1;
            (available_bits / 2) as usize
        }
    }

    fn push(&mut self, op: CigarxOp) {
        assert!(self.len() < Self::CAPACITY, "Cigarx128 overflow");
        if self.0 == 0 {
            self.0 = 1
        }; // default state for bytemuck compatibility
        self.0 = (self.0 << 2) | (op as u128);
    }

    fn pop(&mut self) -> Option<CigarxOp> {
        if self.len() == 0 {
            None
        } else {
            let bits = (self.0 & 0b11) as u8;
            Some(match bits {
                0 => CigarxOp::Match,
                1 => CigarxOp::Mismatch,
                2 => CigarxOp::Deletion,
                _ => CigarxOp::Insertion,
            })
        }
    }

    fn iter(&self) -> impl Iterator<Item = CigarxOp> {
        Cigarx128Iter {
            remaining: self.len(),
            storage: self.0,
        }
    }
}

pub struct Cigarx128Iter {
    remaining: usize,
    storage: u128,
}

impl Iterator for Cigarx128Iter {
    type Item = CigarxOp;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            None
        } else {
            let shift = (self.remaining - 1) * 2;
            let bits = ((self.storage >> shift) & 0b11) as u8;
            self.remaining -= 1;
            Some(match bits {
                0 => CigarxOp::Match,
                1 => CigarxOp::Mismatch,
                2 => CigarxOp::Deletion,
                _ => CigarxOp::Insertion,
            })
        }
    }
}

impl Debug for Cigarx128 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for elem in self.iter() {
            write!(
                f,
                "{}",
                match elem {
                    CigarxOp::Match => '=',
                    CigarxOp::Mismatch => 'X',
                    CigarxOp::Deletion => 'D',
                    CigarxOp::Insertion => 'I',
                }
            )?
        }
        Ok(())
    }
}
