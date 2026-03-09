use pyo3::{pyclass, pymethods};

use crate::sequence::iupac::Iupac;

/// A CRISPR guide sequence encoded as IUPAC 4-bit masks.
///
/// `Guide` is a thin, owned wrapper around a contiguous `Vec<Iupac>`.
/// The representation is intended to be:
/// - **compact** (1 byte per position, assuming `Iupac` is `#[repr(transparent)]` over `u8`)
/// - **FFI-friendly** for passing raw pointers to CUDA/C++
/// - **cheap to iterate** in the hot alignment/mining path
///
/// # Encoding
/// Each `Iupac` stores a 4-bit mask representing allowed bases at a position.
/// This supports degenerate guides (e.g., containing `N`, `R`, etc.).
///
/// # FFI / CUDA notes
/// If you pass `Guide` bytes to CUDA kernels:
/// - Ensure `Iupac` is `#[repr(transparent)] pub struct Iupac(u8)`
/// - Ensure the CUDA side interprets each byte as the same 4-bit mask encoding
/// - Use [`Guide::as_ptr`] and [`Guide::len`] for pointer+length.
///
/// # Invariants
/// - `self.0.len()` is the guide length in bases.
/// - Each element is an IUPAC mask; invalid input should be handled at construction time.
///
/// # Construction policy
/// This module provides both strict and lossy constructors:
/// - **strict**: reject invalid ASCII
/// - **lossy**: map invalid ASCII to `N` (wildcard), useful for robustness in pipelines
#[repr(transparent)]
#[derive(Clone, PartialEq, Eq, Default)]
#[pyclass]
pub struct Guide(Vec<Iupac>);

impl Guide {
    // Construct a guide from already-encoded IUPAC codes
    #[inline]
    pub fn from_iupac(vec: Vec<Iupac>) -> Self {
        Self(vec)
    }

    /// Strictly construct a guide from ASCII bytes.
    ///
    /// Returns `None` if any byte is not a valid IUPAC code.
    ///
    /// Prefer this in validation-heavy paths (CLI parsing, tests, etc.).
    #[inline]
    pub fn try_from_ascii_bytes(bytes: &[u8]) -> Option<Self> {
        let mut v = Vec::with_capacity(bytes.len());
        for &b in bytes {
            v.push(Iupac::try_from_ascii(b)?);
        }
        Some(Self(v))
    }

    /// Lossy construction from ASCII bytes: invalid chars are mapped to `N`.
    ///
    /// Prefer this for high-throughput pipelines where you do not want to fail hard.
    #[inline]
    pub fn from_ascii_bytes_lossy(bytes: &[u8]) -> Self {
        let mut v = Vec::with_capacity(bytes.len());
        for &b in bytes {
            v.push(Iupac::from_ascii_lossy(b));
        }
        Self(v)
    }

    /// Reverse-complement the guide.
    ///
    /// This creates a new `Guide` where:
    /// - the order is reversed, and
    /// - each position is complemented (A<->T, C<->G at bitmask level),
    /// including correct behavior for degenerate IUPAC codes.
    ///
    /// Complexity: O(n) time and O(n) additional memory.
    #[inline]
    pub fn reverse_complement(&self) -> Self {
        let n = self.0.len();
        let mut out = Vec::with_capacity(n);

        // Iterate from end to start and complement each code
        // This is typically slightly cheaper than iter().cloned().map(...).rev().collect()
        for &c in self.0.iter().rev() {
            out.push(c.complement());
        }
        Self(out)
    } 

    /// Render the guide as a UTF-8 string (IUPAC characters).
    ///
    /// Intended for logging/debugging and user output.
    /// For hot paths, prefer `as_ptr()/len()` instead.
    #[inline]
    pub fn to_string_iupac(&self) -> String {
        // Allocate exactly as many bytes as the guide length (ASCII output)
        let mut s = String::with_capacity(self.0.len());
        for &c in &self.0 {
            s.push(c.to_utf8());
        }
        s
    }

    /// Length of the guide in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if the guide is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Immutable view of the encoded guide.
    #[inline]
    pub fn as_slice(&self) -> &[Iupac] {
        &self.0
    }

    /// Pointer to the underlying contiguous `Iupac` buffer.
    ///
    /// Useful for passing into FFI (e.g., `*const u8` after casting).
    #[inline]
    pub fn as_ptr(&self) -> *const Iupac {
        self.0.as_ptr()
    }

    /// Iterate over the guide.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, Iupac> {
        self.0.iter()
    }
}

#[pymethods]
impl Guide {
    #[new]
    pub fn new(s: &str) -> Self {
        Self::from_ascii_bytes_lossy(s.as_bytes())
    }
}

/// Construct from a `&str`.
///
/// IMPORTANT: a Rust `&str` is UTF-8; using `.chars()` can accept non-ASCII and
/// can produce multi-byte characters. Since guides are expected to be ASCII IUPAC,
/// we parse bytes.
///
/// This implementation is **lossy** to avoid panics in pipelines.
/// If you want strict behavior, use `Guide::try_from_ascii_bytes(value.as_bytes())`.
impl From<&str> for Guide {
    #[inline]
    fn from(value: &str) -> Self {
        Self::from_ascii_bytes_lossy(value.as_bytes())
    }
}

/// Indexing inside a guide.
/// Panics if `idx` is out of bounds (standard Rust indexing semantics).
impl std::ops::Index<usize> for Guide {
    type Output = Iupac;
    #[inline]
    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

impl std::fmt::Debug for Guide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Same representation as Display, but keeps Debug trait available
        write!(f, "{}", self)
    }    
}

impl std::fmt::Display for Guide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for &c in &self.0 {
            write!(f, "{}", c.to_utf8())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn reverse_complement_normal() {
        const GUIDE_SRC: &str = "ATTGAGATAGTGGNGG";
        const GUIDE_REV: &str = "CCNCCACTATCTCAAT";

        let guide: Guide = GUIDE_SRC.into();
        let rev = guide.reverse_complement();

        let out = rev.to_string_iupac();
        assert_eq!(out, GUIDE_REV);
    }

    #[test]
    fn reverse_complement_iupac() {
        const GUIDE_SRC: &str = "RSH";
        const GUIDE_REV: &str = "DSY";

        let guide: Guide = GUIDE_SRC.into();
        let out = guide.reverse_complement().to_string_iupac();
        assert_eq!(out, GUIDE_REV, "output = {out}, correct = {GUIDE_REV}");
    }

    #[test]
    fn strict_constructor_rejects_invalid() {
        // 'Z' is not a valid IUPAC code
        let g = Guide::try_from_ascii_bytes(b"ACZ");
        assert!(g.is_none());
    }

    #[test]
    fn lossy_constructor_maps_invalid_to_n() {
        let g = Guide::from_ascii_bytes_lossy(b"ACZ");
        let s = g.to_string_iupac();
        // Z -> N in lossy mode
        assert_eq!(s, "ACN");
    }
}
