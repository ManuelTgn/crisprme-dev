use thiserror::Error;

use crate::error;

#[derive(Debug, Error)]
pub enum AnnotationError {
    #[error("Failed to read BED file: {0}")]
    IoError(String),

    #[error("BED file is malformed at line {line}")]
    MalformedBed { line: usize },

    #[error("Invalid feature ID {0}")]
    InvalidFeatureId(usize),

    #[error("Annotation input is empty")]
    EmptyInput,
}

/// Errors raised while parsing a PAM string or working with its
/// finite set of concrete variants.
///
/// All variants map to a Python exception through
/// [`crate::python::pyerrors`], so a Python caller always receives a
/// descriptive, typed error rather than an opaque Rust panic.
#[derive(Debug, Error)]
pub enum PamError {
    /// The PAM string contained a byte that is not a valid IUPAC code.
    #[error("invalid PAM character at position {position} (ASCII byte {byte})")]
    InvalidCharacter { position: usize, byte: u8 },

    /// The PAM is so degenerate that its concrete-variant count exceeds
    /// the range addressable by the `u16` variant index used downstream.
    ///
    /// This is effectively unreachable for real PAMs (a length-8 all-`N`
    /// PAM already reaches the ceiling) and exists purely to make the
    /// `u16` index representation provably safe.
    #[error("PAM defines {count} concrete variants (plen={plen}); \
             exceeds the {max} addressable by a u16 index")]
    TooManyVariants { count: u64, plen: usize, max: u32 },

    /// A variant index handed to the decoder is out of range.
    #[error("PAM variant index {index} out of range (valid: 0..{count})")]
    IndexOutOfRange { index: u16, count: u32 },

    #[error("PAM has {count} unconstrained (N) positions; the maximum supported is {max}")]
    TooManyWildcards { count: usize, max: usize },
}

/// Errors raised while constructing a candidate target or working with it.
///
/// All variants map to a Python exception through
/// [`crate::python::pyerrors`], so a Python caller always receives a
/// descriptive, typed error rather than an opaque Rust panic.
#[derive(Debug, Error)]
pub enum TargetError {
    #[error("PAM length {plen} is invalid for window size {size} (must be 1..={size})")]
    PamOutOfRange { plen: usize, size: usize },
    #[error("window size {size} exceeds the maximum supported length {max}")]
    WindowTooLong { size: usize, max: usize },
}

