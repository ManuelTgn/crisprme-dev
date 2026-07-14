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
}

/// Errors raised while building the contig id -> name table for the report.
///
/// Both variants map to a Python `ValueError` through
/// [`crate::python::pyerrors`], consistent with the other config-time errors.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ContigLabelsError {
    /// The name table is empty; a report would have no contigs to label.
    #[error("contig name table is empty")]
    Empty,

    /// A name is empty, or contains a byte that would corrupt the CSV row.
    #[error("contig name {name:?} (id {id}) contains illegal byte {byte} \
             (one of ',' '\"' '\\n' '\\r'), which would break the CSV")]
    InvalidName { id: u32, name: String, byte: u8 },
}
