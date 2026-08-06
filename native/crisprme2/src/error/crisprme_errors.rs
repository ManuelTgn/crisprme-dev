use thiserror::Error;

use crate::error;

/// Fixed columns in every report line: `chromosome` .. `aggregate_score`.
pub const FIXED_COLS: usize = 14;
/// Upper bound on appended annotation columns (one bit each in the `u32` slot).
pub const MAX_ANNOTATIONS: usize = 10;

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

    #[error(
        "annotation BED {path} defines more than {max} distinct terms in \
         column 4 (hit term #{count}); each term needs one bit of the u32 \
         annotation slot, so at most {max} are supported per file"
    )]
    TooManyFeatures {
        path: String,
        count: usize,
        max: usize,
    },
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
    #[error(
        "PAM defines {count} concrete variants (plen={plen}); \
             exceeds the {max} addressable by a u16 index"
    )]
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
    #[error(
        "contig name {name:?} (id {id}) contains illegal byte {byte} \
             (one of ',' '\"' '\\n' '\\r'), which would break the CSV"
    )]
    InvalidName { id: u32, name: String, byte: u8 },
}

#[derive(Debug, Error)]
pub enum PartitionError {
    #[error("{context}: {source}")]
    Io {
        context: String,
        #[source]
        source: std::io::Error,
    },

    #[error("intermediate report has no header line")]
    MissingHeader,

    #[error(
        "header has {cols} columns; expected {FIXED_COLS}..={} \
         ({FIXED_COLS} fixed + 0..{MAX_ANNOTATIONS} annotation columns)",
        FIXED_COLS + MAX_ANNOTATIONS
    )]
    BadHeaderShape { cols: usize },

    #[error("line {line} of the intermediate report is not valid UTF-8")]
    NotUtf8 { line: usize },

    #[error("line {line} has {found} columns; the header declared {expected}")]
    FieldCount {
        line: usize,
        expected: usize,
        found: usize,
    },

    #[error("line {line}, column {col}: cannot parse {what} from {value:?}")]
    BadField {
        line: usize,
        col: usize,
        value: String,
        what: &'static str,
    },
}

impl PartitionError {
    /// Build an [`PartitionError::Io`] with a human-readable context.
    pub(crate) fn io(context: impl Into<String>, source: std::io::Error) -> Self {
        Self::Io {
            context: context.into(),
            source,
        }
    }
}

/// Errors raised while parsing a chain file or applying liftover to a report.
///
/// Maps to Python through [`crate::python::pyerrors`]: `Io` -> `IOError`, the
/// rest -> `ValueError`.
#[derive(Debug, Error)]
pub enum LiftoverError {
    #[error("{context}: {source}")]
    Io {
        context: String,
        #[source]
        source: std::io::Error,
    },

    #[error("malformed chain file at line {line}: {msg}")]
    MalformedChain { line: usize, msg: String },

    #[error("intermediate report has no header line")]
    MissingHeader,

    #[error("report is missing required column {0:?}")]
    MissingColumn(String),

    #[error("report line {line}: {msg}")]
    BadRow { line: usize, msg: String },
}

impl LiftoverError {
    pub(crate) fn io(context: impl Into<String>, source: std::io::Error) -> Self {
        Self::Io {
            context: context.into(),
            source,
        }
    }
    pub(crate) fn malformed_chain(line: usize, msg: impl Into<String>) -> Self {
        Self::MalformedChain {
            line,
            msg: msg.into(),
        }
    }
    pub(crate) fn bad_row(line: usize, msg: impl Into<String>) -> Self {
        Self::BadRow {
            line,
            msg: msg.into(),
        }
    }
}
