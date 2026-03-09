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