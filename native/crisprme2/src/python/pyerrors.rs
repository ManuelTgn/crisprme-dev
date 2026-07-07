use crate::error::crisprme_errors::AnnotationError;

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::PyErr;

impl From<AnnotationError> for PyErr {
    fn from(err: AnnotationError) -> PyErr {
        match err {
            AnnotationError::IoError(msg) =>
                PyIOError::new_err(msg),

            AnnotationError::MalformedBed { line } =>
                PyValueError::new_err(
                    format!("Malformed BED file at line {}", line)
                ),
            
            AnnotationError::InvalidFeatureId(id) =>
                PyValueError::new_err(
                    format!("Invalid feature ID {}", id)
                ),

            AnnotationError::EmptyInput =>
                PyValueError::new_err("Annotation input cannot be empty"),
        }
    }
}


impl From<PamError> for PyErr {
    fn from(err: PamError) -> PyErr {
        match err {
            // Bad input / configuration -> ValueError
            PamError::InvalidCharacter { .. } | PamError::TooManyVariants { .. } =>
                PyValueError::new_err(err.to_string()),
            // Out-of-range lookup -> IndexError (more Pythonic)
            PamError::IndexOutOfRange { .. } =>
                PyIndexError::new_err(err.to_string()),
        }
    }
}