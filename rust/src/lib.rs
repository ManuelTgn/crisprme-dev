// modules used by the main function
mod scan;
mod pam; 
mod iupac;
mod target;
mod hashing;

use std::collections::HashMap;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3::exceptions::PyValueError;
use pyo3::{PyResult, PyObject};

// Type alias for the complex value data, just for cleaner code
type OccurrenceData = Vec<(String, usize, bool)>;

/// Finds all potential target candidates (CRISPR gRNAs) within a given sequence.
///
/// This function converts the input sequence and PAM into IUPAC bitmasks and performs a
/// parallelized scan to identify all positions where the target sequence and its 
/// associated PAM match the defined criteria.
///
/// # Arguments
/// * `sequence` (str): The large DNA/RNA sequence to be scanned (e.g., a contig or chromosome).
/// * `contig` (str): The name/identifier of the sequence (e.g., "chr1").
/// * `pam_seq` (str): The Protospacer Adjacent Motif (PAM) sequence (e.g., "NGG").
/// * `k` (usize): The length of the target/protospacer sequence, excluding the PAM.
/// * `right` (bool): If `true`, the PAM is expected to be immediately *right* of the target sequence.
///                   If `false`, the PAM is expected to be immediately *left* of the target sequence.
/// * `threads` (usize): The number of threads to use for parallel scanning.
///
/// # Returns
/// A `list` of `target::Target` objects, where each object contains the position, 
/// orientation, and bitmask sequence of a found target.
///
/// # Errors
/// Returns a `PyValueError` if input constraints are violated (e.g., invalid sizes or PAM sequence).
#[pyfunction]
pub fn find_target_candidates(
    py: Python,
    sequence: &str, 
    contig: &str, 
    pam_seq: &str, 
    k: usize, 
    right: bool,
    path: &str,
    threads: usize,
// ) -> PyResult<PyObject> {
) -> PyResult<()> {
    // --- input validation ---

    if k == 0 {
        return Err(PyValueError::new_err(
            "Size must be greater than 0",
        ));
    }
    
    if threads == 0 {
        return Err(PyValueError::new_err(
            "threads must be greater than 0",
        ));
    }
    
    let seq_len = sequence.len();
    
    if seq_len == 0 {
        return Err(PyValueError::new_err(
            "sequence cannot be empty",
        ));
    }
    
    if k > seq_len {
        return Err(PyValueError::new_err(
            format!("size ({}) cannot be greater than sequence length ({})", k, seq_len),
        ));
    }

    // --- PAM parsing and validation ---

    // parse the PAM sequence string into a ParsedPAM struct (converting to bitmasks)
    // the .map_err converts the internal Rust String error into a Python PyValueError
    let pat = pam::ParsedPAM::new(pam_seq)
        .map_err(|e| PyErr::new::<PyValueError, _>(format!("Invalid PAM sequence: {}", e)))?;

    // -- execution ---

    // // execute the core parallel scanning logic and return the results
    scan::scan_targets(sequence, contig, &pat, k, right, path, threads);

    Ok(())
    // let py_dict = PyDict::new(py);
    
    // // Iterate over the Rust HashMap
    // for (vec_u8_key, occurrences) in targets_map.into_iter() {
        
    //     // Convert Vec<u8> (Rust) to Python bytes (immutable, hashable)
    //     // .as_bytes() creates the Python bytes object from the Rust slice.
    //     let py_key = vec_u8_key.as_slice().into_py(py);
        
    //     // Convert the OccurrenceData (Vec<...>) to a Python list of tuples
    //     // PyO3 handles the Vec -> list and inner tuple conversion automatically here.
    //     let py_value = occurrences.into_py(py);
        
    //     // Insert the key (bytes) and value (list) into the dictionary
    //     // This is where the panic occurred previously, but now the key is guaranteed 'bytes'.
    //     py_dict.set_item(py_key, py_value)?;
    // }

    // // Return the dictionary as a PyObject
    // Ok(py_dict.into())
}

/// Defines the Python module structure and exposes Rust functions
#[pymodule]
fn target_candidates_parser(_py: Python, m : &PyModule) -> PyResult<()> {
    // add the top-level function to the Python module
    m.add_function(wrap_pyfunction!(find_target_candidates, m)?)?;
    
    Ok(())
}
