// modules used by the main function
mod scan;
mod pam; 
mod iupac;
mod target;
mod hashing;

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::PyResult;

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
    sequence: &str, 
    contig: &str, 
    pam_seq: &str, 
    k: usize, 
    right: bool,
    threads: usize,
) -> PyResult<hashing::HashedTargets> {
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

    // execute the core parallel scanning logic and return the results
    Ok(scan::scan_targets(sequence, contig, &pat, k, right, threads))
}

/// Defines the Python module structure and exposes Rust functions
#[pymodule]
fn target_candidates_parser(_py: Python, m : &PyModule) -> PyResult<()> {
    // add the top-level function to the Python module
    m.add_function(wrap_pyfunction!(find_target_candidates, m)?)?;
    
    Ok(())
}
