// modules used by the main function
mod scan;
mod pam; 
mod iupac;
mod target;
mod threadpool;

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
pub fn extract_targets_rs(
    sequence: &str, 
    pam_seq: &str,
    size: usize, 
    right: bool,
    threads: usize,
) -> PyResult<(Vec<usize>, Vec<u8>)> {
    let pat = pam::ParsedPAM::new(pam_seq)
        .map_err(|e| PyErr::new::<PyValueError, _>(format!("Invalid PAM sequence: {e}")))?;

    // Execute the core parallel scanning logic and return the results
    scan::scan_targets(sequence, &pat, size, right, threads)
        .map_err(|e| PyErr::new::<PyValueError, _>(e))
}


/// Defines the Python module structure and exposes Rust functions
#[pymodule]
fn target_candidates_scanner_rs(_py: Python, m : &PyModule) -> PyResult<()> {
    // add the top-level function to the Python module
    m.add_function(wrap_pyfunction!(extract_targets_rs, m)?)?;
    
    Ok(())
}
