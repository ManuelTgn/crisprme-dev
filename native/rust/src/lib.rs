// modules used by the main function
mod bindings;
mod crispr;
mod utils;
mod memory;
mod alignment;
mod sequence;
mod batching;
mod storage;
mod engine;
mod error;
mod annotation;
pub mod python;

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::PyResult;

use crate::{alignment::thresholds::Thresholds, batching::batching::{BatcherStats, FeedStatus, TargetBatcher}, crispr::guide::Guide, engine::{hybrid::HybridEngine, params::AlignmentParams}, python::views::AlignmentBatchView};


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
    let pat = crispr::pam::ParsedPAM::new(pam_seq)
        .map_err(|e| PyErr::new::<PyValueError, _>(format!("Invalid PAM sequence: {e}")))?;

    // Execute the core parallel scanning logic and return the results
    sequence::scanner::scan_targets(sequence, &pat, size, right, threads)
        .map_err(|e| PyErr::new::<PyValueError, _>(e))
}

#[pyfunction]
pub fn initialize_engine_logger() {
    tracing_subscriber::fmt()
            .compact()
            .with_target(false)
            .with_thread_ids(true)
            .with_max_level(tracing::Level::TRACE)
            .init();
}

/// Defines the Python module structure and exposes Rust functions
#[pymodule]
fn _crisprme2_native(_py: Python, m : &PyModule) -> PyResult<()> {
    // add the top-level function to the Python module
    // m.add_function(wrap_pyfunction!(extract_targets_rs, m)?)?;

    m.add_function(wrap_pyfunction!(initialize_engine_logger, m)?)?;

    m.add_class::<TargetBatcher>()?;
    m.add_class::<FeedStatus>()?;
    m.add_class::<BatcherStats>()?;
    m.add_class::<HybridEngine>()?;
    m.add_class::<AlignmentParams>()?;
    m.add_class::<Thresholds>()?;
    m.add_class::<Guide>()?;
    m.add_class::<AlignmentBatchView>()?;
    
    Ok(())
}
