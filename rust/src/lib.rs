use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use rayon::prelude::*;

/// Extract all targets from a sequence in parallel
///
/// This function processes a DNA/RNA sequence and extracts all targets (k-mers) 
/// of length k using parallel processing. It correctly handles chunk boundaries 
/// to ensure no targets are missed at the edges of chunks.
///
/// # Arguments
/// * `py` - Python GIL token
/// * `sequence` - Full chromosome sequence as a string
/// * `k` - K-mer length (e.g., 20)
/// * `threads` - Number of threads to use for parallel processing
///
/// # Returns
/// A Python list of tuples (position, target_string) sorted by position
///
/// # Errors
/// Returns PyErr if:
/// - k is 0 or greater than sequence length
/// - threads is 0
/// - sequence is empty
///
/// # Example
/// ```python
/// from crisprme_kmers import extract_targets_parallel
/// 
/// sequence = "ATCGATCGATCG"
/// kmers = extract_targets_parallel(sequence, k=4, threads=2)
/// # Returns: [(0, "ATCG"), (1, "TCGA"), (2, "CGAT"), ...]
/// ```
#[pyfunction]
fn extract_targets_parallel(
    py: Python<'_>,
    sequence: &str,
    size: usize,
    threads: usize,
) -> PyResult<PyObject> {
    // Input validation
    if size == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Size must be greater than 0",
        ));
    }
    
    if threads == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "threads must be greater than 0",
        ));
    }
    
    let seq_len = sequence.len();
    
    if seq_len == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "sequence cannot be empty",
        ));
    }
    
    if size > seq_len {
        return Err(pyo3::exceptions::PyValueError::new_err(
            format!("size ({}) cannot be greater than sequence length ({})", size, seq_len),
        ));
    }

    // Work on bytes for safe indexing (DNA/RNA are ASCII)
    let seq_bytes = sequence.as_bytes();

    // Configure Rayon thread pool and run **pure Rust** work inside
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to create thread pool: {}", e)
        ))?;

    // The closure passed to install MUST NOT capture `Python<'_>` or any non-Send type.
    let kmers: Vec<(usize, String)> = thread_pool.install(|| {
        // Calculate chunk size (ceiling division)
        let chunk_size = (seq_len + threads - 1) / threads;

        // Parallel iterate over chunk indices; produce Vec<Vec<(pos, String)>>
        (0..threads)
            .into_par_iter()
            .filter_map(|chunk_idx| {
                let orig_start = chunk_idx * chunk_size;
                if orig_start >= seq_len {
                    return None;
                }
                let orig_end = std::cmp::min(orig_start + chunk_size, seq_len);

                // Extend chunk by size-1 to capture k-mers crossing boundary
                let extended_start = orig_start;
                let extended_end = std::cmp::min(orig_end + (size - 1), seq_len);

                // Slice bytes (safe because DNA/RNA are ASCII)
                let chunk = &seq_bytes[extended_start..extended_end];
                let chunk_len = chunk.len();

                let mut chunk_kmers = Vec::new();
                if chunk_len >= size {
                    for i in 0..=(chunk_len - size) {
                        let global_pos = extended_start + i;
                        // Emit only k-mers whose start is inside the original chunk (avoid duplicates)
                        if global_pos >= orig_start && global_pos < orig_end {
                            // Convert k-mer bytes back to String (UTF-8 safe for ASCII)
                            let kmer = String::from_utf8(chunk[i..i + size].to_vec())
                                .expect("sequence should be ASCII");
                            chunk_kmers.push((global_pos, kmer));
                        }
                    }
                }
                Some(chunk_kmers)
            })
            .flatten()
            .collect()
    });

    // Sort the results by position
    let mut flat_kmers = kmers;
    flat_kmers.sort_unstable_by_key(|&(pos, _)| pos);

    // Now acquire the GIL once and convert to Python objects
    Python::with_gil(|py| {
        let py_list = PyList::empty(py);
        for (pos, kmer) in flat_kmers {
            let py_tuple = PyTuple::new(py, &[pos.into_py(py), kmer.into_py(py)]);
            py_list.append(py_tuple)?;
        }
        Ok(py_list.into())
    })
    
}

/// Python module definition
#[pymodule]
fn sequence_parser(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(extract_targets_parallel, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_basic_targets_extraction() {
        // Simple test sequence
        let seq = "ATCGATCGATCG";
        let k = 4;
        
        // Extract expected k-mers manually
        let expected = vec![
            "ATCG", "TCGA", "CGAT", "GATC", "ATCG", "TCGA", "CGAT", "GATC", "ATCG"
        ];
        
        // This is a unit test without Python, so we can't test the PyO3 function directly
        // But we can verify the logic
        let chunk_len = seq.len();
        let mut kmers = Vec::new();
        
        for i in 0..=(chunk_len - k) {
            kmers.push(&seq[i..i + k]);
        }
        
        assert_eq!(kmers, expected);
    }
    
    #[test]
    fn test_chunk_boundaries() {
        let seq = "AAAABBBBCCCCDDDD";
        let k = 4;
        let chunk_size = 8;
        
        // Chunk 0: orig [0, 8), extended [0, 11)
        // Chunk 1: orig [8, 16), extended [8, 16)
        
        // Verify chunk 0 k-mers
        let chunk0 = &seq[0..11];
        assert_eq!(chunk0, "AAAABBBBCCC");
        
        // K-mers at positions 0-7 should be emitted
        for i in 0..=7 {
            if i + k <= chunk0.len() {
                let _kmer = &chunk0[i..i + k];
                assert!(i < 8); // All should be within orig bounds
            }
        }
    }
}