use pyo3::{pyclass, pymethods};

use crate::alignment::thresholds::Thresholds;
use crate::crispr::guide::Guide;

#[derive(Debug, Clone)]
#[pyclass]
pub struct AlignmentParams {
    /// Length of the input sequences
    #[pyo3(get, set)]
    pub sequence_len: usize,
    /// Number of sequences in a batch
    #[pyo3(get, set)]
    pub sequence_batch_size: usize,
    /// Number of alignments in output batch
    #[pyo3(get, set)]
    pub alignment_batch_size: usize,
    /// Miner thresholds
    #[pyo3(get, set)]
    pub thresholds: Thresholds,
    /// Maximum mutation score
    #[pyo3(get, set)]
    pub mutation_max: u32,
    /// Miner guide
    #[pyo3(get, set)]
    pub guide: Guide,
}

#[pymethods]
impl AlignmentParams {
    #[new]
    pub fn new(
        sequence_len: usize,
        sequence_batch_size: usize,
        alignment_batch_size: usize,
        thresholds: Thresholds,
        mutation_max: u32,
        guide: Guide,
    ) -> Self {
        Self {
            sequence_len,
            sequence_batch_size,
            alignment_batch_size,
            thresholds,
            mutation_max,
            guide,
        }
    }
}
