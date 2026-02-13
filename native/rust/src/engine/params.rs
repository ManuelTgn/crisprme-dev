
use crate::alignment::thresholds::Thresholds;
use crate::crispr::guide::Guide;

#[derive(Debug, Clone)]
pub struct AlignmentParams {
    
    /// Length of the input sequences
    pub sequence_len: usize,
    /// Number of sequences in a batch
    pub sequence_batch_size: usize,
    /// Number of alignments in output batch
    pub alignment_batch_size: usize,
    /// Miner thresholds
    pub thresholds: Thresholds,
    /// Maximum mutation score
    pub mutation_max: u32,
    /// Miner guide
    pub guide: Guide,

}