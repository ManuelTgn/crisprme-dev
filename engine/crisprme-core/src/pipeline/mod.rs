use std::path::PathBuf;

use crate::common::{guide::Guide, thresholds::Thresholds};

pub mod engine;
pub mod miner;

#[derive(Debug, Clone)]
pub struct PipelineDescriptor {
    
    /// Length of the input sequences
    pub sequence_len: usize,
    /// Number of sequences in a batch
    pub sequence_batch_size: usize,
    /// Number of alignments in output batch
    pub alignment_batch_size: usize,

    /// File containing all the sequences
    pub sequence_file: PathBuf,
    /// Output file for alignments
    pub output_file: PathBuf,
    
    /// Miner thresholds
    pub thresholds: Thresholds,
    /// Maximum mutation score
    pub mutation_max: u32,
    /// Miner guide
    pub guide: Guide,

}
