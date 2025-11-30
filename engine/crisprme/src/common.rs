use clap::Args;

/*
/// Configuration for the parallelization of the mining process
#[derive(Debug, Args)]
pub struct ThreadsConfig {

    /// Fixed mining chunk size, overwrites everything
    #[arg(long)]
    pub chunk_size: Option<usize>,

    /// Maximum size of the mining chunks
    #[arg(long, default_value_t = 1000)]
    pub max_chunk_size: usize,

    /// Minimum size of the mining chunks
    #[arg(long, default_value_t = 10000)]
    pub min_chunk_size: usize,
}
*/

/// Threholds for the mining process
#[derive(Debug, Args, Clone)]
pub struct CliThresholds {

    /// Max allowed gaps in query
    #[arg(long)]
    pub qgap: u32,

    /// Max allowed gaps in target
    #[arg(long)]
    pub tgap: u32,

    /// Max allowed mismatches
    #[arg(long)]
    pub mism: u32,
}

/// Limit the amout of memory usage
#[derive(Debug, Args)]
pub struct MemoryConfig {

    /// Number of sequences in input batch
    #[arg(short, long)]
    #[clap(default_value_t = 1_000_000)]
    pub sequence_batch_size: usize,

    /// Number of alignments in output batch
    #[arg(short, long)]
    #[clap(default_value_t = 10_000_000)]
    pub alignment_batch_size: usize
}

use crisprme_core::utils::Thresholds;
impl Into<Thresholds> for &CliThresholds {
    fn into(self) -> Thresholds {
        Thresholds { 
            qgap: self.qgap, 
            tgap: self.tgap, 
            mism: self.mism 
        }
    }
}

impl Into<crisprme_core::common::thresholds::Thresholds> for &CliThresholds {
    fn into(self) -> crisprme_core::common::thresholds::Thresholds {
        crisprme_core::common::thresholds::Thresholds { 
            qgap: self.qgap, 
            tgap: self.tgap, 
            mism: self.mism 
        }
    }
}
