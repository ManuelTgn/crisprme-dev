use columnar::{MemoryPool, pipeline::{Emit, Stage, StageError}, Share};
use crossbeam_channel::Receiver;
use itertools::izip;
use rand::Rng;

use crate::model::{alignment::{SeqMinedBatch, SeqMinedFrame}, cigarx::{Cigarx, Cigarx64, CigarxOp}, input::SeqBatch};

// ---------------------------------------------------------------------------
// Fake miner
// ---------------------------------------------------------------------------

/// Fake miner for now, it only checks for match/mismatch
pub struct Miner { pool: MemoryPool }

impl Miner {
    pub fn new(pool: &MemoryPool) -> Self {
        Self { pool: pool.clone() }
    }
}

impl Stage for Miner {

    type I = SeqBatch;
    type O = SeqMinedBatch;

    fn name() -> &'static str { "Miner" }

    #[tracing::instrument(name = "pipeline:miner", skip_all)]
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError> {
        
        let mut sequences  = input.sequences.share();
        let mut occurences = input.occurences.share();

        input.sequences.with_cols(|cols| {

            let rows = cols.content.rows();
            tracing::info!("received {} rows to mine", rows);

            let mut mined = SeqMinedFrame::alloc(&self.pool, rows);
            mined.with_cols(|mut mined| {
                let zipped = izip!(
                    cols.content.iter(),
                    mined.seq_row_idx.iter_mut(),
                    mined.cigarx.iter_mut(),
                    mined.offset.iter_mut()
                );

                for (sequence, mined_seq_row, cigarx, offset) in zipped {
                     
                    *mined_seq_row = rand::rng().random_range(0..rows) as u32;

                    *cigarx = Cigarx64::default();
                    *offset = 2;

                    for j in 0..input.guide.len() {
                        if j % 3 == 0 {
                            cigarx.push(CigarxOp::Deletion);
                        } else if j % 5 == 0 {
                            cigarx.push(CigarxOp::Insertion);
                        } else if sequence[j].matches(input.guide[j]) {
                            cigarx.push(CigarxOp::Match);
                        } else {
                            cigarx.push(CigarxOp::Mismatch);
                        }
                    }
                }

            });

            emitter.emit(SeqMinedBatch {
                guide: input.guide.clone(),
                sequences: sequences.share(),
                occurences: occurences.share(),
                mined,
            })
        })
    }
}

// ---------------------------------------------------------------------------
// GPU miner
// ---------------------------------------------------------------------------

pub struct GpuMiner {
    pool: MemoryPool,
    gpu: usize,
}

impl GpuMiner {
    pub fn new(pool: &MemoryPool, gpu: usize) -> Self {
        Self { 
            pool: pool.clone(), 
            gpu 
        }
    }
}

impl Stage for GpuMiner {

    type I = SeqBatch;
    type O = SeqMinedBatch;

    fn name() -> &'static str { "GpuMiner" }

    // Overload run function
    #[tracing::instrument(name = "pipeline:gpu_miner", skip_all)]
    fn process(&mut self, input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError> {
        todo!()
    }
}