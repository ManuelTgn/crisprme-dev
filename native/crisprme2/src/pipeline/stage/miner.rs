use columnar::{MemoryPool, pipeline::{Emit, Stage, StageError}, Share};
use itertools::izip;
use rand::Rng;

use crate::model::{alignment::{SeqMinedBatch, SeqMinedFrame}, cigarx::{Cigarx, Cigarx64, CigarxOp}, input::SeqBatch};

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


/*
/// Fake miner for now, it only checks for match/mismatch
pub struct MineScanner {
    /// Pool for buffers of mined schema
    pool: Arc<Pool<MinedSchema>>
}

impl MineScanner {
    pub fn new(pool: Arc<Pool<MinedSchema>>) -> Self {
        Self { pool }
    }
}

impl Stage for MineScanner {
    
    type Input  = BatchRef<SeqSchema, SeqBatchMetadata>;
    type Output = BatchMut<MinedSchema, MinedBatchMetadata>;

    fn process<E>(&mut self, input: Self::Input, emitter: &mut E) -> Result<(), StageError>
    where
        E: Emit<Self::Output> 
    {
        use crate::model::input::sequences::schema as ss;
        use crate::model::alignment::mined::schema as ms;

        println!("[MineScanner] received buffer");
        let (seq_ids, seq_contents) = input.columns((ss::id, ss::content));
        let guide = input.metadata.guide.as_slice();

        assert!(guide.len() <= input.metadata.seq_len);
        
        let mut remaining = input.len();
        while remaining > 0 {
            let mut result = self.pool.acquire()
                .map_err(|_| StageError)?;

            let rows = remaining.min(input.len());
            result.set_len(rows);
            result.mutate(
                (ms::seq_id, ms::cigarx, ms::offset),
                |(mined_seq_ids, mined_cigarxs, mined_offsets)| {
                    for i in 0..rows {

                        mined_seq_ids[i] = seq_ids[i];
                        mined_offsets[i] = i as u8;

                        let mut cigarx = Cigarx64::default();
                        for j in 0..guide.len() {
                            if seq_contents[i][j].matches(guide[j]) {
                                cigarx.push(CigarxOp::Match);
                            } else {
                                cigarx.push(CigarxOp::Mismatch);
                            }
                        }

                        mined_cigarxs[i] = cigarx;
                    }
                }
            );

            remaining -= result.len();
            println!("[MineScanner] submitted output buffer with {} rows", result.len());
            emitter.emit(result.with_metadata(
                MinedBatchMetadata { 
                    sequences: input.clone() 
                }
            ))?;
        }
        Ok(())
    }
}
 */