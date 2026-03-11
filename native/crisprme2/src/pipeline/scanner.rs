use std::sync::Arc;

use columnar::{pipeline::{Emit, Stage, StageError}, pool::{BatchMut, BatchRef, Pool}};

use crate::model::{alignment::{MinedBatchMetadata, MinedSchema}, cigarx::{Cigarx, Cigarx64, CigarxOp}, input::{SeqBatchMetadata, SeqSchema}};

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