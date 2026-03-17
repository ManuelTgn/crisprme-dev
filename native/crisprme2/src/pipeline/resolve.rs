
use std::collections::HashMap;

use columnar::{MemoryPool, pipeline::{Emit, Stage, StageError}};
use itertools::izip;
use crate::model::{alignment::{SeqMinedBatch, SeqResolvedBatch, SeqResolvedFrame}, cigarx::{Cigarx, CigarxOp}};

/// Resolve mined alignments using the present cigarx
pub struct Resolver {
    map: HashMap<u32, (usize, usize)>,
    pool: MemoryPool,
}

impl Resolver {
    pub fn new(pool: &MemoryPool) -> Self {
        Self { 
            pool: pool.clone(),
            map: HashMap::new()
        }
    }
}

impl Stage for Resolver {

    type I = SeqMinedBatch;
    type O = SeqResolvedBatch;

    fn name() -> &'static str { "Resolver" }
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError> {
        
        // Create a map from seq_id to position in memory
        // TODO: We can store this at the beginning of the pipeline (?)
        input.sequences.with_cols(|cols| {
            self.map.clear();
            for (i, seq_id) in cols.id.iter().enumerate() {
                self.map.insert(*seq_id, 
                    cols.content.index(i));
            }
        });

        // mined --1:1--> resolved
        let guide = input.guide.as_slice();
        input.sequences.with_cols(|sequences| {
            input.mined.with_cols(|mut mined| {

                let rows = mined.seq_id.rows();
                let mut resolved = SeqResolvedFrame::empty();
                resolved.with_cols(|mut resolved| {

                    // Share columns (seq_id, offset)
                    resolved.seq_id.shared(&mut mined.seq_id);
                    resolved.offset.shared(&mut mined.offset);

                    // Allocate columns (rguide, rseq)
                    resolved.rguide.alloc(&self.pool, rows);
                    resolved.rseq.alloc(&self.pool, rows);

                    // Zipped iterator over all used columns
                    let zipper = izip!(
                        resolved.seq_id.iter(),
                        resolved.rguide.iter_mut(),
                        resolved.rseq.iter_mut(),
                        mined.cigarx.iter(),
                        mined.offset.iter()
                    );

                    // Resolve the guide and sequence
                    for (seq_id, rguide, rseq, cigarx, offset) in zipper {

                        // Indirect look-up to sequence content
                        let (seq_chunk, seq_offset) = self.map[seq_id];
                        let sequence = sequences.content.get_fast(seq_chunk, seq_offset);

                        let mut gpos = 0usize;
                        let mut spos = *offset as usize;  // start at alignment position in sequence
                        let mut opos = 0usize;

                        for op in cigarx.iter() {
                            match op {
                                CigarxOp::Match | CigarxOp::Mismatch => {
                                    rguide[opos] = guide[gpos].to_ascii();
                                    rseq[opos]   = sequence[spos].to_ascii();
                                    gpos += 1; spos += 1;
                                }
                                CigarxOp::Deletion => {
                                    rguide[opos] = b'-';
                                    rseq[opos]   = sequence[spos].to_ascii();
                                    spos += 1;
                                }
                                CigarxOp::Insertion => {
                                    rguide[opos] = guide[gpos].to_ascii();
                                    rseq[opos]   = b'-';
                                    gpos += 1;
                                }
                            }
                            opos += 1;
                        }

                        // Null-terminate both resolved arrays
                        rguide[opos] = 0;
                        rseq[opos]   = 0;

                    }
                });

                emitter.emit(SeqResolvedBatch { 
                    occurences: input.occurences, 
                    resolved
                })
            })
        })
    }
}