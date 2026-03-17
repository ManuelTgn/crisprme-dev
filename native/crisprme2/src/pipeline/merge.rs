use std::sync::Arc;

use bumpalo::collections::Vec;
use columnar::{
    pipeline::{Emit, Stage, StageError},
    MemoryPool,
};
use itertools::izip;

use crate::model::alignment::{AlignmentFrame, SeqResolvedBatch};

/// Broadcasts resolved alignments to sequence occurences
pub struct MergeJoinSorted {
    pool: MemoryPool,
}

/// Invariant: all ids are monotonically increasing
impl MergeJoinSorted {
    pub fn new(pool: &MemoryPool) -> Self {
        Self { pool: pool.clone() }
    }
}

impl Stage for MergeJoinSorted {
    type I = SeqResolvedBatch;
    type O = AlignmentFrame;

    fn name() -> &'static str { "MergeJoinSorted" }
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError> {

        input.resolved.with_cols(|resolved| {
            input.occurences.with_cols(|mut occurence| {

                let rows = occurence.seq_id.rows();
                let mut alignment = AlignmentFrame::empty();
                alignment.with_cols(|mut alignment| {

                    // Shared columns
                    alignment.seq_id.shared(&mut occurence.seq_id);
                    alignment.occurence.shared(&mut occurence.occurence);

                    // New columns required
                    alignment.offset.alloc(&self.pool, rows);
                    alignment.rguide.alloc(&self.pool, rows);
                    alignment.rseq.alloc(&self.pool, rows);

                    // Merged iterators on all columns of resolved
                    let mut resolved_iter = izip!(
                        resolved.seq_id.iter(),         // 0 (Primary key)
                        resolved.offset.iter(),         // 1
                        resolved.rguide.iter(),         // 2
                        resolved.rseq.iter()            // 3
                    );

                    // Merged iterators on all columns of alignment
                    let alignment_iter = izip!(
                        alignment.seq_id.iter(),        // 0 (Foreign key)
                        alignment.offset.iter_mut(),    // 1
                        alignment.rguide.iter_mut(),    // 2
                        alignment.rseq.iter_mut()       // 3
                    );

                    // Cursors at the resolved elements
                    let mut resolved_curr = resolved_iter.next().unwrap();
                    for alignment_curr in alignment_iter {
                        
                        // Advance the resolved cursors
                        while *resolved_curr.0 < *alignment_curr.0 {
                            resolved_curr = resolved_iter.next().unwrap();
                        }

                        *alignment_curr.1 = *resolved_curr.1;
                        *alignment_curr.2 = *resolved_curr.2;
                        *alignment_curr.3 = *resolved_curr.3;
                    }
                });
                emitter.emit(alignment)
            })
        })
    }
}