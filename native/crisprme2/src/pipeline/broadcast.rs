use std::sync::Arc;

use columnar::{arena::Arena, pipeline::{Emit, StageError, Stage}, pool::{BatchMut, BatchRef, Pool}};

use crate::model::alignment::{ALIGN_MAX_FEATURES, ALIGN_MAX_SCORES, AlignmentSchema, ResolvedBatchMetadata, ResolvedSchema, resolved};

/// Broadcasts resolved alignments to sequence occurences
/// It uses a `gather` approach to keep writes sequential and reads random.
pub struct AlignmentBroadcast {
    
    /// Pool for buffers of alignment schema
    pool: Arc<Pool<AlignmentSchema>>,
    /// Temporary buffer
    arena: Arena,
}

impl AlignmentBroadcast {
    pub fn new(pool: Arc<Pool<AlignmentSchema>>, memory: usize) -> Self {
        Self {
            arena: Arena::with_capacity(memory),
            pool
        }
    }
}

impl Stage for AlignmentBroadcast {

    type Input  = BatchMut<ResolvedSchema, ResolvedBatchMetadata>;
    type Output = BatchMut<AlignmentSchema, ()>;

    fn process<E>(&mut self, input: Self::Input, emitter: &mut E) -> Result<(), StageError>
    where
        E: Emit<Self::Output>
    {
        use crate::model::input::occurences::schema   as os;
        use crate::model::alignment::resolved::schema as rs;
        use crate::model::alignment::aligned::schema  as bs;

        let _span = tracing::debug_span!("alignment-broadcast")
            .entered();

        self.arena.scoped(|m| {
        
            // Map from seq_id to index inside the input batch
            // NOTE: seq_id is always from 0 to input.len()
            let mut index = m.alloc_slice_fill(input.len(), 0u32);
            let (seq_ids,) = input.columns((resolved::schema::seq_id,));
            for (i, id) in seq_ids.iter().enumerate() { 
                index[*id as usize] = i as u32;
            }

            let occurence_count = input.metadata.occurences.len();
            tracing::debug!("broadcast {} rows to {} occurence batches", 
                input.len(), occurence_count);

            for (i, batch) in input.metadata.occurences.iter().enumerate() {

                let (occ_seq_ids, occ_occurences) = batch.columns((os::seq_id, os::occurence));
                let (res_rguides, res_rseqs, res_resolved_lens, res_offsets) = 
                    input.columns((rs::rguide, rs::rseq, rs::resolved_len, rs::offset));

                let mut remaining = batch.len();
                tracing::debug!("{} rows to emit in total", remaining);
                while remaining > 0 {

                    // Acquire a new result batch
                    let mut result = self.pool.acquire()
                        .map_err(|_| StageError)?;

                    // How many rows do we need to fill in the result batch
                    let occ_offset = batch.len() - remaining;  // must be before decrement
                    let rows = remaining.min(result.capacity());
                    result.set_len(rows);

                    result.mutate(
                        (bs::id, bs::rguide, bs::rseq, bs::resolved_len, bs::offset, bs::occurence, bs::features, bs::scores), 
                        |(ids, rguides, rseqs, resolved_lens, offsets, occurences, features, scores)| {

                            // Gather all data to the result batch
                            for i in 0..rows {

                                let occ_row = occ_offset + i;
                                let res_row = index[occ_seq_ids[occ_row] as usize] as usize;

                                ids[i]           = 0; // filled by a later stage
                                rguides[i]       = res_rguides[res_row];
                                rseqs[i]         = res_rseqs[res_row];
                                resolved_lens[i] = res_resolved_lens[res_row];
                                offsets[i]       = res_offsets[res_row];
                                occurences[i]    = occ_occurences[occ_row];
                            }

                            // To be computed at a later stage
                            for f in 0..ALIGN_MAX_FEATURES { features[f].fill(0); }
                            for s in 0..ALIGN_MAX_SCORES   { scores[s].fill(0.0); }
                        });

                    remaining -= result.len();
                    emitter.emit(result)?;
                    tracing::debug!("emit alignment batch with {} rows ({}/{} remaining)", 
                        rows, remaining, batch.len());
                }
            }

            Ok(())
        })
    }
}

#[cfg(test)]
mod tests {
    use columnar::pool::connector_mut;
    use crate::model::{alignment::aligned, input::{SeqOccSchema, occurences}};
    use super::*;

    /// Helper: build a stage with given output pool capacity (rows per batch).
    fn make_stage(out_capacity: usize) -> AlignmentBroadcast {
        AlignmentBroadcast {
            pool: Arc::new(Pool::<AlignmentSchema>::new(16, out_capacity)),
            arena: Arena::with_capacity(1 << 16),
        }
    }

    /// Verify that the gather correctly maps each occurrence to the right resolved row.
    ///
    /// Resolved:    seq_id=0 → resolved_len=20, offset=0
    ///              seq_id=1 → resolved_len=25, offset=5
    /// Occurrences: [(seq_id=1, occ=100), (seq_id=0, occ=200), (seq_id=1, occ=300)]
    ///
    /// Expected output rows:
    ///   row 0: resolved_len=25, offset=5, occurence=100
    ///   row 1: resolved_len=20, offset=0, occurence=200
    ///   row 2: resolved_len=25, offset=5, occurence=300
    #[test]
    fn gather_maps_resolved_to_occurrences() {

        let mut stage     = make_stage(16);
        let resolved_pool = Pool::<ResolvedSchema>::new(4, 16);
        let occ_pool      = Pool::<SeqOccSchema>::new(4, 16);

        let mut resolved = resolved_pool.acquire().unwrap();
        resolved.set_len(2);
        resolved.mutate(
            (resolved::schema::seq_id, resolved::schema::resolved_len, resolved::schema::offset),
            |(seq_ids, lens, offsets)| {
                seq_ids[0] = 0; lens[0] = 20; offsets[0] = 0;
                seq_ids[1] = 1; lens[1] = 25; offsets[1] = 5;
            },
        );

        let mut occ = occ_pool.acquire().unwrap();
        occ.set_len(3);
        occ.mutate(
            (occurences::schema::seq_id, occurences::schema::occurence),
            |(seq_ids, occs)| {
                seq_ids[0] = 1; occs[0] = 100;
                seq_ids[1] = 0; occs[1] = 200;
                seq_ids[2] = 1; occs[2] = 300;
            },
        );

        let input = resolved
            .with_metadata(ResolvedBatchMetadata { occurences: vec![occ.freeze()] });

        let (mut tx, rx) = connector_mut::<AlignmentSchema, ()>(1);
        stage.process(input, &mut tx).unwrap();
        drop(tx);

        let mut out_lens    = Vec::new();
        let mut out_offsets = Vec::new();
        let mut out_occs    = Vec::new();

        while let Ok(batch) = rx.recv() {
            let (lens, offsets, occs) = batch.columns((
                aligned::schema::resolved_len,
                aligned::schema::offset,
                aligned::schema::occurence,
            ));

            for i in 0..batch.len() {
                out_lens.push(lens[i]);
                out_offsets.push(offsets[i]);
                out_occs.push(occs[i]);
            }
        }

        assert_eq!(out_lens.len(), 3);
        assert_eq!(out_lens[0], 25); assert_eq!(out_offsets[0], 5); assert_eq!(out_occs[0], 100);
        assert_eq!(out_lens[1], 20); assert_eq!(out_offsets[1], 0); assert_eq!(out_occs[1], 200);
        assert_eq!(out_lens[2], 25); assert_eq!(out_offsets[2], 5); assert_eq!(out_occs[2], 300);
    }

    /// Verify that occurrences are correctly chunked into multiple output batches
    /// when their count exceeds the output pool's rows-per-batch capacity.
    ///
    /// 5 occurrences with capacity=3 must produce 2 batches (3 + 2).
    #[test]
    fn chunked_output_when_occurrences_exceed_batch_capacity() {

        let mut stage     = make_stage(3); // capacity = 3 forces chunking
        let resolved_pool = Pool::<ResolvedSchema>::new(4, 16);
        let occ_pool      = Pool::<SeqOccSchema>::new(4, 16);

        let mut resolved = resolved_pool.acquire().unwrap();
        resolved.set_len(1);
        resolved.mutate(
            (resolved::schema::seq_id, resolved::schema::resolved_len),
            |(seq_ids, lens)| { seq_ids[0] = 0; lens[0] = 10; },
        );

        let mut occ = occ_pool.acquire().unwrap();
        occ.set_len(5);
        occ.mutate(
            (occurences::schema::seq_id, occurences::schema::occurence),
            |(seq_ids, occs)| {
                for i in 0..5 { seq_ids[i] = 0; occs[i] = i as u64; }
            },
        );

        let input = resolved
            .with_metadata(ResolvedBatchMetadata { occurences: vec![occ.freeze()] });

        let (mut tx, rx) = connector_mut::<AlignmentSchema, ()>(16);
        stage.process(input, &mut tx).unwrap();
        drop(tx);

        let batches: Vec<_> = std::iter::from_fn(|| rx.recv().ok()).collect();
        let total: usize = batches.iter().map(|b| b.len()).sum();

        assert_eq!(total, 5);
        assert!(batches.len() == 2, "expected 2 batches, got {}", batches.len());
    }

    /// Verify that all occurrence batches attached to the resolved metadata are processed,
    /// and that rows from each batch are independently mapped through the same resolved index.
    ///
    /// Resolved:   seq_id=0 → resolved_len=20, offset=0
    ///             seq_id=1 → resolved_len=25, offset=5
    ///
    /// Occ batch A: [(seq_id=0, occ=10), (seq_id=1, occ=20)]
    /// Occ batch B: [(seq_id=1, occ=30), (seq_id=0, occ=40), (seq_id=0, occ=50)]
    ///
    /// Expected 5 output rows in order:
    ///   0: len=20, off=0, occ=10   (from A, seq 0)
    ///   1: len=25, off=5, occ=20   (from A, seq 1)
    ///   2: len=25, off=5, occ=30   (from B, seq 1)
    ///   3: len=20, off=0, occ=40   (from B, seq 0)
    ///   4: len=20, off=0, occ=50   (from B, seq 0)
    #[test]
    fn multiple_occurrence_batches_are_all_processed() {

        let mut stage     = make_stage(16);
        let resolved_pool = Pool::<ResolvedSchema>::new(4, 16);
        let occ_pool      = Pool::<SeqOccSchema>::new(4, 16);

        let mut resolved = resolved_pool.acquire().unwrap();
        resolved.set_len(2);
        resolved.mutate(
            (resolved::schema::seq_id, resolved::schema::resolved_len, resolved::schema::offset),
            |(seq_ids, lens, offsets)| {
                seq_ids[0] = 0; lens[0] = 20; offsets[0] = 0;
                seq_ids[1] = 1; lens[1] = 25; offsets[1] = 5;
            },
        );

        let mut occ_a = occ_pool.acquire().unwrap();
        occ_a.set_len(2);
        occ_a.mutate(
            (occurences::schema::seq_id, occurences::schema::occurence),
            |(seq_ids, occs)| {
                seq_ids[0] = 0; occs[0] = 10;
                seq_ids[1] = 1; occs[1] = 20;
            },
        );

        let mut occ_b = occ_pool.acquire().unwrap();
        occ_b.set_len(3);
        occ_b.mutate(
            (occurences::schema::seq_id, occurences::schema::occurence),
            |(seq_ids, occs)| {
                seq_ids[0] = 1; occs[0] = 30;
                seq_ids[1] = 0; occs[1] = 40;
                seq_ids[2] = 0; occs[2] = 50;
            },
        );

        let input = resolved
            .with_metadata(ResolvedBatchMetadata { 
                occurences: vec![occ_a.freeze(), occ_b.freeze()] 
            });

        let (mut tx, rx) = connector_mut::<AlignmentSchema, ()>(16);
        stage.process(input, &mut tx).unwrap();
        drop(tx);

        let mut out_lens    = Vec::new();
        let mut out_offsets = Vec::new();
        let mut out_occs    = Vec::new();

        while let Ok(batch) = rx.recv() {
            let (lens, offsets, occs) = batch.columns((
                aligned::schema::resolved_len,
                aligned::schema::offset,
                aligned::schema::occurence,
            ));
            
            for i in 0..batch.len() {
                out_lens.push(lens[i]);
                out_offsets.push(offsets[i]);
                out_occs.push(occs[i]);
            }
        }

        assert_eq!(out_lens.len(), 5);
        assert_eq!(out_lens[0], 20); assert_eq!(out_offsets[0], 0); assert_eq!(out_occs[0], 10);
        assert_eq!(out_lens[1], 25); assert_eq!(out_offsets[1], 5); assert_eq!(out_occs[1], 20);
        assert_eq!(out_lens[2], 25); assert_eq!(out_offsets[2], 5); assert_eq!(out_occs[2], 30);
        assert_eq!(out_lens[3], 20); assert_eq!(out_offsets[3], 0); assert_eq!(out_occs[3], 40);
        assert_eq!(out_lens[4], 20); assert_eq!(out_offsets[4], 0); assert_eq!(out_occs[4], 50);
    }
}
