
use columnar::{
    Column, MemoryPool, pipeline::{Emit, Stage, PipelineError}
};
use itertools::izip;

use crate::model::alignment::{AlignmentFrame, SeqResolvedBatch};

/// Cross-product scatter: each occurrence is paired with every resolved entry sharing its
/// seq_row_idx. Occurrences with no resolved match produce no output (inner join).
///
/// Invariants:
///   - every resolved.seq_row_idx has at least one matching occurrence (not vice versa)
///   - all seq_row_idx values < source_seq_count
/// 
pub struct Broadcast {

    /// Memory pool used to allocate outputs
    pool: MemoryPool,

    /// table[seq_row_idx] = (end, count)
    ///   - count: number of resolved rows for this seq_row_idx
    ///   - end:   one-past-the-last index into flat_rows (i.e. start = end - count)
    table: Vec<(usize, usize)>,

    /// Flat buffer of resolved row indices, grouped by seq_row_idx.
    /// Reused across batches; resized to total resolved count each batch.
    flat_rows: Vec<usize>,

    /// Tracks which table slots were written in the current batch for selective reset.
    /// Kept on the struct to avoid per-batch allocation.
    written: Vec<usize>,
}

impl Broadcast {

    pub fn new(pool: &MemoryPool) -> Self {
        Self { pool: pool.clone(), table: vec![], flat_rows: vec![], written: vec![] }
    }

    /// Populates `table` and `flat_rows` from the resolved `seq_row_idx` column.
    /// Two passes: count per slot, then fill flat_rows using table[slot].0 as a cursor.
    /// After this call: table[slot] = (end, count), range = flat_rows[end-count..end].
    fn build(&mut self, seq_row_idx: &Column<'_, u32>) {

        self.written.clear();

        for idx in seq_row_idx.iter() {
            let slot = *idx as usize;

            // Mark this slot as dirty
            if self.table[slot].1 == 0 { 
                self.written.push(slot); 
            }

            self.table[slot].1 += 1;
        }

        tracing::debug!("source sequence for each alignment (0..A): {:?}",
            seq_row_idx.iter().collect::<Vec<_>>());

        tracing::debug!("alignments count for each source sequence (0..S): {:?}", 
            self.table.iter().map(|(_, count)| count).collect::<Vec<_>>());

        self.flat_rows.resize(seq_row_idx.rows(), 0);

        let mut offset = 0usize;
        for &slot in &self.written {
            self.table[slot].0 = offset;
            offset += self.table[slot].1;
        }

        for (i, idx) in seq_row_idx.iter().enumerate() {
            let slot = *idx as usize;
            self.flat_rows[self.table[slot].0] = i;
            self.table[slot].0 += 1;
        }
    }

    /// Returns the number of output rows: each occurrence contributes
    /// as many rows as there are resolved entries for its seq_row_idx.
    fn count(&self, seq_row_idx: &Column<'_, u32>) -> usize {

        tracing::debug!("alignments count for each occurence (0..C): {:?}",
            seq_row_idx.iter()
                .map(|idx| self.table[*idx as usize].1)
                .collect::<Vec<_>>()
        );

        seq_row_idx.iter()
            .map(|idx| self.table[*idx as usize].1)
            .sum()
    }

    /// Resets only the slots written in the current batch back to (0, 0).
    fn reset(&mut self) {
        for &slot in &self.written {
            self.table[slot] = (0, 0);
        }
    }
}

impl Stage for Broadcast {

    type I = SeqResolvedBatch;
    type O = AlignmentFrame;

    fn name() -> &'static str { "Broadcast" }

    #[tracing::instrument(name = "pipeline:broadcast", skip_all)]
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), PipelineError> {

        if self.table.len() < input.source_seq_count {
            self.table.resize(input.source_seq_count, (0, 0));
        }

        input.resolved.with_cols(|resolved| {
            input.occurences.with_cols(|occurence| {

                tracing::info!("starting broadcast of {} resolved rows into {} occurences",
                    resolved.seq_row_idx.rows(), occurence.seq_row_idx.rows());

                self.build(&resolved.seq_row_idx);

                let out_rows = self.count(&occurence.seq_row_idx);
                let mut alignment = AlignmentFrame::alloc(&self.pool, out_rows);
                alignment.with_cols(|mut alignment| {

                    let src_iter = izip!(
                        occurence.seq_row_idx.iter(),
                        occurence.occurence.iter(),
                    );
                    
                    let mut dst_iter = izip!(
                        alignment.seq_row_idx.iter_mut(),
                        alignment.occurence.iter_mut(),
                        alignment.offset.iter_mut(),
                        alignment.rguide.iter_mut(),
                        alignment.rseq.iter_mut(),
                    );

                    for (seq_idx, src_occ) in src_iter {
                        
                        let (end, count) = self.table[*seq_idx as usize];
                        for &row in &self.flat_rows[end - count..end] {

                            let (
                                dst_id, 
                                dst_occ, 
                                dst_offset, 
                                dst_rguide, 
                                dst_rseq
                            ) = dst_iter.next().unwrap();

                            *dst_id     = *seq_idx;
                            *dst_occ    = *src_occ;
                            *dst_offset = *resolved.offset.get(row);
                            *dst_rguide = *resolved.rguide.get(row);
                            *dst_rseq   = *resolved.rseq.get(row);
                        }
                    }
                });

                tracing::info!("total broadcasted rows: {}", out_rows);

                self.reset();
                emitter.emit(alignment)
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;

    use columnar::pipeline::Stage;

    use crate::{
        model::{
            alignment::{AlignmentFrame, SeqResolvedBatch, SeqResolvedFrame},
            input::SeqOccFrame,
            occurence::Occurence,
        },
        pipeline::test::{Collector, make_pool},
    };

    use super::*;

    // Convenience: make a [u8; 32] with `b` in the first byte
    fn arr(b: u8) -> [u8; 32] {
        let mut a = [0u8; 32];
        a[0] = b;
        a
    }

    fn make_occs(pool: &MemoryPool, pairs: &[(u32, u64)]) -> SeqOccFrame {
        let mut frame = SeqOccFrame::alloc(pool, pairs.len());
        frame.with_cols(|mut cols| {
            for (i, (idx, occ)) in pairs.iter().enumerate() {
                *cols.seq_row_idx.get_mut(i) = *idx;
                *cols.occurence.get_mut(i) = Occurence(*occ);
            }
        });
        frame
    }

    fn make_resolved(pool: &MemoryPool, rows: &[(u32, u8, [u8; 32], [u8; 32])]) -> SeqResolvedFrame {
        let mut frame = SeqResolvedFrame::alloc(pool, rows.len());
        frame.with_cols(|mut cols| {
            for (i, (idx, offset, rguide, rseq)) in rows.iter().enumerate() {
                *cols.seq_row_idx.get_mut(i) = *idx;
                *cols.offset.get_mut(i) = *offset;
                *cols.rguide.get_mut(i) = *rguide;
                *cols.rseq.get_mut(i) = *rseq;
            }
        });
        frame
    }

    // 1 resolved × 1 occurrence = 1 output row
    #[test]
    fn single_resolved_single_occurrence() {
        let pool = make_pool();
        let collector = Collector(RefCell::new(vec![]));
        let mut broadcast = Broadcast::new(&pool);

        broadcast.process(SeqResolvedBatch {
            source_seq_count: 1,
            occurences: make_occs(&pool, &[(0, 100)]),
            resolved: make_resolved(&pool, &[(0, 42, arr(1), arr(2))]),
        }, &collector).unwrap();

        let mut outputs = collector.into_inner();
        assert_eq!(outputs.len(), 1);
        outputs[0].with_cols(|cols| {
            assert_eq!(cols.seq_row_idx.rows(), 1);
            assert_eq!(*cols.offset.get(0), 42);
            assert_eq!(cols.rguide.get(0)[0], 1);
            assert_eq!(cols.rseq.get(0)[0], 2);
        });
    }

    // 2 resolved × 2 occurrences = 4 output rows (full cross product)
    #[test]
    fn cross_product_two_resolved_two_occurrences() {
        let pool = make_pool();
        let collector = Collector(RefCell::new(vec![]));
        let mut broadcast = Broadcast::new(&pool);

        broadcast.process(SeqResolvedBatch {
            source_seq_count: 1,
            occurences: make_occs(&pool, &[(0, 100), (0, 200)]),
            resolved: make_resolved(&pool, &[
                (0, 10, arr(10), arr(11)),
                (0, 20, arr(20), arr(21)),
            ]),
        }, &collector).unwrap();

        // Scatter order: for each occurrence, for each resolved row
        // occ[0]×res[0], occ[0]×res[1], occ[1]×res[0], occ[1]×res[1]
        let mut outputs = collector.into_inner();
        assert_eq!(outputs.len(), 1);
        outputs[0].with_cols(|cols| {
            assert_eq!(cols.seq_row_idx.rows(), 4);
            assert_eq!(*cols.offset.get(0), 10);
            assert_eq!(*cols.offset.get(1), 20);
            assert_eq!(*cols.offset.get(2), 10);
            assert_eq!(*cols.offset.get(3), 20);
        });
    }

    // Occurrences with no resolved match are dropped (inner join)
    #[test]
    fn unmatched_occurrences_are_dropped() {
        let pool = make_pool();
        let collector = Collector(RefCell::new(vec![]));
        let mut broadcast = Broadcast::new(&pool);

        // seq 0 has resolved; seq 1 does not
        broadcast.process(SeqResolvedBatch {
            source_seq_count: 2,
            occurences: make_occs(&pool, &[(0, 100), (0, 200), (1, 300), (1, 400)]),
            resolved: make_resolved(&pool, &[(0, 42, arr(1), arr(2))]),
        }, &collector).unwrap();

        let mut outputs = collector.into_inner();
        assert_eq!(outputs.len(), 1);
        outputs[0].with_cols(|cols| {
            assert_eq!(cols.seq_row_idx.rows(), 2);
            assert!(cols.seq_row_idx.iter().all(|idx| *idx == 0));
        });
    }

    // resolved.seq_row_idx can arrive in any order
    #[test]
    fn random_resolved_order_matches_correctly() {
        let pool = make_pool();
        let collector = Collector(RefCell::new(vec![]));
        let mut broadcast = Broadcast::new(&pool);

        // Resolved arrives in reverse order: seq 2, 0, 1
        broadcast.process(SeqResolvedBatch {
            source_seq_count: 3,
            occurences: make_occs(&pool, &[(0, 100), (1, 200), (2, 300)]),
            resolved: make_resolved(&pool, &[
                (2, 30, arr(30), arr(31)),
                (0, 10, arr(10), arr(11)),
                (1, 20, arr(20), arr(21)),
            ]),
        }, &collector).unwrap();

        let mut outputs = collector.into_inner();
        outputs[0].with_cols(|cols| {
            assert_eq!(cols.seq_row_idx.rows(), 3);
            assert_eq!(*cols.offset.get(0), 10); // occ seq 0 → resolved offset 10
            assert_eq!(*cols.offset.get(1), 20); // occ seq 1 → resolved offset 20
            assert_eq!(*cols.offset.get(2), 30); // occ seq 2 → resolved offset 30
        });
    }

    // Table must be fully reset between batches; stale entries from batch N must not
    // contaminate batch N+1
    #[test]
    fn table_is_reset_between_batches() {
        let pool = make_pool();
        let collector = Collector(RefCell::new(vec![]));
        let mut broadcast = Broadcast::new(&pool);

        // Batch 1: seq 0 has a resolved entry
        broadcast.process(SeqResolvedBatch {
            source_seq_count: 2,
            occurences: make_occs(&pool, &[(0, 100)]),
            resolved: make_resolved(&pool, &[(0, 10, arr(1), arr(2))]),
        }, &collector).unwrap();

        // Batch 2: only seq 1 has a resolved entry; seq 0 does NOT
        // If the table were not reset, seq 0's occurrence would pick up the stale entry
        broadcast.process(SeqResolvedBatch {
            source_seq_count: 2,
            occurences: make_occs(&pool, &[(0, 200), (1, 300)]),
            resolved: make_resolved(&pool, &[(1, 20, arr(3), arr(4))]),
        }, &collector).unwrap();

        let mut outputs = collector.into_inner();
        outputs[0].with_cols(|cols| assert_eq!(cols.seq_row_idx.rows(), 1)); // batch 1
        outputs[1].with_cols(|cols| {
            assert_eq!(cols.seq_row_idx.rows(), 1); // only seq 1 survives
            assert_eq!(*cols.seq_row_idx.get(0), 1);
        });
    }
}