
use columnar::{
    Column, MemoryPool, pipeline::{Emit, Stage, StageError}
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
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError> {

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

                            let (dst_id, dst_occ, dst_offset, dst_rguide, dst_rseq) = dst_iter.next().unwrap();

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
