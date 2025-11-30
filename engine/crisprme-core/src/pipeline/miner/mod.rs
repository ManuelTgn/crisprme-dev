use crate::{
    common::{
        alignment::{Alignment, AlignmentOp, AlignmentState},
        cigarx::{Cigarx, CigarxOp},
        guide::Guide,
        iupac::Iupac,
        thresholds::Thresholds,
    },
    memory::arena::Memory,
};

/*
pub struct SimpleMiner<'mem> {
    batch: &'mem SequenceBatch<'mem>,
    guide: Guide,
    thresholds: Thresholds,

    /// Cache for beamlight search, allocated only once
    cache_min_ed: &'mem mut [u32],
    cache_seq_idx: usize,

    stack: Vec<AlignmentState>,
    cigarx: Cigarx<u64>,

    stored_seq_idx: usize,
    stored_offset: usize,
}

impl<'mem> SimpleMiner<'mem> {
    pub fn new(
        mem: &'mem Memory,
        batch: &'mem SequenceBatch<'mem>,
        guide: &Guide,
        thresholds: &Thresholds,
    ) -> Self {
        let cache = mem.alloc_slice_fill((guide.len() + 1) * (batch.seq_len() + 1), 0);
        SimpleMiner {
            batch,
            thresholds: *thresholds,
            guide: guide.clone(),
            stack: Vec::with_capacity(64),
            cigarx: Cigarx::default(),
            cache_min_ed: cache.into_mut(),
            cache_seq_idx: usize::MAX,
            stored_seq_idx: 0,
            stored_offset: 0,
        }
    }
}

impl<'mem> SimpleMiner<'mem> {
    fn invalid_state(&self, state: AlignmentState) -> bool {
        !state.below_thresholds(&self.thresholds)
            || state.tidx() >= self.batch.seq_len() as u32
            || state.qidx() >= self.guide.len() as u32
    }

    fn push(&mut self, op: CigarxOp, state: AlignmentState) {
        let next_state = match op {
            CigarxOp::Match => state.op_match(),
            CigarxOp::Mismatch => state.op_mismatch(),
            CigarxOp::Deletion => state.op_delete(),
            CigarxOp::Insertion => state.op_insert(),
        };

        //println!("\tnext_state: {next_state:?}");
        self.stack.push(next_state);
        self.cigarx.push(op);
    }

    /// Recompute the minimum edit distance cache
    fn fill_min_ed(&mut self, q: &[Iupac], t: &[Iupac]) {
        // Solution, all of the bottom row is a solution
        for ti in 0..=self.batch.seq_len() {
            *self.min_ed_mut(self.guide.len(), ti) = 0;
        }

        // Insert remaining T as gaps
        for qi in 0..=self.guide.len() {
            let value = (self.guide.len() - qi) as u32;
            *self.min_ed_mut(qi, self.batch.seq_len()) = value;
        }

        for qi in (0..self.guide.len()).rev() {
            for ti in (0..self.batch.seq_len()).rev() {
                // NOTE: This also handles wildcards
                let m = if q[qi].matches(t[ti]) { 0 } else { 1 };

                let a = self.min_ed(qi + 1, ti + 1) + m;
                let b = self.min_ed(qi + 1, ti + 0) + 1;
                let c = self.min_ed(qi + 0, ti + 1) + 1;
                let value = a.min(b).min(c);

                *self.min_ed_mut(qi, ti) = value;
            }
        }

        /*
        // Print matrix
        for qi in 0..=self.config.query_len {
            println!();
            for ti in 0..=self.config.target_len {
                print!("{:>3}", self.min_ed(qi, ti));
            }
        }
        */
    }

    fn min_ed(&self, q: usize, t: usize) -> u32 {
        self.cache_min_ed[q * (self.batch.seq_len() + 1) + t]
    }

    fn min_ed_mut(&mut self, q: usize, t: usize) -> &mut u32 {
        &mut self.cache_min_ed[q * (self.batch.seq_len() + 1) + t]
    }
}

/// Generate all alignments
impl<'mem> Iterator for SimpleMiner<'mem> {
    type Item = (usize, Alignment);
    fn next(&mut self) -> Option<Self::Item> {

        /*
        // Create cache for current sequence if necesary
        if self.cache_seq_idx != self.stored_seq_idx {
            let s = self.batch.sequence_at(self.stored_seq_idx);
            let g = self.guide.clone();
            self.fill_min_ed(&g, &s);
        }
        */

        // Process the alignments from the last sequence
        for curr_seq_idx in self.stored_seq_idx..self.batch.len() {
            // Process from the last offset
            for curr_offset in self.stored_offset..self.batch.seq_len() {
                //println!("section, seq-idx: {curr_seq_idx}, offset: {curr_offset}");

                // Add state if stack is empty, this is always true if we are at the beginning
                // of a new pair of sequence idx and offset
                if self.stack.is_empty() {
                    //println!("stack empty, adding initial state at {curr_offset}");
                    self.stack.push(AlignmentState::initial(curr_offset as u32));
                    self.cigarx = Cigarx::default();
                }

                // Get state at the top of the stack
                while let Some(state) = self.stack.last() {
                    let state = *state;

                    //println!("current state [{}]: {state:?}", self.stack.len() - 1);
                    //println!("current cigarx: {:?}", self.cigarx);

                    // Return solution
                    if state.qidx() as usize == self.guide.len() && !state.invalid(&self.thresholds)
                    {
                        let result_cigarx = self.cigarx;

                        // Remove solution state, it doesn't make sense to continue
                        // as you can only add tgaps at the end of the alignment in our use-case
                        self.cigarx.pop();
                        self.stack.pop();

                        self.stored_seq_idx = curr_seq_idx;
                        self.stored_offset = curr_offset;

                        // Generate solution
                        let solution = Alignment::new(
                            self.batch.id_at(curr_seq_idx),
                            result_cigarx,
                            curr_offset as u8,
                        );

                        //println!("\tsolution! {solution:?}");
                        return Some((curr_seq_idx, solution));
                    }

                    // Remove invalid states
                    if self.invalid_state(state) {
                        //println!("\tinvalid!");
                        self.cigarx.pop();
                        self.stack.pop();
                        continue;
                    }

                    /*
                    // Prune if minimum edit distance to solution is incompatible
                    if state.ed() + self.min_ed(state.qidx() as usize, state.tidx() as usize)
                        > self.thresholds.ed()
                    {
                        println!("\tpruned!");
                        self.cigarx.pop();
                        self.stack.pop();
                        continue;
                    }
                    */

                    // Continue exploration
                    let (travel_state, travel) = state.travel();
                    let last_idx = self.stack.len() - 1;
                    self.stack[last_idx] = travel_state;

                    //println!("\ttravel: {travel:?}");

                    match travel {
                        AlignmentOp::Match => {
                            let t = self.batch.sequence_at(curr_seq_idx).as_slice()
                                [state.tidx() as usize];
                            let g = self.guide[state.qidx() as usize];

                            //println!("\tt:{t:?} vs g:{g:?} -> {}", g.matches(t));
                            match g.matches(t) {
                                false => self.push(CigarxOp::Mismatch, state),
                                true => self.push(CigarxOp::Match, state),
                            }
                        }
                        AlignmentOp::Deletion => self.push(CigarxOp::Deletion, state),
                        AlignmentOp::Insertion => self.push(CigarxOp::Insertion, state),
                        AlignmentOp::Exhausted => {
                            self.cigarx.pop();
                            self.stack.pop();
                        }
                    }

                    // Reset offset for the new sequence
                    self.stored_offset = 0;
                }
            }
        }

        // Nothing more to mine
        None
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use super::SimpleMiner;
    use crate::{
        common::{
            alignment::Alignment, cigarx::Cigarx, guide::Guide, iupac::Iupac,
            thresholds::Thresholds,
        },
        memory::{
            arena::Arena,
            batch::{AlignmentBatch, SequenceBatch},
        },
    };

    #[test]
    fn simple_miner_1() {
        let guide = Guide::from("AN");
        let ids = &mut [0, 1];
        let sequences = &mut [
            Iupac::from_utf8('A'),
            Iupac::from_utf8('N'),
            Iupac::from_utf8('T'),
            Iupac::from_utf8('A'),
            Iupac::from_utf8('N'),
            Iupac::from_utf8('T'),
        ];

        let correct = HashSet::from([
            Alignment::new(0, Cigarx::from("=="), 0),
            Alignment::new(0, Cigarx::from("=="), 1),
            Alignment::new(1, Cigarx::from("=="), 0),
            Alignment::new(1, Cigarx::from("=="), 1),
        ]);

        let mut arena = Arena::alloc(1024);
        arena.scoped(|memory| {
            let batch = SequenceBatch::new(sequences, ids, 3, false);
            let miner = SimpleMiner::new(
                &memory,
                &batch,
                &guide,
                &Thresholds {
                    qgap: 0,
                    tgap: 0,
                    mism: 1,
                },
            );

            let mut alignments = AlignmentBatch::new_in(&memory, 100, false);
            for (_, align) in miner {
                assert!(alignments.push(align));
            }

            //println!("{:?}", alignments);
            for align in alignments.alignments() {
                assert!(correct.contains(align));
            }
        });
    }

    #[test]
    fn simple_miner_2() {
        let guide = Guide::from("AT");
        let ids = &mut [0];
        let sequences = &mut [
            Iupac::from_utf8('A'),
            Iupac::from_utf8('C'),
            Iupac::from_utf8('T'),
        ];

        let correct = HashSet::from([
            Alignment::new(0, Cigarx::from("=="), 0),
            Alignment::new(0, Cigarx::from("=="), 1),
            Alignment::new(1, Cigarx::from("=="), 0),
            Alignment::new(1, Cigarx::from("=="), 1),
        ]);

        let mut arena = Arena::alloc(1024);
        arena.scoped(|memory| {
            let batch = SequenceBatch::new(sequences, ids, 3, false);
            let miner = SimpleMiner::new(
                &memory,
                &batch,
                &guide,
                &Thresholds {
                    qgap: 1,
                    tgap: 1,
                    mism: 0,
                },
            );

            let mut alignments = AlignmentBatch::new_in(&memory, 100, false);
            for (_, align) in miner {
                assert!(alignments.push(align));
            }

            println!("{:?}", alignments);
            for align in alignments.alignments() {
                assert!(correct.contains(align));
            }
        });
    }
}
*/
