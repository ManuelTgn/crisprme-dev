use crate::model::{
    alignment::{SeqMinedBatch, SeqResolvedBatch, SeqResolvedFrame},
    cigarx::{Cigarx, CigarxOp},
};
use columnar::{
    pipeline::{Emit, Stage, StageError},
    MemoryPool,
};
use itertools::izip;

/// Resolve mined alignments using the present cigarx
pub struct Resolver {
    pool: MemoryPool,
}

impl Resolver {
    pub fn new(pool: &MemoryPool) -> Self {
        Self { pool: pool.clone() }
    }
}

impl Stage for Resolver {
    type I = SeqMinedBatch;
    type O = SeqResolvedBatch;

    fn name() -> &'static str {
        "Resolver"
    }

    #[tracing::instrument(name = "pipeline:resolver", skip_all)]
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError> {
        let guide = input.guide.as_slice();

        // mined --1:1--> resolved
        input.sequences.with_cols(|sequences| {
            input.mined.with_cols(|mut mined| {
                let source_seq_count = sequences.content.rows();
                let rows = mined.seq_row_idx.rows();

                tracing::info!("received {} rows to resolve", rows);

                let mut resolved = SeqResolvedFrame::empty();
                resolved.with_cols(|mut resolved| {
                    // Share columns (seq_id, offset)
                    resolved.seq_row_idx.shared(&mut mined.seq_row_idx);
                    resolved.offset.shared(&mut mined.offset);

                    // Allocate columns (rguide, rseq)
                    resolved.rguide.alloc(&self.pool, rows);
                    resolved.rseq.alloc(&self.pool, rows);

                    // Zipped iterator over all used columns
                    let zipper = izip!(
                        resolved.seq_row_idx.iter(),
                        resolved.rguide.iter_mut(),
                        resolved.rseq.iter_mut(),
                        mined.cigarx.iter(),
                        mined.offset.iter()
                    );

                    // Resolve the guide and sequence
                    for (seq_row_idx, rguide, rseq, cigarx, offset) in zipper {
                        // Indirect look-up to sequence content
                        // NOTE: it should be fast enough
                        let sequence = sequences.content.get(*seq_row_idx as usize);

                        let mut gpos = 0usize;
                        let mut spos = *offset as usize; // start at alignment position in sequence
                        let mut opos = 0usize;

                        for op in cigarx.iter() {
                            match op {
                                CigarxOp::Match | CigarxOp::Mismatch => {
                                    rguide[opos] = guide[gpos].to_ascii();
                                    rseq[opos] = sequence[spos].to_ascii();
                                    gpos += 1;
                                    spos += 1;
                                }
                                CigarxOp::Deletion => {
                                    rguide[opos] = b'-';
                                    rseq[opos] = sequence[spos].to_ascii();
                                    spos += 1;
                                }
                                CigarxOp::Insertion => {
                                    rguide[opos] = guide[gpos].to_ascii();
                                    rseq[opos] = b'-';
                                    gpos += 1;
                                }
                            }
                            opos += 1;
                        }

                        // Null-terminate both resolved arrays
                        rguide[opos] = 0;
                        rseq[opos] = 0;
                    }
                });

                emitter.emit(SeqResolvedBatch {
                    source_seq_count,
                    occurences: input.occurences,
                    resolved,
                })
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;

    use columnar::pipeline::Stage;

    use crate::{
        crispr::guide::Guide,
        model::{
            alignment::{SeqMinedBatch, SeqMinedFrame, SeqResolvedBatch},
            cigarx::{Cigarx, Cigarx64, CigarxOp},
            input::{SeqFrame, SeqOccFrame},
            occurence::Occurence,
        },
        pipeline::test::{make_pool, Collector},
        sequence::iupac::Iupac,
    };

    use super::*;

    // Build a SeqFrame from ASCII strings (one per sequence, padded with zeros)
    fn make_seqs(pool: &MemoryPool, seqs: &[&str]) -> SeqFrame {
        let mut frame = SeqFrame::alloc(pool, seqs.len());
        frame.with_cols(|mut cols| {
            for (i, s) in seqs.iter().enumerate() {
                let content = cols.content.get_mut(i);
                for (j, b) in s.bytes().enumerate() {
                    content[j] = Iupac::from_ascii_lossy(b);
                }
            }
        });
        frame
    }

    // Build a SeqMinedFrame from (seq_row_idx, cigarx, offset) tuples
    fn make_mined(pool: &MemoryPool, rows: &[(u32, Cigarx64, u8)]) -> SeqMinedFrame {
        let mut frame = SeqMinedFrame::alloc(pool, rows.len());
        frame.with_cols(|mut cols| {
            for (i, (idx, cigarx, offset)) in rows.iter().enumerate() {
                *cols.seq_row_idx.get_mut(i) = *idx;
                *cols.cigarx.get_mut(i) = *cigarx;
                *cols.offset.get_mut(i) = *offset;
            }
        });
        frame
    }

    // Build a single Cigarx64 from a slice of ops
    fn cigar(ops: &[CigarxOp]) -> Cigarx64 {
        let mut c = Cigarx64::default();
        for &op in ops {
            c.push(op);
        }
        c
    }

    fn make_batch(
        pool: &MemoryPool,
        guide: &str,
        seqs: &[&str],
        rows: &[(u32, Cigarx64, u8)],
    ) -> SeqMinedBatch {
        SeqMinedBatch {
            guide: Guide::new(guide),
            sequences: make_seqs(pool, seqs),
            occurences: SeqOccFrame::alloc(pool, 0), // not used by resolver
            mined: make_mined(pool, rows),
        }
    }

    // All-match cigar: both rguide and rseq get the bases as-is
    #[test]
    fn all_match_copies_guide_and_sequence() {
        let pool = make_pool();
        let collector = Collector(RefCell::new(vec![]));

        let mut resolver = Resolver::new(&pool);

        let ops = cigar(&[CigarxOp::Match; 4]);
        resolver
            .process(
                make_batch(&pool, "ACGT", &["ACGT"], &[(0, ops, 0)]),
                &collector,
            )
            .unwrap();

        let mut outputs = collector.into_inner();
        outputs[0].resolved.with_cols(|cols| {
            let rguide = cols.rguide.get(0);
            let rseq = cols.rseq.get(0);
            assert_eq!(&rguide[..5], b"ACGT\0");
            assert_eq!(&rseq[..5], b"ACGT\0");
        });
    }

    // Deletion op: guide gets '-', sequence base is copied
    #[test]
    fn deletion_inserts_gap_in_guide() {
        let pool = make_pool();
        let collector = Collector(RefCell::new(vec![]));

        let mut resolver = Resolver::new(&pool);

        // guide "ACT" (3 bases), sequence "ACGT" — cigar M M D M
        // pos 0: M → guide[0]=A, seq[0]=A
        // pos 1: M → guide[1]=C, seq[1]=C
        // pos 2: D → guide='-',  seq[2]=G
        // pos 3: M → guide[2]=T, seq[3]=T
        let ops = cigar(&[
            CigarxOp::Match,
            CigarxOp::Match,
            CigarxOp::Deletion,
            CigarxOp::Match,
        ]);
        resolver
            .process(
                make_batch(&pool, "ACT", &["ACGT"], &[(0, ops, 0)]),
                &collector,
            )
            .unwrap();

        let mut outputs = collector.into_inner();
        outputs[0].resolved.with_cols(|cols| {
            assert_eq!(&cols.rguide.get(0)[..5], b"AC-T\0");
            assert_eq!(&cols.rseq.get(0)[..5], b"ACGT\0");
        });
    }

    // Insertion op: sequence gets '-', guide base is copied
    #[test]
    fn insertion_inserts_gap_in_sequence() {
        let pool = make_pool();
        let collector = Collector::new();

        let mut resolver = Resolver::new(&pool);

        // guide "ACGT" (4 bases), sequence "AGT" — cigar M I M M
        // pos 0: M → guide[0]=A, seq[0]=A
        // pos 1: I → guide[1]=C, seq='-'
        // pos 2: M → guide[2]=G, seq[1]=G
        // pos 3: M → guide[3]=T, seq[2]=T
        let ops = cigar(&[
            CigarxOp::Match,
            CigarxOp::Insertion,
            CigarxOp::Match,
            CigarxOp::Match,
        ]);
        resolver
            .process(
                make_batch(&pool, "ACGT", &["AGT"], &[(0, ops, 0)]),
                &collector,
            )
            .unwrap();

        let mut outputs = collector.into_inner();
        outputs[0].resolved.with_cols(|cols| {
            assert_eq!(&cols.rguide.get(0)[..5], b"ACGT\0");
            assert_eq!(&cols.rseq.get(0)[..5], b"A-GT\0");
        });
    }

    // Offset shifts the start position in the sequence
    #[test]
    fn offset_shifts_sequence_start() {
        let pool = make_pool();
        let collector = Collector::new();

        let mut resolver = Resolver::new(&pool);

        // sequence "TTACGT", guide "AC", offset=2 → aligned region is "ACGT"[2..4]
        let ops = cigar(&[CigarxOp::Match, CigarxOp::Match]);
        resolver
            .process(
                make_batch(&pool, "AC", &["TTACGT"], &[(0, ops, 2)]),
                &collector,
            )
            .unwrap();

        let mut outputs = collector.into_inner();
        outputs[0].resolved.with_cols(|cols| {
            assert_eq!(&cols.rguide.get(0)[..3], b"AC\0");
            assert_eq!(&cols.rseq.get(0)[..3], b"AC\0");
        });
    }

    // Multiple mined rows referencing different sequences are resolved correctly
    #[test]
    fn multiple_rows_look_up_correct_sequence() {
        let pool = make_pool();
        let collector = Collector::new();

        let mut resolver = Resolver::new(&pool);

        // seq 0 = "AAAA", seq 1 = "CCCC", guide = "AC"
        // row 0: seq_row_idx=1, cigar=[M,M] → rseq = "CC"
        // row 1: seq_row_idx=0, cigar=[M,M] → rseq = "AA"
        let ops = cigar(&[CigarxOp::Match, CigarxOp::Match]);
        resolver
            .process(
                make_batch(&pool, "AC", &["AAAA", "CCCC"], &[(1, ops, 0), (0, ops, 0)]),
                &collector,
            )
            .unwrap();

        let mut outputs = collector.into_inner();
        outputs[0].resolved.with_cols(|cols| {
            assert_eq!(&cols.rseq.get(0)[..3], b"CC\0"); // seq 1
            assert_eq!(&cols.rseq.get(1)[..3], b"AA\0"); // seq 0
        });
    }

    // source_seq_count equals the number of sequences in the SeqFrame
    #[test]
    fn source_seq_count_is_sequence_frame_len() {
        let pool = make_pool();
        let collector = Collector::new();

        let mut resolver = Resolver::new(&pool);

        let ops = cigar(&[CigarxOp::Match]);
        resolver
            .process(
                make_batch(&pool, "A", &["A", "C", "G"], &[(0, ops, 0)]),
                &collector,
            )
            .unwrap();

        let outputs = collector.into_inner();
        assert_eq!(outputs[0].source_seq_count, 3);
    }

    // Occurences are passed through unchanged
    #[test]
    fn occurences_are_passed_through() {
        let pool = make_pool();
        let collector = Collector::new();

        let mut resolver = Resolver::new(&pool);

        let mut occs = SeqOccFrame::alloc(&pool, 2);
        occs.with_cols(|mut cols| {
            *cols.seq_row_idx.get_mut(0) = 0;
            *cols.occurence.get_mut(0) = Occurence::new(1, 100, 0);
            *cols.seq_row_idx.get_mut(1) = 0;
            *cols.occurence.get_mut(1) = Occurence::new(2, 200, 1);
        });

        let ops = cigar(&[CigarxOp::Match]);
        let batch = SeqMinedBatch {
            guide: Guide::new("A"),
            sequences: make_seqs(&pool, &["A"]),
            occurences: occs,
            mined: make_mined(&pool, &[(0, ops, 0)]),
        };

        resolver.process(batch, &collector).unwrap();

        let mut outputs = collector.into_inner();
        outputs[0].occurences.with_cols(|cols| {
            assert_eq!(cols.seq_row_idx.rows(), 2);
            assert_eq!(*cols.seq_row_idx.get(0), 0);
            assert_eq!(*cols.seq_row_idx.get(1), 0);
        });
    }
}
