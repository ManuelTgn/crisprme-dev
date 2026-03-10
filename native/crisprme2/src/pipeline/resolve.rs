use std::sync::Arc;

use columnar::{arena::Arena, pipeline::{Emit, Stage, StageError}, pool::{BatchMut, Pool}};

use crate::model::{alignment::{MinedBatchMetadata, MinedSchema, ResolvedBatchMetadata, ResolvedSchema}, cigarx::{Cigarx, CigarxOp}};
use crate::sequence::iupac::Iupac;

/// Resolve mined alignments using the present cigarx
pub struct AlignmentSimpleResolve {
    
    /// Pool for buffers of resolved alignment schema
    pool: Arc<Pool<ResolvedSchema>>,
    /// Temporary buffer
    arena: Arena,
}

impl Stage for AlignmentSimpleResolve {

    type Input  = BatchMut<MinedSchema, MinedBatchMetadata>;
    type Output = BatchMut<ResolvedSchema, ResolvedBatchMetadata>;

    fn process<E>(&mut self, input: Self::Input, emitter: &mut E) -> Result<(), StageError>
    where
        E: Emit<Self::Output>
    {
        use crate::model::alignment::mined::schema    as ms;
        use crate::model::alignment::resolved::schema as rs;
        use crate::model::input::sequences::schema    as ss;

        let source_seq_batch = &input.metadata.sequences;
        let guide = &source_seq_batch.metadata.guide;

        let (sequences,) = source_seq_batch.columns((ss::content,));
        let (mined_seq_ids, mined_offsets, mined_cigarxs) = 
                input.columns((ms::seq_id, ms::offset, ms::cigarx));

        self.arena.scoped(|m| {

            // Map from mined_seq_id to index inside the sequence batch
            // NOTE: seq_id is always from 0 to input.len()
            let mut index = m.alloc_slice_fill(input.len(), 0u32);
            for (i, id) in mined_seq_ids.iter().enumerate() { 
                index[*id as usize] = i as u32; 
            }

            let mut remaining = input.len();
            while remaining > 0 {

                // Acquire a new result batch
                let mut result = self.pool.acquire()
                    .map_err(|_| StageError)?;

                let rows = remaining.min(result.capacity());
                result.set_len(rows);
                result.mutate(
                    (rs::seq_id, rs::resolved_len, rs::offset, rs::rguide, rs::rseq), 
                    |(seq_ids, lens, offsets, rguides, rseqs)| {
                        for i in 0..rows {
                            
                            let cigarx = mined_cigarxs[i];
                            let offset = mined_offsets[i];

                            seq_ids[i] = mined_seq_ids[i]; 
                            lens[i]    = cigarx.len() as u8;
                            offsets[i] = offset;

                            // Indirect look-up to sequence content
                            let sequence = &sequences[index[i] as usize];
                            let guide = guide.as_slice();
                            
                            // We can fill rguide and rseq
                            let mut gpos = 0usize;
                            let mut spos = offset as usize;
                            let mut opos = 0usize;

                            for op in cigarx.iter() {
                                match op {
                                    CigarxOp::Match | CigarxOp::Mismatch => {
                                        rguides[i][opos] = guide[gpos];
                                        rseqs[i][opos]   = sequence[spos];
                                        gpos += 1;
                                        spos += 1;
                                    }
                                    CigarxOp::Deletion => {
                                        // Sequence has a base, guide has a gap
                                        rguides[i][opos] = Iupac::default();
                                        rseqs[i][opos]   = sequence[spos];
                                        spos += 1;
                                    }
                                    CigarxOp::Insertion => {
                                        // Guide has a base, sequence has a gap
                                        rguides[i][opos] = guide[gpos];
                                        rseqs[i][opos]   = Iupac::default();
                                        gpos += 1;
                                    }
                                }
                                opos += 1;
                            }
                        }
                    });

                remaining -= result.len();
                emitter.emit(result.with_metadata(
                    ResolvedBatchMetadata {
                        occurences: source_seq_batch.metadata.occurences.clone()
                    }
                ))?;
            }

            Ok(())
        })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use columnar::pool::{connector_mut, Pool};
    use super::*;

    use crate::{
        crispr::guide::Guide,
        model::{
            alignment::{mined, resolved, MinedBatchMetadata, MinedSchema, ResolvedBatchMetadata, ResolvedSchema},
            cigarx::{Cigarx, CigarxOp, Cigarx64},
            input::{SeqBatchMetadata, SeqSchema, sequences, SEQ_MAX_LEN},
        },
        sequence::iupac::Iupac,
    };

    fn make_stage(capacity: usize) -> AlignmentSimpleResolve {
        AlignmentSimpleResolve {
            pool:  Arc::new(Pool::<ResolvedSchema>::new(16, capacity)),
            arena: Arena::with_capacity(1 << 16),
        }
    }

    /// Build a `[Iupac; SEQ_MAX_LEN]` from an ASCII string, padding remainder with `Iupac::default()`.
    fn make_seq(s: &str) -> [Iupac; SEQ_MAX_LEN] {
        let mut arr = [Iupac::default(); SEQ_MAX_LEN];
        for (i, c) in s.chars().take(SEQ_MAX_LEN).enumerate() {
            arr[i] = Iupac::from_utf8(c);
        }
        arr
    }

    /// Build a `Cigarx64` from a slice of operations.
    fn make_cigarx(ops: &[CigarxOp]) -> Cigarx64 {
        let mut cig = Cigarx64::default();
        for &op in ops { cig.push(op); }
        cig
    }

    /// Build a single-row input batch.
    fn make_input(guide: Guide, content: [Iupac; SEQ_MAX_LEN], cigarx: Cigarx64, offset: u8) 
    -> BatchMut<MinedSchema, MinedBatchMetadata> 
    {
        let seq_pool   = Pool::<SeqSchema>::new(4, 16);
        let mined_pool = Pool::<MinedSchema>::new(4, 16);

        let mut seqs = seq_pool.acquire().unwrap();
        seqs.set_len(1);
        seqs.mutate(
            (sequences::schema::id, sequences::schema::content),
            |(ids, contents)| { 
                ids[0] = 0; contents[0] = content; 
            },
        );

        let seqs_ref = seqs
            .with_metadata(SeqBatchMetadata { 
                seq_len: cigarx.len() as u32, 
                guide, 
                occurences: vec![] 
            })
            .freeze();

        let mut mined = mined_pool.acquire().unwrap();
        mined.set_len(1);
        mined.mutate(
            (mined::schema::seq_id, mined::schema::offset, mined::schema::cigarx),
            |(seq_ids, offsets, cigarxs)| { 
                seq_ids[0] = 0; offsets[0] = offset; cigarxs[0] = cigarx; 
            },
        );
        mined.with_metadata(MinedBatchMetadata { sequences: seqs_ref })
    }

    /// Collect the single output batch from the stage.
    fn run(stage: &mut AlignmentSimpleResolve, input: BatchMut<MinedSchema, MinedBatchMetadata>)
        -> BatchMut<ResolvedSchema, ResolvedBatchMetadata>
    {
        let (mut tx, rx) = connector_mut::<ResolvedSchema, ResolvedBatchMetadata>(4);
        stage.process(input, &mut tx).unwrap();
        drop(tx);
        rx.recv().unwrap()
    }

    /// All Match ops: verify guide and sequence bases are copied at the right positions.
    ///
    ///   Guide:    "ACGT"
    ///   Sequence: "XACGT..."  (offset=1 skips the leading X)
    ///   Cigarx:   ====
    ///
    ///   Expected rguide = "ACGT",  rseq = "ACGT"
    #[test]
    fn match_with_offset_copies_correct_bases() {

        let cigarx = make_cigarx(&[CigarxOp::Match; 4]);
        let input  = make_input(Guide::from("ACGT"), make_seq("NACGT"), cigarx, 1);
        let batch  = run(&mut make_stage(16), input);

        let (rguides, rseqs, lens) = batch.columns((
            resolved::schema::rguide, resolved::schema::rseq, resolved::schema::resolved_len,
        ));

        assert_eq!(lens[0], 4);

        let rg: String = rguides[0][..4].iter().map(|b| b.to_utf8()).collect();
        assert_eq!(rg, "ACGT", "rguide mismatch");

        let rs: String = rseqs[0][..4].iter().map(|b| b.to_utf8()).collect();
        assert_eq!(rs, "ACGT", "rseq mismatch");
    }

    /// A Deletion op inserts a gap in rguide while consuming a sequence base.
    ///
    ///   Guide:    "ACT"   (3 ops without the extra base)
    ///   Sequence: "ACGT"  (sequence has an extra G that the guide doesn't)
    ///   Cigarx:   = = D =
    ///
    ///   Expected rguide: [A, C, gap, T]
    ///   Expected rseq:   [A, C, G,   T]
    #[test]
    fn deletion_inserts_gap_in_rguide() {

        let cigarx = make_cigarx(&[CigarxOp::Match, CigarxOp::Match, CigarxOp::Deletion, CigarxOp::Match]);
        let input  = make_input(Guide::from("ACT"), make_seq("ACGT"), cigarx, 0);
        let batch  = run(&mut make_stage(16), input);

        let (g, s, lens) = batch.columns((
            resolved::schema::rguide, resolved::schema::rseq, resolved::schema::resolved_len,
        ));

        assert_eq!(lens[0], 4);

        assert_eq!(g[0][0], Iupac::from_utf8('A'));
        assert_eq!(g[0][1], Iupac::from_utf8('C'));
        assert_eq!(g[0][2], Iupac::default(), "expected gap in rguide at deletion");
        assert_eq!(g[0][3], Iupac::from_utf8('T'));
        assert_eq!(s[0][0], Iupac::from_utf8('A'));
        assert_eq!(s[0][1], Iupac::from_utf8('C'));
        assert_eq!(s[0][2], Iupac::from_utf8('G'));
        assert_eq!(s[0][3], Iupac::from_utf8('T'));
    }

    /// An Insertion op inserts a gap in rseq while consuming a guide base.
    ///
    ///   Guide:    "ACGT"  (guide has an extra G not present in the sequence)
    ///   Sequence: "ACT"   (sequence has no base at the insertion position)
    ///   Cigarx:   = = I =
    ///
    ///   Expected rguide: [A, C, G,   T]
    ///   Expected rseq:   [A, C, gap, T]   (T is seq[2])
    #[test]
    fn insertion_inserts_gap_in_rseq() {

        let cigarx = make_cigarx(&[CigarxOp::Match, CigarxOp::Match, CigarxOp::Insertion, CigarxOp::Match]);
        let input  = make_input(Guide::from("ACGT"), make_seq("ACT"), cigarx, 0);
        let batch  = run(&mut make_stage(16), input);

        let (g, s, lens) = batch.columns((
            resolved::schema::rguide, resolved::schema::rseq, resolved::schema::resolved_len,
        ));

        assert_eq!(lens[0], 4);

        assert_eq!(g[0][0], Iupac::from_utf8('A'));
        assert_eq!(g[0][1], Iupac::from_utf8('C'));
        assert_eq!(g[0][2], Iupac::from_utf8('G'));
        assert_eq!(g[0][3], Iupac::from_utf8('T'));
        assert_eq!(s[0][0], Iupac::from_utf8('A'));
        assert_eq!(s[0][1], Iupac::from_utf8('C'));
        assert_eq!(s[0][2], Iupac::default(), "expected gap in rseq at insertion");
        assert_eq!(s[0][3], Iupac::from_utf8('T')); // reads seq[2]
    }
}