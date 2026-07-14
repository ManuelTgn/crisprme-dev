use std::{fs::File, io::BufReader, io::Read, path::PathBuf};

use columnar::{
    pipeline::{PipelineError, Source},
    MemoryPool,
};
use itertools::izip;

use crate::{
    alignment::thresholds::Thresholds,
    crispr::guide::Guide,
    model::{
        input::{SeqBatch, SeqFrame, SeqOccFrame},
        occurence::Occurence,
    },
    sequence::sequence::Sequence,
};

use crate::model::occurence::Strand;

/// Allocates and fills a `SeqFrame` with exactly `n` sequences read from the binary stream.
/// Returns `None` on clean EOF.
pub trait ReadSeqFrame {
    fn read(&mut self, pool: &MemoryPool, n: usize) -> Option<SeqFrame>;
}

/// Allocates and fills a `SeqOccFrame` from a length-prefixed-per-sequence position stream.
/// Also returns the actual number of sequences read, to synchronise with `ReadSeqFrame`.
pub trait ReadOccFrame {
    /// Read positions for up to `n` sequences, allocating a frame sized to the actual total.
    /// Returns `None` on clean EOF, otherwise `(frame, actual_n_seqs)`.
    fn read(&mut self, pool: &MemoryPool, n: usize) -> Option<SeqOccFrame>;
}

pub struct BinarySequenceReader {
    reader: BufReader<File>,
    pub sequence_len: usize,
    total: usize,
    read: usize,
}

impl BinarySequenceReader {
    pub fn open(path: PathBuf, sequence_len: usize) -> Self {
        let file = File::open(path).expect("File not found");
        let file_size = file
            .metadata()
            .expect("File metadata could not be read")
            .len();
        assert!(
            file_size as usize % sequence_len == 0,
            "sequence file size ({file_size} bytes) is not divisible by sequence_len ({sequence_len})"
        );
        Self {
            reader: BufReader::new(file),
            total: file_size as usize / sequence_len,
            sequence_len,
            read: 0,
        }
    }
}

impl ReadSeqFrame for BinarySequenceReader {
    fn read(&mut self, pool: &MemoryPool, n: usize) -> Option<SeqFrame> {
        assert!(
            self.sequence_len <= 32,
            "seq_len must be <= SEQ_MAX_LEN (32)"
        );

        if self.read >= self.total {
            return None;
        }

        let to_read = (self.total - self.read).min(n);
        let mut frame = SeqFrame::alloc(pool, to_read);
        frame.with_cols(|mut cols| {
            for row in cols.content.iter_mut() {
                let bytes = bytemuck::cast_slice_mut(&mut row[..self.sequence_len]);
                self.reader
                    .read_exact(bytes)
                    .expect("error reading sequence");
            }
        });

        self.read += to_read;
        Some(frame)
    }
}

pub struct BinaryPositionReader {
    reader: BufReader<File>,
    /// Scratch buffer for counts; reused across batches.
    counts: Vec<u32>,
    /// Scratch buffer for positions; reused across batches
    positions: Vec<u32>,
}

impl BinaryPositionReader {
    pub fn open(path: PathBuf) -> Self {
        let file = File::open(path).expect("position file not found");
        Self {
            reader: BufReader::new(file),
            counts: vec![],
            positions: vec![],
        }
    }
}

impl ReadOccFrame for BinaryPositionReader {
    fn read(&mut self, pool: &MemoryPool, n: usize) -> Option<SeqOccFrame> {
        self.positions.clear();
        self.counts.clear();

        for _ in 0..n {
            let mut len_buf = [0u8; 4];
            if let Err(e) = self.reader.read_exact(&mut len_buf) {
                if e.kind() == std::io::ErrorKind::UnexpectedEof {
                    break;
                } else {
                    panic!("unable to read position count: {:?}", e);
                }
            }

            let len = u32::from_le_bytes(len_buf) as usize;
            let offset = self.positions.len();

            self.positions.resize(offset + len, 0);
            let bytes = bytemuck::cast_slice_mut(&mut self.positions[offset..]);
            self.reader
                .read_exact(bytes)
                .expect("unable to read positions");

            self.counts.push(len as u32);
        }

        if self.counts.is_empty() {
            return None;
        }

        let total: usize = self.counts.iter().map(|&c| c as usize).sum();
        let mut frame = SeqOccFrame::alloc(pool, total);
        let counts = &self.counts;
        let positions = &self.positions;

        frame.with_cols(|mut cols| {
            let mut pos_offset = 0;
            let mut occ_iter = izip!(cols.seq_row_idx.iter_mut(), cols.occurence.iter_mut());

            for (seq_idx, &count) in counts.iter().enumerate() {
                let count = count as usize;
                for &pos in &positions[pos_offset..pos_offset + count] {
                    let (idx, occ) = occ_iter.next().expect("occ frame size mismatch");

                    *occ = Occurence::new(0, pos, Strand::from_bit(0));
                    *idx = seq_idx as u32;
                }
                pos_offset += count;
            }
        });

        Some(frame)
    }
}

pub struct Reader {
    pool: MemoryPool,
    guide: Guide,
    thresholds: Thresholds,
    batch_size: usize,
    sequences: BinarySequenceReader,
    positions: BinaryPositionReader,
    total_sequences: usize,
}

impl Reader {
    pub fn open(
        seq_path: PathBuf,
        pos_path: PathBuf,
        sequence_len: usize,
        batch_size: usize,
        guide: Guide,
        thresholds: Thresholds,
        pool: MemoryPool,
    ) -> std::io::Result<Self> {
        Ok(Self {
            pool,
            guide,
            thresholds,
            batch_size,
            sequences: BinarySequenceReader::open(seq_path, sequence_len),
            positions: BinaryPositionReader::open(pos_path),
            total_sequences: 0,
        })
    }
}

impl Source for Reader {
    type O = SeqBatch;

    fn name() -> &'static str {
        "Reader"
    }

    fn next(&mut self) -> Result<Option<Self::O>, PipelineError> {
        let mut occs = match self.positions.read(&self.pool, self.batch_size) {
            None => return Ok(None),
            Some(frame) => frame,
        };
        let n_seqs = self.positions.counts.len();
        let mut seqs = self
            .sequences
            .read(&self.pool, n_seqs)
            .expect("position/sequence count mismatch");

        /*
        seqs.with_cols(|cols| {
            for (i, element) in cols.content.iter().enumerate() {
                let element = Sequence::new(element);
                println!("{element:?}");
                occs.with_cols(|o| {
                    for (seq_row_idx, occ) in o.seq_row_idx.iter().zip(o.occurence.iter()) {
                        if *seq_row_idx == i as u32 {
                            println!("  occurs at position {}, seq_row_idx {}", occ.position(), *seq_row_idx);
                        }
                    }
                });
            }
        });
        */

        self.total_sequences += n_seqs;
        println!(
            "Read batch: {} sequences, total so far: {}",
            n_seqs, self.total_sequences
        );

        Ok(Some(SeqBatch {
            seq_len: self.sequences.sequence_len,
            guide: self.guide.clone(),
            thresholds: self.thresholds,
            sequences: seqs,
            occurences: occs,
        }))
    }
}

#[cfg(test)]
mod test {
    use super::{BinaryPositionReader, BinarySequenceReader, ReadOccFrame, ReadSeqFrame, Reader};
    use crate::{
        alignment::thresholds::Thresholds, crispr::guide::Guide, pipeline::test::make_pool,
    };
    use columnar::pipeline::Source;
    use std::{io::Write, path::PathBuf};

    /// Owns a temp file path and deletes it on drop.
    struct TempFile(PathBuf);

    impl TempFile {
        fn new() -> Self {
            let path =
                std::env::temp_dir().join(format!("crisprme_test_{}", rand::random::<u64>()));
            Self(path)
        }

        fn path(&self) -> PathBuf {
            self.0.clone()
        }

        fn writer(&self) -> std::fs::File {
            std::fs::File::create(&self.0).unwrap()
        }
    }

    impl Drop for TempFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    /// Write sequences as flat raw bytes (one seq_len-byte row per sequence).
    fn write_seqs(f: &mut impl Write, seqs: &[Vec<u8>]) {
        for s in seqs {
            f.write_all(s).unwrap();
        }
    }

    /// Write positions as length-prefixed u32 LE arrays (one per sequence).
    fn write_positions(f: &mut impl Write, data: &[Vec<u32>]) {
        for pos_array in data {
            f.write_all(&(pos_array.len() as u32).to_le_bytes())
                .unwrap();
            for &p in pos_array {
                f.write_all(&p.to_le_bytes()).unwrap();
            }
        }
    }

    /// Build a seq of `seq_len` bytes where byte[i] = (seq_idx + i) % 256.
    fn make_seq(seq_idx: usize, seq_len: usize) -> Vec<u8> {
        (0..seq_len).map(|i| ((seq_idx + i) % 256) as u8).collect()
    }

    /// Decode the genomic position from an Occurence (encoded as `position << 1`).
    fn occ_position(occ: &crate::model::occurence::Occurence) -> u32 {
        ((occ.0 >> 1) & 0xFFFF_FFFF) as u32
    }

    #[test]
    fn seq_reader_reads_all_in_one_batch() {
        let seq_len = 4;
        let seqs: Vec<Vec<u8>> = (0..3).map(|i| make_seq(i, seq_len)).collect();
        let tmp = TempFile::new();
        write_seqs(&mut tmp.writer(), &seqs);

        let pool = make_pool();
        let mut reader = BinarySequenceReader::open(tmp.path(), seq_len);
        assert_eq!(reader.total, 3);

        let mut frame = reader.read(&pool, 10).unwrap();
        frame.with_cols(|cols| {
            for (i, row) in cols.content.iter().enumerate() {
                let got: Vec<u8> = bytemuck::cast_slice(&row[..seq_len]).to_vec();
                assert_eq!(got, make_seq(i, seq_len), "sequence {i} mismatch");
            }
        });

        assert!(reader.read(&pool, 10).is_none());
    }

    #[test]
    fn seq_reader_batches_correctly() {
        let seq_len = 4;
        let seqs: Vec<Vec<u8>> = (0..5).map(|i| make_seq(i, seq_len)).collect();
        let tmp = TempFile::new();
        write_seqs(&mut tmp.writer(), &seqs);

        let pool = make_pool();
        let mut reader = BinarySequenceReader::open(tmp.path(), seq_len);

        let mut b1 = reader.read(&pool, 3).unwrap();
        let mut b2 = reader.read(&pool, 3).unwrap(); // only 2 remain
        assert!(reader.read(&pool, 3).is_none());

        b1.with_cols(|cols| assert_eq!(cols.content.iter().count(), 3));
        b2.with_cols(|cols| assert_eq!(cols.content.iter().count(), 2));
    }

    #[test]
    fn pos_reader_assigns_seq_row_idx_and_position_correctly() {
        // seq 0 → [10, 20],  seq 1 → [30],  seq 2 → [40, 50, 60]
        let data: Vec<Vec<u32>> = vec![vec![10, 20], vec![30], vec![40, 50, 60]];
        let tmp = TempFile::new();
        write_positions(&mut tmp.writer(), &data);

        let pool = make_pool();
        let mut reader = BinaryPositionReader::open(tmp.path());
        let mut frame = reader.read(&pool, 10).unwrap();

        let expected_idx: &[u32] = &[0, 0, 1, 2, 2, 2];
        let expected_pos: &[u32] = &[10, 20, 30, 40, 50, 60];
        frame.with_cols(|cols| {
            let idxs: Vec<u32> = cols.seq_row_idx.iter().copied().collect();
            let pos: Vec<u32> = cols.occurence.iter().map(occ_position).collect();
            assert_eq!(idxs, expected_idx);
            assert_eq!(pos, expected_pos);
        });
    }

    #[test]
    fn pos_reader_batches_and_eof() {
        let data: Vec<Vec<u32>> = vec![vec![1], vec![2, 3], vec![4]];
        let tmp = TempFile::new();
        write_positions(&mut tmp.writer(), &data);

        let pool = make_pool();
        let mut reader = BinaryPositionReader::open(tmp.path());

        let mut b1 = reader.read(&pool, 2).unwrap(); // seqs 0+1 → 3 occs
        b1.with_cols(|cols| assert_eq!(cols.seq_row_idx.iter().count(), 3));

        let mut b2 = reader.read(&pool, 2).unwrap(); // seq 2 → 1 occ
        b2.with_cols(|cols| assert_eq!(cols.seq_row_idx.iter().count(), 1));

        assert!(reader.read(&pool, 2).is_none());
    }

    #[test]
    fn reader_full_pipeline_two_batches() {
        let seq_len = 4;
        let seqs: Vec<Vec<u8>> = (0..4).map(|i| make_seq(i, seq_len)).collect();
        let data: Vec<Vec<u32>> = vec![vec![100], vec![200, 201], vec![300], vec![400, 401, 402]];

        let seq_tmp = TempFile::new();
        let pos_tmp = TempFile::new();
        write_seqs(&mut seq_tmp.writer(), &seqs);
        write_positions(&mut pos_tmp.writer(), &data);

        let pool = make_pool();
        let mut reader = Reader::open(
            seq_tmp.path(),
            pos_tmp.path(),
            seq_len,
            2,
            Guide::new("ACGT"),
            Thresholds::new(1, 1, 2),
            pool,
        )
        .unwrap();

        // Batch 1: seqs 0,1 — 1+2 = 3 occs
        let mut b1 = reader.next().unwrap().unwrap();
        assert_eq!(b1.seq_len, seq_len);
        b1.occurences
            .with_cols(|cols| assert_eq!(cols.seq_row_idx.iter().count(), 3));

        // Batch 2: seqs 2,3 — 1+3 = 4 occs
        let mut b2 = reader.next().unwrap().unwrap();
        b2.occurences
            .with_cols(|cols| assert_eq!(cols.seq_row_idx.iter().count(), 4));

        assert!(reader.next().unwrap().is_none());
    }
}
