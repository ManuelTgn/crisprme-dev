use std::fmt::Write as FmtWrite;
use std::os::unix::fs::FileExt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{path::PathBuf, sync::Arc};
use std::fs::{File, OpenOptions};

use columnar::pipeline::{PipelineError, Sink};
use itertools::izip;

use crate::model::alignment::AlignmentFrame;

/// Lock-free multi-threaded CSV writer.
///
/// Each `CsvWriterSink` formats a batch into its own `buffer` (no contention),
/// then atomically claims a byte-range in the file with `fetch_add` and writes
/// it at that offset via `pwrite`.
pub struct CsvWriter {
    offset: AtomicUsize,
    file: File,
}

impl CsvWriter {
    pub fn open(path: PathBuf) -> Arc<Self> {
        Arc::new(Self {
            offset: AtomicUsize::new(0),
            file: OpenOptions::new()
                .truncate(true)
                .create(true)
                .write(true)
                .open(path)
                .expect("unable to open write file")
        })
    }
}

pub struct CsvWriterSink {
    inner: Arc<CsvWriter>,
    /// Per-sink row buffer — formatted here, then pwrite'd atomically.
    buffer: String,
}

impl CsvWriterSink {
    pub fn new(writer: &Arc<CsvWriter>) -> Self {
        Self { inner: writer.clone(), buffer: String::new() }
    }
}

impl Sink for CsvWriterSink {
    type I = AlignmentFrame;

    fn name() -> &'static str { "CsvWriter" }

    fn consume(&mut self, mut item: Self::I) -> Result<(), PipelineError> {
        self.buffer.clear();

        item.with_cols(|mut cols| {

            // split() gives [Column<T>; N]; from_fn builds an array of iterators
            // that can be advanced one element at a time in a simple loop per row.
            let features = cols.features.split();
            let scores    = cols.scores.split();

            let mut feat_iters  = std::array::from_fn::<_, 10, _>(|i| features[i].iter());
            let mut score_iters  = std::array::from_fn::<_, 4, _>(|i| scores[i].iter());

            for (occ, offset, rguide, rseq) in izip!(
                cols.occurence.iter(),
                cols.offset.iter(),
                cols.rguide.iter(),
                cols.rseq.iter(),
            ) {
                let contig   = (occ.0 >> 33) as u32;
                let position = ((occ.0 >> 1) & 0xFFFF_FFFF) as u32;
                let strand    = (occ.0 & 1) as u8;

                write!(self.buffer, "{},{},{},{}", contig, position, strand, offset).unwrap();

                self.buffer.push(',');
                for &b in rguide.iter().take_while(|&&b| b != 0) {
                    self.buffer.push(b as char);
                }

                self.buffer.push(',');
                for &b in rseq.iter().take_while(|&&b| b != 0) {
                    self.buffer.push(b as char);
                }

                for it in &mut feat_iters {
                    write!(self.buffer, ",{}", it.next().unwrap()).unwrap();
                }
                
                for it in &mut score_iters {
                    write!(self.buffer, ",{:.6}", it.next().unwrap()).unwrap();
                }

                self.buffer.push('\n');
            }
        });

        if !self.buffer.is_empty() {
            let offset = self.inner.offset.fetch_add(self.buffer.len(), Ordering::Relaxed);
            self.inner.file
                .write_at(self.buffer.as_bytes(), offset as u64)
                .expect("csv write failed");
        }

        Ok(())
    }
}

