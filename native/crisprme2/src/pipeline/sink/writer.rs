use std::fmt::Write as FmtWrite;
use std::os::unix::fs::FileExt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{path::{Path, PathBuf}, sync::Arc};
use std::fs::{File, OpenOptions};
use std::io;

use columnar::pipeline::{PipelineError, Sink};
use itertools::izip;

use crate::error::crisprme_errors::ContigLabelsError;
use crate::model::alignment::AlignmentFrame;
use crate::crispr::pam::PAM;

/// Where the PAM sits relative to the protospacer.
///
/// Named rather than a bare `bool` so call sites can't silently invert it —
/// see the `right`/`upstream` mismatch in `sequence::scanner`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PamPlacement {
    /// `<PAM><protospacer>` - e.g. Cas12a `TTTV`.
    Upstream,
    /// `<protospacer><PAM>` — e.g. SpCas9 `NGG`.
    Downstream,
}

impl PamPlacement {
    #[inline]
    pub const fn from_upstream(upstream: bool) -> Self {
        if upstream { Self::Upstream } else { Self::Downstream }
    }
}

/// Run-level, immutable rendering context for the guide column.
///
/// The placement branch is resolved **once**, at pipeline construction, into a
/// prefix/suffix pair of which exactly one is empty. The row loop then does two
/// unconditional `push_str` calls — no per-row branch, no per-row `ParsedPAM`
/// lookup, and the empty side compiles to a length-0 memcpy.
#[derive(Debug, Clone)]
pub struct PamContext {
    prefix: Box<str>,
    suffix: Box<str>,
}

impl PamContext {
    pub fn new(pam: &PAM, placement: PamPlacement) -> Self {
        let motif = pam.motif();
        match placement {
            PamPlacement::Upstream   => Self { prefix: motif.into(), suffix: Box::from("") },
            PamPlacement::Downstream => Self { prefix: Box::from(""), suffix: motif.into() },
        }
    }

    /// Render `<prefix><aligned-guide><suffix>` into `buf`.
    ///
    /// `rguide` is a fixed-width, NUL-terminated ASCII row written by
    /// `Resolver` (bases + `-` for bulges). Bytes past the first NUL are stale.
    #[inline]
    pub(crate) fn render_guide(&self, buf: &mut String, rguide: &[u8]) {
        buf.push_str(&self.prefix);
        for &b in rguide.iter().take_while(|&&b| b != 0) {
            buf.push(b as char);
        }
        buf.push_str(&self.suffix);
    }
}

/// `contig` column renderer.
///
/// Contig ids are dense (`0..len`), assigned by `search._compute_contig_ids`
/// and packed into every `Occurence` by `TargetBatcher::feed_chunk`. The
/// name table is therefore a slice indexed by id (avoid using hashing).
#[derive(Debug, Clone)]
pub enum ContigLabels {
    /// `names[id]` is the FASTA contig name.
    Names(Box<[Box<str>]>),
    /// No mapping supplied (e.g. `dataset_pipeline`): emit the raw id.
    Ids,
}

impl ContigLabels {
    /// Build a name table from `names` **in contig-id order**.
    ///
    /// # Errors
    /// * [`ContigLabelsError::Empty`] — no names supplied.
    /// * [`ContigLabelsError::InvalidName`] — a name is empty, or contains
    ///   `,` `"` `\n` `\r`, which would break the CSV row.
    pub fn from_names(names: Vec<String>) -> Result<Self, ContigLabelsError> {
        if names.is_empty() {
            let err = ContigLabelsError::Empty;
            tracing::error!("{err}");
            return Err(err);
        }

        for (id, name) in names.iter().enumerate() {
            let bad = if name.is_empty() {
                Some(0u8)
            } else {
                name.bytes().find(|b| matches!(b, b',' | b'"' | b'\n' | b'\r'))
            };
            if let Some(byte) = bad {
                let err = ContigLabelsError::InvalidName {
                    id: id as u32,
                    name: name.clone(),
                    byte,
                };
                tracing::error!("{err}");
                return Err(err);
            }
        }

        tracing::info!("contig labels: {} names", names.len());
        Ok(Self::Names(
            names.into_iter().map(String::into_boxed_str).collect(),
        ))
    }

    /// The name for `id`, or `None` when unmapped / out of range.
    #[inline(always)]
    pub fn name(&self, id: u32) -> Option<&str> {
        match self {
            Self::Names(names) => names.get(id as usize).map(|n| &**n),
            Self::Ids => None,
        }
    }

    #[inline(always)]
    pub fn is_named(&self) -> bool {
        matches!(self, Self::Names(_))
    }
}

/// Lock-free multi-threaded CSV writer.
///
/// Each `CsvWriterSink` formats a batch into its own `buffer` (no contention),
/// then atomically claims a byte-range in the file with `fetch_add` and writes
/// it at that offset via `pwrite`.
pub struct CsvWriter {
    offset: AtomicUsize,
    file: File,
    /// Shared, immutable: cloned into each sink at construction.
    pam: PamContext,
    contigs: ContigLabels,
}

impl CsvWriter {
    /// Open the report file, truncating any previous run.
    ///
    /// Returns `io::Error` instead of panicking so the PyO3 layer can surface a
    /// descriptive `OSError` (bad path, no permission, read-only mount).
    pub fn open(path: impl AsRef<Path>, pam: PamContext, contigs: ContigLabels) -> io::Result<Arc<Self>> {
        let path = path.as_ref();
        let file = OpenOptions::new()
            .truncate(true)
            .create(true)
            .write(true)
            .open(path)?;

        tracing::info!("CSV report -> {}", path.display());
        Ok(Arc::new(Self { offset: AtomicUsize::new(0), file, pam, contigs }))
    }

    /// Atomically reserve `bytes` of file space; returns the claimed offset.
    #[inline]
    fn claim(&self, bytes: usize) -> u64 {
        self.offset.fetch_add(bytes, Ordering::Relaxed) as u64
    }
}

pub struct CsvWriterSink {
    inner: Arc<CsvWriter>,
    /// Sink-local copy — keeps the row loop off the shared `Arc` cache line.
    pam: PamContext,
    contigs: ContigLabels,
    /// Fires the unmapped-contig error once per sink, not once per row.
    warned_unmapped: bool,
    /// Per-sink row buffer — formatted here, then pwrite'd atomically.
    buffer: String,
}

impl CsvWriterSink {
    pub fn new(writer: &Arc<CsvWriter>) -> Self {
        Self { 
            inner: writer.clone(), 
            pam: writer.pam.clone(), 
            contigs: writer.contigs.clone(),
            warned_unmapped: false,
            buffer: String::new(),
        }
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
            let features  = cols.features.split();
            let scores    = cols.scores.split();

            let mut feat_iters  = std::array::from_fn::<_, 10, _>(|i| features[i].iter());
            let mut score_iters  = std::array::from_fn::<_, 4, _>(|i| scores[i].iter());

            for (occ, offset, rguide, rseq) in izip!(
                cols.occurence.iter(),
                cols.offset.iter(),
                cols.rguide.iter(),
                cols.rseq.iter(),
            ) {
                let contig_id = occ.contig();

                // contig column
                match self.contigs.name(contig_id) {
                    Some(name) => self.buffer.push_str(name),
                    None => {
                        if self.contigs.is_named() && !self.warned_unmapped {
                            // Invariant break: Python handed rust a name table 
                            // shorter than the ids the batcher packed. Degrade
                            // to the id rather than panicking! 
                            // a panic here kills the sink worker silently and
                            // the report ends up empty.
                            tracing::error!(
                                "contig id {} has no name (table holds {} entries); \
                                falling back to numeric ids for this sink",
                                contig_id,
                                match &self.contigs {
                                    ContigLabels::Names(n) => n.len(),
                                    ContigLabels::Ids => 0,
                                }
                            );
                            self.warned_unmapped = true;
                        }
                        write!(self.buffer, "{contig_id}").unwrap();
                    }
                }

                // Layout is owned by `Occurence`; never unpack the u64 here.
                write!(
                    self.buffer,
                    ",{},{},{}",
                    occ.position(),
                    occ.strand(),
                    offset,
                ).expect("fmt::Write for String is infallible");

                // Aligned guide (PAM decorated) columns
                self.buffer.push(',');
                self.pam.render_guide(&mut self.buffer, rguide);

                // Aligned target columns
                self.buffer.push(',');
                for &b in rseq.iter().take_while(|&&b| b != 0) {
                    self.buffer.push(b as char);
                }

                /*
                for it in &mut feat_iters {
                    write!(self.buffer, ",{}", it.next().unwrap()).unwrap();
                }
                
                
                for it in &mut score_iters {
                    write!(self.buffer, ",{:.6}", it.next().unwrap()).unwrap();
                }
                */

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

