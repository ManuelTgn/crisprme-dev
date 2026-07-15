use std::fmt::Write as FmtWrite;
use std::fs::{File, OpenOptions};
use std::io;
use std::os::unix::fs::FileExt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{
    path::{Path, PathBuf},
    sync::Arc,
};

use columnar::pipeline::{PipelineError, Sink};
use itertools::izip;

use crate::crispr::pam::PAM;
use crate::error::crisprme_errors::ContigLabelsError;
use crate::model::alignment::AlignmentFrame;
use crate::model::occurence::Strand;

/// `Occurence::pam_id()` sentinel meaning "no concrete PAM was recorded" — the
/// target is then decorated with the degenerate motif instead. Produced by the
/// raw-reader path and by `feed_chunk` when the reference base under the PAM is
/// ambiguous.
pub const PAM_ID_NONE: u16 = u16::MAX;

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
        if upstream {
            Self::Upstream
        } else {
            Self::Downstream
        }
    }

    #[inline(always)]
    pub const fn is_upstream(self) -> bool {
        matches!(self, Self::Upstream)
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
    /// Guide-column decoration: the degenerate motif, pre-split so exactly one
    /// of prefix/suffix is non-empty (unchanged behaviour).
    prefix: Box<str>,
    suffix: Box<str>,
    placement: PamPlacement,
    /// Target-column decoration: concrete PAM variants indexed by pam id
    /// (reported-strand orientation), e.g. `variants[3] == "TGG"` for `NGG`.
    variants: Box<[Box<str>]>,
    /// Fallback when an occurrence carries no concrete PAM (`PAM_ID_NONE`, or an
    /// id past the table): the degenerate motif.
    motif: Box<str>,
}

impl PamContext {
    pub fn new(pam: &PAM, placement: PamPlacement) -> Self {
        let motif: Box<str> = pam.motif().into();

        // Guide decoration (degenerate motif) — resolved once into a
        // prefix/suffix pair, exactly as before.
        let (prefix, suffix) = match placement {
            PamPlacement::Upstream => (motif.clone(), Box::<str>::from("")),
            PamPlacement::Downstream => (Box::<str>::from(""), motif.clone()),
        };

        // Target decoration (concrete variants) — one small string per pam id.
        // `variant_count <= 65_536` and every id in range is a valid variant.
        let variants: Box<[Box<str>]> = (0..pam.variant_count())
            .map(|id| {
                pam.pam_variant_ascii(id as u16)
                    .expect("id < variant_count is always a valid variant")
                    .into_boxed_str()
            })
            .collect();

        Self {
            prefix,
            suffix,
            placement,
            variants,
            motif,
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

    /// Render the aligned target `rseq`, decorated with its **concrete** matched
    /// PAM, into `buf`.
    ///
    /// Mirrors [`render_guide`], but the decoration is the concrete PAM variant
    /// (e.g. `TGG`) selected by `pam_id`, not the degenerate motif. An
    /// out-of-range id — [`PAM_ID_NONE`], or any occurrence with no concrete PAM
    /// (raw-reader path / ambiguous reference base) — falls back to the motif so
    /// the column width stays consistent. Placement matches the guide: the PAM
    /// is prepended when upstream, appended when downstream.
    ///
    /// `pam_id` is stored by `feed_chunk` in reported-strand orientation, so the
    /// variant needs no strand handling here.
    #[inline]
    pub(crate) fn render_target(&self, buf: &mut String, rseq: &[u8], pam_id: u16) {
        let pam: &str = self
            .variants
            .get(pam_id as usize)
            .map(|v| &**v)
            .unwrap_or(&self.motif);

        if self.placement.is_upstream() {
            buf.push_str(pam);
        }
        for &b in rseq.iter().take_while(|&&b| b != 0) {
            buf.push(b as char);
        }
        if !self.placement.is_upstream() {
            buf.push_str(pam);
        }
    }

    #[inline(always)]
    pub fn target_start(&self, window_fwd_left: u32, offset: u8, strand: Strand) -> u32 {
        if strand.scanned_on_revcomp(self.placement.is_upstream()) {
            window_fwd_left
        } else {
            window_fwd_left + offset as u32
        }
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
                name.bytes()
                    .find(|b| matches!(b, b',' | b'"' | b'\n' | b'\r'))
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
    pub fn name(&self, id: u16) -> Option<&str> {
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
    pub fn open(
        path: impl AsRef<Path>,
        pam: PamContext,
        contigs: ContigLabels,
    ) -> io::Result<Arc<Self>> {
        let path = path.as_ref();
        let file = OpenOptions::new()
            .truncate(true)
            .create(true)
            .write(true)
            .open(path)?;

        tracing::info!("CSV report -> {}", path.display());
        Ok(Arc::new(Self {
            offset: AtomicUsize::new(0),
            file,
            pam,
            contigs,
        }))
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

    fn name() -> &'static str {
        "CsvWriter"
    }

    fn consume(&mut self, mut item: Self::I) -> Result<(), PipelineError> {
        self.buffer.clear();

        item.with_cols(|mut cols| {
            // split() gives [Column<T>; N]; from_fn builds an array of iterators
            // that can be advanced one element at a time in a simple loop per row.
            let features = cols.features.split();
            let scores = cols.scores.split();

            let mut feat_iters = std::array::from_fn::<_, 10, _>(|i| features[i].iter());
            let mut score_iters = std::array::from_fn::<_, 4, _>(|i| scores[i].iter());

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

                let strand = occ.strand();
                let position = self.pam.target_start(occ.position(), *offset, strand) + 1;

                // Layout is owned by `Occurence`; never unpack the u64 here.
                write!(self.buffer, ",{},{}", position, strand)
                    .expect("fmt::Write for String is infallible");

                // Aligned guide (PAM decorated) columns
                self.buffer.push(',');
                self.pam.render_guide(&mut self.buffer, rguide);

                // Aligned target columns
                self.buffer.push(',');
                self.pam.render_target(&mut self.buffer, rseq, occ.pam());

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
            let offset = self
                .inner
                .offset
                .fetch_add(self.buffer.len(), Ordering::Relaxed);
            self.inner
                .file
                .write_at(self.buffer.as_bytes(), offset as u64)
                .expect("csv write failed");
        }

        Ok(())
    }
}
