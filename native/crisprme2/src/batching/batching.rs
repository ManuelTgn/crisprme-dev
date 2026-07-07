//! Window batching: dedup of protospacer windows and accumulation of their
//! genomic occurrences.
//!
//! [`TargetBatcher`] drives the per-chunk flow: encode the chunk to IUPAC
//! bitmasks, run the parallel PAM/target scanner, then for every accepted hit
//! split the window into a **canonical protospacer** (the map key) and a
//! **PAM variant index**, and accumulate `protospacer -> [occurrence]` into an
//! in-memory map. When a size threshold is crossed the map is flushed to a
//! [`WindowBatch`] for the downstream aligner.
//!
//! # What changed with PAM stripping
//!
//! The map key is now the protospacer **without** the PAM (length
//! `size - plen`), produced in canonical 5'→3' orientation by
//! [`TargetExtractor`]. Two consequences:
//! * occurrences that differ only in their PAM variant, or only in strand,
//!   collapse onto the same key (more dedup), and
//! * each occurrence records which PAM variant it used, as a `u16` index,
//!   alongside its packed `(contig, pos, strand)` — see [`OccRecord`].
//!
//! Performance note: the protospacer bytes are built into a single reusable
//! buffer per chunk; a boxed map key is allocated only when a genuinely new
//! unique window is inserted (repeat hits take the allocation-free
//! `get_mut(&[u8])` path).

use std::sync::atomic::{AtomicUsize, Ordering};

use crate::crispr::guide::Guide;
use crate::crispr::target::TargetExtractor;
use crate::crispr::{guide, pam};
use crate::memory::batch::AlignmentRingBatch;
use crate::sequence::{iupac, scanner};

use ahash::AHashMap;

use crossbeam_channel::Receiver;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Key: owned protospacer bytes (IUPAC bitmasks), length == `size - plen`,
/// in canonical 5'→3' orientation.
type WindowKey = Box<[u8]>;

/// Packed `(contig_id, pos, strand_bit)` in a single `u64`.
///
/// Layout: `[ contig_id : 31 ][ pos : 32 ][ strand : 1 ]`, i.e.
/// `occ = (contig_id << 33) | (pos << 1) | strand`.
type Occ = u64;

/// A single genomic occurrence of a protospacer window.
///
/// Wraps the packed positional data together with the PAM variant index for
/// that occurrence. Kept deliberately small (`u64 + u16`) and `Copy` so the
/// per-window occurrence vectors stay cache-friendly.
#[derive(Clone, Copy, Debug)]
pub struct OccRecord {
    /// Packed `(contig_id, pos, strand)`; see [`pack_occ`] / [`unpack_occ`].
    pub packed: Occ,
    /// PAM variant index in canonical orientation (see [`ParsedPAM::pam_index`]).
    pub pam_idx: u16,
}

#[inline(always)]
fn pack_occ(contig_id: u32, pos: u32, strand_bit: u8) -> Occ {
    ((contig_id as u64) << 33) | ((pos as u64) << 1) | ((strand_bit as u64) & 1)
}

#[inline(always)]
pub fn unpack_occ(occ: Occ) -> (u32, u32, u8) {
    let contig_id = (occ >> 33) as u32;
    let pos = ((occ >> 1) & 0xFFFF_FFFF) as u32;
    let strand_bit = (occ & 1) as u8;
    (contig_id, pos, strand_bit)
}

#[pyclass]
#[derive(Clone)]
pub struct BatcherStats {
    #[pyo3(get)]
    pub hits_in_batch: usize,
    #[pyo3(get)]
    pub unique_windows: usize,
}

#[pyclass]
#[derive(Clone)]
pub struct FeedStatus {
    #[pyo3(get)]
    pub flushed: bool,
    #[pyo3(get)]
    pub stats: BatcherStats,
}

static TARGET_BATCHER_NEXT_ID: AtomicUsize = AtomicUsize::new(0);

/// Accumulates unique protospacer windows and their occurrences, flushing to
/// the pipeline when a size threshold is crossed.
#[pyclass]
pub struct TargetBatcher {
    #[pyo3(get)]
    id: usize,

    // config
    size: usize,
    right: bool,
    threads: usize,
    batch_hits: usize,
    max_unique: usize,
    overlap_left: usize,

    // Stream of completed alignment batches
    alignment_rx: Option<Receiver<AlignmentRingBatch>>,

    // parsed PAM
    pam: pam::ParsedPAM,

    // guide
    guide: guide::Guide,

    // precomputed protospacer/PAM split geometry (built once from pam + size)
    extractor: TargetExtractor,

    // state
    map: AHashMap<WindowKey, Vec<OccRecord>>,
    hits_in_batch: usize,
}

#[pymethods]
impl TargetBatcher {
    #[new]
    pub fn new(
        pam_seq: &str,
        guide_seq: &str,
        size: usize,
        right: bool,
        threads: usize,
        batch_hits: usize,
        max_unique: usize,
        overlap_left: usize,
    ) -> PyResult<Self> {
        let pam = pam::ParsedPAM::new(pam_seq)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("Invalid PAM sequence: {e}")))?;

        let guide = guide::Guide::from(guide_seq);

        if size > 0 && overlap_left < size.saturating_sub(1) {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Invalid overlap_left={overlap_left}: must be >= size-1={} to avoid losing kmers at chunk boundaries",
                size.saturating_sub(1)
            )));
        }

        // Build the protospacer/PAM split geometry once. Maps TargetError to a
        // descriptive Python exception via `From<TargetError> for PyErr`.
        let extractor = TargetExtractor::new(pam.plen(), size, right)?;

        let id = TARGET_BATCHER_NEXT_ID.fetch_add(1, Ordering::SeqCst);

        tracing::info!(
            "TargetBatcher #{id} ready: pam={pam_seq:?}, size={size}, \
             proto_len={}, pam_variants={}, threads={threads}",
            extractor.proto_len(),
            pam.variant_count(),
        );

        Ok(Self {
            id,
            alignment_rx: None,
            size,
            right,
            threads,
            batch_hits,
            max_unique,
            overlap_left,
            pam,
            guide,
            extractor,
            map: AHashMap::new(),
            hits_in_batch: 0,
        })
    }

    pub fn feed_chunk(
        &mut self,
        contig_id: u32,
        chunk_start: u32,
        chunk_seq: &str,
        valid_len: usize,
    ) -> PyResult<FeedStatus> {
        let seq_bitmask: Vec<u8> = iupac::sequence_encoder(chunk_seq);

        let (pos_local, strand) = scanner::scan_targets_bitmask(
            &seq_bitmask,
            &self.pam,
            self.size,
            self.right,
            self.threads,
        )
        .map_err(|e| PyErr::new::<PyValueError, _>(e))?;

        debug_assert_eq!(pos_local.len(), strand.len());

        let chunk_len = seq_bitmask.len();
        if self.size == 0 || chunk_len < self.size {
            return Ok(self.feed_status(false));
        }

        let max_start_excl = chunk_len - self.size + 1;
        let core_len = valid_len;

        // Accept only hits whose window start falls in this chunk's "core"
        // region, so overlapping chunk boundaries do not double-count kmers.
        let (accept_lo, mut accept_hi) = if chunk_start == 0 {
            (0usize, core_len)
        } else {
            let ov = self.overlap_left;
            let recovery = self.size.saturating_sub(1);
            let lo = ov.saturating_sub(recovery);
            let hi = ov + core_len;
            (lo, hi)
        };

        if accept_hi > max_start_excl {
            accept_hi = max_start_excl;
        }

        if accept_hi <= accept_lo {
            return Ok(self.feed_status(self.should_flush()));
        }

        // Copy of the (Copy) geometry so the hot loop borrows only `self.pam`
        // and `self.map`, not the whole `self`.
        let extractor = self.extractor;
        let mut proto: Vec<u8> = Vec::with_capacity(extractor.proto_len());
        let mut accepted = 0usize;

        for i in 0..pos_local.len() {
            let p = pos_local[i];
            if p < accept_lo || p >= accept_hi {
                continue;
            }

            let pos_global = chunk_start as usize + p;
            if pos_global > (u32::MAX as usize) {
                return Err(PyErr::new::<PyValueError, _>("Position overflow"));
            }

            let strand_bit = strand[i]; // 1 = fwd (+), 0 = rev (-)
            let window = &seq_bitmask[p..p + self.size];

            // Split into canonical protospacer (into `proto`) + PAM index.
            let pam_idx = extractor.extract(&self.pam, window, strand_bit, &mut proto);
            let rec = OccRecord {
                packed: pack_occ(contig_id, pos_global as u32, strand_bit),
                pam_idx,
            };

            // Allocate a boxed key only for a genuinely new unique window.
            if let Some(v) = self.map.get_mut(proto.as_slice()) {
                v.push(rec);
            } else {
                self.map
                    .insert(proto.as_slice().to_vec().into_boxed_slice(), vec![rec]);
            }

            self.hits_in_batch += 1;
            accepted += 1;
        }

        tracing::debug!(
            "feed_chunk contig={contig_id} start={chunk_start} raw_hits={} \
             accepted={accepted} unique_windows={} hits_in_batch={}",
            pos_local.len(),
            self.map.len(),
            self.hits_in_batch,
        );

        Ok(self.feed_status(self.should_flush()))
    }

    pub fn flush_and_align(&mut self, _max_mm: usize, _bdna: usize, _brna: usize) -> PyResult<()> {
        // Collect window batch on flush.
        let _batch: WindowBatch = self.flush_to_batch();
        // TODO: dispatch `_batch` to the aligner.
        Ok(())
    }

    /// Flush remaining data at end of genome. Returns stats of what was
    /// flushed (and clears internal state).
    pub fn finalize(&mut self) -> PyResult<BatcherStats> {
        let stats = BatcherStats {
            hits_in_batch: self.hits_in_batch,
            unique_windows: self.map.len(),
        };
        tracing::debug!(
            "finalize batcher #{}: clearing {} residual windows ({} hits)",
            self.id,
            stats.unique_windows,
            stats.hits_in_batch,
        );
        self.clear_batch();
        Ok(stats)
    }

    /// Introspection (debug).
    pub fn stats(&self) -> PyResult<BatcherStats> {
        Ok(BatcherStats {
            hits_in_batch: self.hits_in_batch,
            unique_windows: self.map.len(),
        })
    }
}

impl TargetBatcher {
    pub fn id(&self) -> usize {
        self.id
    }

    pub fn get_window_count(&self) -> usize {
        self.map.len()
    }

    pub fn get_window_keys(&self) -> impl Iterator<Item = &WindowKey> {
        self.map.keys()
    }

    pub fn extract_alignment_rx(&mut self) -> Option<Receiver<AlignmentRingBatch>> {
        self.alignment_rx.take()
    }

    /// Build a [`FeedStatus`] snapshot from the current counters.
    #[inline(always)]
    fn feed_status(&self, flushed: bool) -> FeedStatus {
        FeedStatus {
            flushed,
            stats: BatcherStats {
                hits_in_batch: self.hits_in_batch,
                unique_windows: self.map.len(),
            },
        }
    }

    /// Convert the current batch (unique windows + occurrences) into a
    /// [`WindowBatch`] and clear the batch's hit counter.
    pub fn flush_to_batch(&mut self) -> WindowBatch {
        let map: AHashMap<WindowKey, Vec<OccRecord>> = std::mem::take(&mut self.map);

        let unique = map.len();
        let mut windows: Vec<WindowKey> = Vec::with_capacity(unique);
        let mut occs: Vec<Vec<OccRecord>> = Vec::with_capacity(unique);

        let mut total_hits = 0usize;
        for (k, v) in map {
            total_hits += v.len();
            windows.push(k);
            occs.push(v);
        }

        self.hits_in_batch = 0;

        tracing::info!(
            "flushing batch: {unique} unique windows, {total_hits} occurrences"
        );

        WindowBatch {
            windows,
            occs,
            total_hits,
        }
    }

    #[inline(always)]
    fn should_flush(&self) -> bool {
        self.hits_in_batch >= self.batch_hits || self.map.len() >= self.max_unique
    }

    #[inline(always)]
    fn clear_batch(&mut self) {
        self.map.clear();
        self.hits_in_batch = 0;
    }

    pub fn set_alignment_stream(&mut self, rx: Receiver<AlignmentRingBatch>) {
        self.alignment_rx = Some(rx);
    }

    /// Length of the stored sequence (the protospacer, i.e. `size - plen`).
    ///
    /// This is what downstream frames should size rows to — the PAM has been
    /// stripped out of the stored window.
    pub fn get_sequence_len(&self) -> usize {
        self.extractor.proto_len()
    }

    pub fn get_guide(&self) -> Guide {
        self.guide.clone()
    }
}

/// A flushed batch: unique protospacer windows plus their occurrences.
#[derive(Debug)]
pub struct WindowBatch {
    /// Unique protospacer windows, each length == `size - plen`.
    pub windows: Vec<WindowKey>,
    /// Occurrences for each window (parallel to `windows`).
    pub occs: Vec<Vec<OccRecord>>,
    /// Total occurrences across all windows.
    pub total_hits: usize,
}

impl WindowBatch {
    #[inline]
    pub fn len(&self) -> usize {
        self.windows.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }
}