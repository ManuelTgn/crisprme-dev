use crate::crispr::guide::Guide;
use crate::crispr::{guide, pam};
use crate::memory::batch::AlignmentRingBatch;
use crate::sequence::{iupac, scanner};
use std::sync::atomic::{AtomicUsize, Ordering};

use ahash::AHashMap;

use crossbeam_channel::Receiver;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Key: owned window bytes (IUPAC bitmasks), length == size.
type WindowKey = Box<[u8]>;

/// Occurrence: packed (contig_id, pos, strand_bit) into u64.
/// Layout: [ contig_id:31.. ] [ pos:32 bits ] [ strand:1 bit ]
/// occ = (contig_id << 33) | (pos << 1) | strand
type Occ = u64;

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

/// TargetBatcher class
#[pyclass]
pub struct TargetBatcher {
    #[pyo3(get)]
    id: usize,

    // config
    size: usize,
    upstream: bool,
    threads: usize,
    batch_hits: usize,
    max_unique: usize,
    overlap_left: usize,

    // Stream of completed alignment batches
    alignment_rx: Option<Receiver<AlignmentRingBatch>>,

    // parsed PAM
    pam: pam::PAM,

    // guide
    guide: guide::Guide,

    // state
    map: AHashMap<WindowKey, Vec<Occ>>,
    hits_in_batch: usize,
}

#[pymethods]
impl TargetBatcher {
    #[new]
    pub fn new(
        pam_seq: &str,
        guide_seq: &str,
        size: usize,
        upstream: bool,
        threads: usize,
        batch_hits: usize,
        max_unique: usize,
        overlap_left: usize,
    ) -> PyResult<Self> {
        let pam = pam::PAM::new(pam_seq)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("Invalid PAM sequence: {e}")))?;

        let guide = guide::Guide::from(guide_seq);

        if size > 0 && overlap_left < size.saturating_sub(1) {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Invalid overlap_left={overlap_left}: must be >= size-1={} to avoid losing kmers at chunk boundaries",
                size.saturating_sub(1)
            )));
        }

        Ok(Self {
            id: TARGET_BATCHER_NEXT_ID.fetch_add(1, Ordering::SeqCst),
            alignment_rx: None,
            size,
            upstream,
            threads,
            batch_hits,
            max_unique,
            overlap_left,
            pam: pam,
            guide: guide,
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
            self.upstream,
            self.threads,
        )
        .map_err(|e| PyErr::new::<PyValueError, _>(e))?;

        debug_assert_eq!(pos_local.len(), strand.len());

        if cfg!(debug_assertions) {
            eprintln!(
                "[DEBUG] contig_id={} chunk_start={} size={} raw_hits={}",
                contig_id,
                chunk_start,
                self.size,
                pos_local.len()
            );
            for i in 0..pos_local.len().min(20) {
                eprintln!(
                    "  -> local_pos={} strand={}",
                    pos_local[i],
                    if strand[i] == 1 { '+' } else { '-' }
                );
            }
        }

        let chunk_len = seq_bitmask.len();
        if self.size == 0 || chunk_len < self.size {
            return Ok(FeedStatus {
                flushed: false,
                stats: BatcherStats {
                    hits_in_batch: self.hits_in_batch,
                    unique_windows: self.map.len(),
                },
            });
        }

        let max_start_excl = chunk_len - self.size + 1;
        let core_len = valid_len;

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
            let flushed = self.should_flush();
            return Ok(FeedStatus {
                flushed,
                stats: BatcherStats {
                    hits_in_batch: self.hits_in_batch,
                    unique_windows: self.map.len(),
                },
            });
        }

        for i in 0..pos_local.len() {
            let p = pos_local[i];
            if p < accept_lo || p >= accept_hi {
                continue;
            }

            let pos_global = chunk_start as usize + p;
            if pos_global > (u32::MAX as usize) {
                return Err(PyErr::new::<PyValueError, _>("Position overflow"));
            }

            let strand_bit = strand[i]; // 1=fwd (+), 0=rev (-)

            let start = p;
            let end = start + self.size;
            let window = &seq_bitmask[start..end];
            let key: WindowKey = window.to_vec().into_boxed_slice();

            let occ = pack_occ(contig_id, pos_global as u32, strand_bit);

            if let Some(v) = self.map.get_mut(&key) {
                v.push(occ);
            } else {
                self.map.insert(key, vec![occ]);
            }

            self.hits_in_batch += 1;
        }

        Ok(FeedStatus {
            flushed: self.should_flush(),
            stats: BatcherStats {
                hits_in_batch: self.hits_in_batch,
                unique_windows: self.map.len(),
            },
        })
    }

    pub fn flush_and_align(&mut self, max_mm: usize, bdna: usize, brna: usize) -> PyResult<()> {
        // Collect window batches on flush
        let batch: WindowBatch = self.flush_to_batch();

        println!("aligning");
        Ok(())
    }

    /// Flush remaining data at end of genome. Returns stats of what was flushed (and clears).
    pub fn finalize(&mut self) -> PyResult<BatcherStats> {
        let stats = BatcherStats {
            hits_in_batch: self.hits_in_batch,
            unique_windows: self.map.len(),
        };
        self.clear_batch();
        Ok(stats)
    }

    /// Introspection (debug)
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

    // TODO: Check if this is the best way to do it
    pub fn get_window_keys(&self) -> impl Iterator<Item = &WindowKey> {
        self.map.keys()
    }

    pub fn extract_alignment_rx(&mut self) -> Option<Receiver<AlignmentRingBatch>> {
        self.alignment_rx.take()
    }

    /// Convert the current batch (unique windows + occurrences) into a `WindowBatch`
    /// and clear internal state.
    pub fn flush_to_batch(&mut self) -> WindowBatch {
        let cap = self.max_unique; // invariant: max_unique <= miner src capacity

        // Fast path: whole map fits
        if self.map.len() <= cap {
            let map = std::mem::take(&mut self.map);
            let unique = map.len();
            let mut windows = Vec::with_capacity(unique);
            let mut occs = Vec::with_capacity(unique);
            let mut total_hits = 0usize;
            for (k, v) in map {
                total_hits += v.len();
                windows.push(k);
                occs.push(v);
            }
            self.hits_in_batch = 0;
            return WindowBatch {
                windows,
                occs,
                total_hits,
            };
        }

        // Overshoot path: emit exactly `cap` windows, keep the rest for the next submit
        let mut windows = Vec::with_capacity(cap);
        let mut occs = Vec::with_capacity(cap);
        let mut total_hits = 0usize;
        let take: Vec<WindowKey> = self.map.keys().take(cap).cloned().collect();
        for k in take {
            if let Some(v) = self.map.remove(&k) {
                total_hits += v.len();
                windows.push(k);
                occs.push(v);
            }
        }
        self.hits_in_batch -= total_hits; // retained windows stay counted
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

    pub fn get_sequence_len(&self) -> usize {
        self.size
    }
    pub fn get_guide(&self) -> Guide {
        self.guide.clone()
    }
}

/// WindowBatch
#[derive(Debug)]
pub struct WindowBatch {
    /// Unique windows, each length == sequence_len (aka `size` used in scanning/aligning)
    pub windows: Vec<WindowKey>,
    /// Occurrences for each window (parallel to `windows`)
    pub occs: Vec<Vec<Occ>>,
    /// Total occurrences across all windows
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
