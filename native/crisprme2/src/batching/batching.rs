use ahash::AHashMap;
use crossbeam_channel::Receiver;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::crispr::guide::Guide;
use crate::crispr::{guide, pam};
use crate::memory::batch::AlignmentRingBatch;
use crate::model::occurence::Strand;
use crate::sequence::iupac::Iupac;
use crate::sequence::{iupac, scanner};

/// Key: the window's IUPAC bitmask bytes, right-padded with zeros to
/// `WINDOW_MAX_BASES`
pub type WindowKey = [u8; WINDOW_MAX_BASES];

/// Bases per window that fit in a `WindowKey`. Must equal `SEQ_MAX_LEN`
pub const WINDOW_MAX_BASES: usize = 32;

#[inline(always)]
fn pack_window(window: &[u8]) -> WindowKey {
    debug_assert!(window.len() <= WINDOW_MAX_BASES);
    let mut key = [0u8; WINDOW_MAX_BASES];
    key[..window.len()].copy_from_slice(window);
    key
}

/// Occurrence: packed (contig_id, pos, strand_bit) into u64.
/// Layout: [ contig_id:31.. ] [ pos:32 bits ] [ strand:1 bit ]
/// occ = (contig_id << 33) | (pos << 1) | strand
type Occ = u64;

#[inline(always)]
fn pack_occ(contig_id: u16, pam_id: u16, pos: u32, strand_bit: u8) -> Occ {
    ((contig_id as u64) << 49)
        | ((pam_id as u64) << 33)
        | ((pos as u64) << 1)
        | ((strand_bit as u64) & 1)
}

#[inline(always)]
pub fn unpack_occ(occ: Occ) -> (u16, u16, u32, u8) {
    let contig_id = (occ >> 49) as u16;
    let pam_id = (occ >> 33) as u16;
    let pos = ((occ >> 1) & 0xFFFF_FFFF) as u32;
    let strand_bit = (occ & 1) as u8;
    (contig_id, pam_id, pos, strand_bit)
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
                "Invalid overlap_left={overlap_left}: must be >= size-1={} to avoid \
                 losing kmers at chunk boundaries",
                size.saturating_sub(1)
            )));
        }

        if size > WINDOW_MAX_BASES {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Window size {size} exceeds the maximum of {WINDOW_MAX_BASES} bases \
                 (guide + PAM + max bulge must fit in a WindowKey)"
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
            map: AHashMap::with_capacity(max_unique),
            hits_in_batch: 0,
        })
    }

    pub fn feed_chunk(
        &mut self,
        contig_id: u16,
        chunk_start: u32,
        strand: u8,
        chunk_seq: &str,
        valid_len: usize,
    ) -> PyResult<FeedStatus> {
        let seq_bitmask: Vec<u8> = iupac::sequence_encoder(chunk_seq);

        let pos_local = scanner::scan_targets_bitmask(
            &seq_bitmask,
            &self.pam,
            self.size,
            self.upstream,
            self.threads,
        )
        .map_err(|e| PyErr::new::<PyValueError, _>(e))?;

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
                    if strand == 1 { '+' } else { '-' }
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

        // Per-chunk, not per-hit: which physical orientation did the scanner see?
        let scanned_on_rc = Strand::from_bit(strand).scanned_on_revcomp(self.upstream);
        let plen = self.pam.bytes.len();

        for i in 0..pos_local.len() {
            let p = pos_local[i];
            if p < accept_lo || p >= accept_hi {
                continue;
            }

            let start = p;
            let end = start + (self.size - plen) + 1;

            // Left-most FORWARD coordinate of this window.
            //
            // `chunk_start` is a forward coordinate; `p` indexes the *scanned*
            // sequence. On a forward chunk the two frames agree and they simply
            // add. On an RC chunk the scanned frame runs against the forward
            // strand, so index `p` sits `chunk_len - p - size` bases from the
            // chunk's forward start — the old `chunk_start + p` mixed the two
            // frames.
            //
            //   forward chunk: chunk_start + p
            //   RC chunk:      chunk_start + chunk_len - p - size
            let window_fwd_left = if scanned_on_rc {
                // `end <= chunk_len` holds because `p < max_start_excl`;
                // checked anyway so a future change to the accept-window can't
                // silently wrap
                let back = chunk_len.checked_sub(end).ok_or_else(|| {
                    PyErr::new::<PyValueError, _>(format!(
                        "window [{start},{end}) escapes chunk (len={chunk_len})"
                    ))
                })?;
                chunk_start as usize + back - plen + 1
            } else {
                chunk_start as usize + start
            };

            if window_fwd_left > u32::MAX as usize {
                return Err(PyErr::new::<PyValueError, _>("Position overflow"));
            }

            // Read candidate target sequence
            let window = &seq_bitmask[start..end];
            let key: WindowKey = pack_window(window);

            // Read candidate target PAM sequence
            let pstart = end - 1;
            let wpam = &seq_bitmask[pstart..pstart + plen];
            let pam_id = self.pam.pam_index(wpam);

            let occ = pack_occ(contig_id, pam_id as u16, window_fwd_left as u32, strand);

            self.map.entry(key).or_default().push(occ);
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

    pub fn get_window_keys(&self) -> impl Iterator<Item = WindowKey> + '_ {
        self.map.keys().copied()
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
            let unique = self.map.len();
            let mut windows = Vec::with_capacity(unique);
            let mut occs = Vec::with_capacity(unique);
            let mut total_hits = 0usize;
            // drain() empties the map but keeps its table allocated, so the
            // next batch does not rehash from zero capacity
            for (k, v) in self.map.drain() {
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
        let take: Vec<WindowKey> = self.map.keys().take(cap).copied().collect();
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
    pub fn get_pam_len(&self) -> usize {
        self.pam.plen()
    }
    pub fn get_guide(&self) -> Guide {
        self.guide.clone()
    }
}

/// WindowBatch
#[derive(Debug)]
pub struct WindowBatch {
    /// Unique windows as packed `u128` keys, in emission order. Expand with
    /// `unpack_window(key, size, out)`; `size` comes from
    /// `TargetBatcher::get_sequence_len`
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
