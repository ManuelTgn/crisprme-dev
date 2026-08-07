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
use smallvec::SmallVec;

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

/// Occurrences of one unique window.
///
/// `SmallVec<[Occ; 2]>` occupies the same 24 bytes as `Vec<Occ>` — the inline
/// `[u64; 2]` fits in the space `Vec` uses for its pointer and capacity — so
/// this costs nothing in map-entry size. In a genome-wide 25-mer scan the
/// overwhelming majority of windows occur once or twice, and those now never
/// touch the heap: no `grow_one` on insert, and no `free` when `WindowBatch`
/// drops.
type OccList = SmallVec<[Occ; 2]>;

/// Result of scanning one physical orientation of a chunk.
///
/// The two variants exist because they map to different `FeedStatus.flushed`
/// values in the current public API: a chunk too short to hold a window reports
/// `false` unconditionally, whereas a chunk that was actually scanned reports
/// `should_flush()`. Collapsing them would be a silent behaviour change.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FeedOutcome {
    /// `chunk_len < size` — no window could be extracted, nothing was touched.
    TooShort,
    /// The accept window was established and any hits were recorded. Also
    /// covers the empty-accept-window case, which today reports `should_flush()`.
    Scanned,
}

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
    map: AHashMap<WindowKey, OccList>,
    hits_in_batch: usize,

    /// Reusable IUPAC-encoding buffer for `feed_chunk`. Held here so the
    /// per-chunk encode does not allocate and free a chunk-sized `Vec` on
    /// every call (twice per chunk, once per strand).
    scratch: Vec<u8>,
    /// Reverse-complement bitmask buffer, sibling to `scratch`. Held across
    /// chunks so the RC pass allocates only on the first chunk.
    scratch_rc: Vec<u8>,
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
            scratch: Vec::new(),
            scratch_rc: Vec::new(),
        })
    }

    pub fn feed_chunk(
        &mut self,
        py: Python<'_>,
        contig_id: u16,
        chunk_start: u32,
        strand: u8,
        chunk_seq: &str,
        valid_len: usize,
    ) -> PyResult<FeedStatus> {
        // `feed_bitmask` takes `&mut self`, so the scratch buffer is moved out
        // for the duration of the call and moved back afterwards — on the error
        // path too — to keep its allocation alive across chunks.
        let mut buf = std::mem::take(&mut self.scratch);
        iupac::sequence_encoder_into(chunk_seq, &mut buf);
        let outcome =
            py.detach(|| self.feed_bitmask(&buf, contig_id, chunk_start, strand, valid_len));
        self.scratch = buf;

        let flushed = match outcome? {
            FeedOutcome::TooShort => false,
            FeedOutcome::Scanned => self.should_flush(),
        };
        Ok(FeedStatus {
            flushed,
            stats: BatcherStats {
                hits_in_batch: self.hits_in_batch,
                unique_windows: self.map.len(),
            },
        })
    }

    /// Feed both physical orientations of one chunk from a single encode.
    ///
    /// `chunk_seq` is the forward chunk; the reverse-complement bitmask is
    /// derived from the encoded forward bitmask rather than from a second
    /// ASCII string, so the caller no longer builds one and the chunk is
    /// encoded once instead of twice.
    ///
    /// `strand_fwd` and `strand_rc` are the strand bits stamped on each
    /// orientation's occurrences. They are supplied by the caller rather than
    /// derived here because an upstream PAM swaps them — the RC pass is the
    /// one labelled strand 1, which forces the PAM downstream of the target.
    ///
    /// `chunk_start` is a forward contig coordinate and applies unchanged to
    /// both passes; `feed_bitmask` converts frames per hit.
    ///
    /// Unlike two separate `feed_chunk` calls, the flush signal is evaluated
    /// once, after both orientations. Batch boundaries therefore differ from
    /// the two-call form; the set of emitted rows does not.
    pub fn feed_chunk_both(
        &mut self,
        py: Python<'_>,
        contig_id: u16,
        chunk_start: u32,
        strand_fwd: u8,
        strand_rc: u8,
        chunk_seq: &str,
        valid_len: usize,
    ) -> PyResult<FeedStatus> {
        // `feed_bitmask` takes `&mut self`, so both buffers are moved out for
        // the duration and moved back afterwards — on the error path too — to
        // keep their allocations alive across chunks.
        let mut fwd = std::mem::take(&mut self.scratch);
        let mut rc = std::mem::take(&mut self.scratch_rc);

        // `chunk_seq` borrows Python-owned memory, so the encode must stay
        // inside the GIL — it is the only step that touches it. Everything
        // after this point is pure Rust over owned buffers.
        iupac::sequence_encoder_into(chunk_seq, &mut fwd);

        // Release the GIL for the RC build, both rayon scans, and the batching
        // loop. None of it touches the interpreter, and it is the bulk of the
        // call.
        let (res_fwd, res_rc) = py.detach(|| {
            iupac::revcomp_bitmask_into(&fwd, &mut rc);
            let a = self.feed_bitmask(&fwd, contig_id, chunk_start, strand_fwd, valid_len);
            let b = if a.is_ok() {
                self.feed_bitmask(&rc, contig_id, chunk_start, strand_rc, valid_len)
            } else {
                Ok(FeedOutcome::TooShort) // discarded; the forward error is returned
            };
            (a, b)
        });

        self.scratch = fwd;
        self.scratch_rc = rc;

        let outcome_fwd = res_fwd?;
        let outcome_rc = res_rc?;

        // Matches the two-call form: a chunk too short for a window reported
        // `false` unconditionally, anything else reported `should_flush()`.
        let flushed = if outcome_fwd == FeedOutcome::TooShort && outcome_rc == FeedOutcome::TooShort
        {
            false
        } else {
            self.should_flush()
        };
        Ok(FeedStatus {
            flushed,
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

    /// Scan one already-encoded orientation of a chunk and record its hits.
    ///
    /// `bitmask` is the IUPAC-encoded chunk in the orientation to scan, and
    /// `strand` is the strand bit stamped onto every occurrence it produces.
    /// `chunk_start` is always a *forward* contig coordinate regardless of
    /// which orientation `bitmask` holds; the frame conversion happens in
    /// `window_fwd_left` below.
    ///
    /// Split out of `feed_chunk` so that a chunk encoded once can be scanned in
    /// both orientations without re-crossing the FFI boundary.
    fn feed_bitmask(
        &mut self,
        bitmask: &[u8],
        contig_id: u16,
        chunk_start: u32,
        strand: u8,
        valid_len: usize,
    ) -> PyResult<FeedOutcome> {
        let pos_local = scanner::scan_targets_bitmask(
            bitmask,
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

        let chunk_len = bitmask.len();
        if self.size == 0 || chunk_len < self.size {
            return Ok(FeedOutcome::TooShort);
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
            return Ok(FeedOutcome::Scanned);
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
            let window = &bitmask[start..end];
            let key: WindowKey = pack_window(window);

            // Read candidate target PAM sequence
            let pstart = end - 1;
            let wpam = &bitmask[pstart..pstart + plen];
            let pam_id = self.pam.pam_index(wpam);

            let occ = pack_occ(contig_id, pam_id as u16, window_fwd_left as u32, strand);

            self.map.entry(key).or_default().push(occ);
            self.hits_in_batch += 1;
        }

        Ok(FeedOutcome::Scanned)
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
    pub occs: Vec<OccList>,
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
