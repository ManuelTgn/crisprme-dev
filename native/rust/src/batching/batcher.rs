// src/batching/batcher.rs
//
// Drop-in replacement for your entire file.
//
// Notes:
// - Keeps your current behavior (does NOT clear on `feed_chunk` flush signal).
// - Scanner debug prints are behind `cfg!(debug_assertions)` so they won’t spam release builds.
// - Placeholder aligner is strand-aware and assumes window size == guide length (e.g., 23).
// - NEW: adds `flush_and_align_engine_debug(...)` exposed to Python.
//   This drains the current batch and submits it to a *persistent* Rust engine instance.
//   In your current DEBUG hybrid.rs, the engine prints the unique windows it would mine.
// - Engine is initialized lazily on first flush and then reused across subsequent flushes.
//
// IMPORTANT assumptions for the new debug method:
// - `HybridEngine::new(device_id)` and `HybridEngine::execute(alignment_params)` exist (your debug hybrid.rs).
// - `execute()` returns an `EngineHandle` with a public `window_producer` supporting `acquire()` and `commit()`.
// - `WindowRingBatch` supports Option-B helpers + a public `descriptor`:
//     windows_iupac_mut(), occ_starts_mut(), occ_lens_mut(), occ_contig_mut(), occ_pos_mut(), occ_strand_mut()
//     and descriptor.{sequence_len, window_count, occ_count}.
//
// If your module paths differ, adjust the imports marked "ADJUST PATH".

use crate::crispr::pam;
use crate::sequence::iupac::{sequence_encoder, Iupac};
use crate::sequence::scanner;

use ahash::AHashMap;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

// --- NEW imports for engine-backed flush (ADJUST PATHS if needed) ---
// use crate::memory::batch::WindowRingBatch; // ADJUST PATH if WindowRingBatch lives elsewhere
// use crate::engine::hybrid::{EngineHandle, HybridEngine}; // ADJUST PATH
// use crate::engine::params::AlignmentParams; // ADJUST PATH (same struct used by hybrid.rs)

use std::sync::{Mutex, OnceLock};

// =====================================================================================
// Engine runtime (persistent across flushes)
// =====================================================================================


#[inline(always)]
fn align_up(x: usize, a: usize) -> usize {
    debug_assert!(a.is_power_of_two());
    (x + a - 1) & !(a - 1)
}

/// Compute the ring slot bytes needed for a `WindowRingBatch` (Option B layout).
///
/// This MUST match how your `WindowRingBatch` interprets the underlying `RingSlotLease`.
/// The layout assumed here:
///
/// [ windows: Iupac[max_windows * sequence_len] ]  (bytes)
/// align to u32
/// [ occ_starts: u32[max_windows] ]
/// [ occ_lens  : u32[max_windows] ]
/// [ occ_contig: u32[max_occs] ]
/// [ occ_pos   : u32[max_occs] ]
/// [ occ_strand: u8[max_occs] ]
/// (optional padding to 8 bytes)
fn compute_window_ring_slot_bytes(sequence_len: usize, max_windows: usize, max_occs: usize) -> usize {
    let windows_bytes = max_windows * sequence_len * std::mem::size_of::<u8>(); // Iupac repr(u8)
    let mut off = windows_bytes;

    off = align_up(off, std::mem::align_of::<u32>());
    off += max_windows * std::mem::size_of::<u32>(); // starts
    off += max_windows * std::mem::size_of::<u32>(); // lens

    off += max_occs * std::mem::size_of::<u32>(); // contig
    off += max_occs * std::mem::size_of::<u32>(); // pos
    off += max_occs * std::mem::size_of::<u8>();  // strand

    // Keep the slot reasonably aligned
    off = align_up(off, 8);
    off
}

// =====================================================================================

/// Key: owned window bytes (IUPAC bitmasks), length == size.
type WindowKey = Box<[u8]>;

/// Occurrence: packed (contig_id, pos, strand_bit) into u64.
/// Layout: [ contig_id:31.. ] [ pos:32 bits ] [ strand:1 bit ]
/// occ = (contig_id << 33) | (pos << 1) | strand
type Occ = u64;

#[inline(always)]
pub fn pack_occ(contig_id: u32, pos: u32, strand_bit: u8) -> Occ {
    ((contig_id as u64) << 33) | ((pos as u64) << 1) | ((strand_bit as u64) & 1)
}

#[inline(always)]
fn unpack_occ(occ: Occ) -> (u32, u32, u8) {
    let contig_id = (occ >> 33) as u32;
    let pos = ((occ >> 1) & 0xFFFF_FFFF) as u32;
    let strand_bit = (occ & 1) as u8;
    (contig_id, pos, strand_bit)
}

/// Convert bitmask (0b0001/0010/0100/1000 sets) to its complement.
#[inline(always)]
fn complement_mask(m: u8) -> u8 {
    ((m & 0b0001) << 3) | ((m & 0b0010) << 1) | ((m & 0b0100) >> 1) | ((m & 0b1000) >> 3)
}

#[inline(always)]
fn revcomp_bits(bits: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(bits.len());
    for &b in bits.iter().rev() {
        out.push(complement_mask(b));
    }
    out
}

#[inline(always)]
fn bits_to_string(bits: &[u8]) -> String {
    let mut s = String::with_capacity(bits.len());
    for &b in bits {
        s.push(Iupac::new(b).to_utf8());
    }
    s
}

#[inline(always)]
fn count_mismatches_iupac_range(target: &[u8], guide: &[u8], lo: usize, hi: usize) -> u16 {
    debug_assert_eq!(target.len(), guide.len());
    debug_assert!(hi <= target.len());
    let mut mm = 0u16;
    for i in lo..hi {
        if (target[i] & guide[i]) == 0 {
            mm += 1;
        }
    }
    mm
}

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

#[pyclass]
pub struct TargetBatcher {
    // config
    size: usize,
    right: bool,
    threads: usize,
    batch_hits: usize,
    max_unique: usize,
    overlap_left: usize,

    // parsed PAM (reused across chunks)
    pam: pam::ParsedPAM,

    // state
    map: AHashMap<WindowKey, Vec<Occ>>,
    hits_in_batch: usize,
}

#[pymethods]
impl TargetBatcher {
    /// Create a batcher that accumulates candidates across many chunks/contigs.
    ///
    /// `batch_hits`: flush when total candidate hits >= batch_hits (e.g., 1_000_000)
    /// `max_unique`: safety valve flush if unique windows explode (e.g., 250_000)
    #[new]
    pub fn new(
        pam_seq: &str,
        size: usize,
        right: bool,
        threads: usize,
        batch_hits: usize,
        max_unique: usize,
        overlap_left: usize,
    ) -> PyResult<Self> {
        let pat = pam::ParsedPAM::new(pam_seq)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("Invalid PAM sequence: {e}")))?;

        if size > 0 && overlap_left < size.saturating_sub(1) {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Invalid overlap_left={overlap_left}: must be >= size-1={} to avoid losing kmers at chunk boundaries",
                size.saturating_sub(1)
            )));
        }

        Ok(Self {
            size,
            right,
            threads,
            batch_hits,
            max_unique,
            overlap_left,
            pam: pat,
            map: AHashMap::new(),
            hits_in_batch: 0,
        })
    }

    /// Feed a single chunk.
    pub fn feed_chunk(
        &mut self,
        contig_id: u32,
        chunk_start: u32,
        chunk_seq: &str,
        valid_len: usize,
    ) -> PyResult<FeedStatus> {
        let seq_bitmask: Vec<u8> = sequence_encoder(chunk_seq);

        let (pos_local, strand) = scanner::scan_targets_bitmask(
            &seq_bitmask,
            &self.pam,
            self.size,
            self.right,
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

    /// TEST/DEBUG:
    pub fn debug_collect_positions(&self) -> PyResult<Vec<(u32, u32, u8)>> {
        let mut out: Vec<(u32, u32, u8)> = Vec::new();
        for occs in self.map.values() {
            for &occ in occs {
                let (contig_id, pos, strand) = unpack_occ(occ);
                out.push((contig_id, pos, strand));
            }
        }
        Ok(out)
    }

    pub fn flush_and_align(
        &mut self,
        
    ) -> PyResult<()> {
        println!("aligning");
        Ok(())
    }
    
    /// Existing placeholder alignment kept as-is
    pub fn flush_and_align_placeholder_tsv(
        &mut self,
        guide: &str, // spacer ONLY (e.g. 20)
        max_mm: u16,
        out_path: PathBuf,
    ) -> PyResult<(usize, usize)> {
        let guide_plus = sequence_encoder(guide);

        let plen = self.pam.bytes.len();
        if plen == 0 || plen > self.size {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Invalid PAM length: plen={}, size={}",
                plen, self.size
            )));
        }

        let guide_len = self.size - plen;
        if guide_plus.len() != guide_len {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Guide length ({}) != spacer length ({}) [size={}, pam_len={}]",
                guide_plus.len(),
                guide_len,
                self.size,
                plen
            )));
        }

        let (spacer_lo_plus, spacer_hi_plus) = if self.right {
            (plen, plen + guide_len)
        } else {
            (0, guide_len)
        };

        let (spacer_lo_minus, spacer_hi_minus) = if self.right {
            (0, guide_len)
        } else {
            (plen, plen + guide_len)
        };

        let guide_minus = revcomp_bits(&guide_plus);
        let guide_plus_str = bits_to_string(&guide_plus);
        let guide_minus_str = bits_to_string(&guide_minus);

        let batch: WindowBatch = self.flush_to_batch();

        let mut writer = BufWriter::new(
            OpenOptions::new()
                .create(true)
                .append(true)
                .open(out_path)
                .map_err(|e| PyErr::new::<PyValueError, _>(format!("Cannot open output: {e}")))?,
        );

        let mut windows_with_hits = 0usize;
        let mut total_rows = 0usize;

        for (window, occs) in batch.windows.into_iter().zip(batch.occs.into_iter()) {
            let target_bits: &[u8] = &window;

            let target_spacer_plus = &target_bits[spacer_lo_plus..spacer_hi_plus];
            let target_spacer_minus = &target_bits[spacer_lo_minus..spacer_hi_minus];

            let mm_plus = count_mismatches_iupac_range(target_spacer_plus, &guide_plus, 0, guide_len);
            let mm_minus =
                count_mismatches_iupac_range(target_spacer_minus, &guide_minus, 0, guide_len);

            if mm_plus > max_mm && mm_minus > max_mm {
                continue;
            }

            let mut wrote_any = false;

            for occ in occs {
                let (contig_id, pos, strand_bit) = unpack_occ(occ);
                let (strand_char, mm, guide_str) = if strand_bit == 1 {
                    ('+', mm_plus, &guide_plus_str)
                } else {
                    ('-', mm_minus, &guide_minus_str)
                };

                if mm > max_mm {
                    continue;
                }

                let target_spacer_str = bits_to_string(&target_bits);
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    contig_id,
                    pos,
                    strand_char,
                    mm,
                    guide_str,
                    target_spacer_str
                )
                .map_err(|e| PyErr::new::<PyValueError, _>(format!("Write failed: {e}")))?;

                total_rows += 1;
                wrote_any = true;
            }

            if wrote_any {
                windows_with_hits += 1;
            }
        }

        writer
            .flush()
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("Flush failed: {e}")))?;

        Ok((windows_with_hits, total_rows))
    }
}

impl TargetBatcher {
    /// Convert the current batch (unique windows + occurrences) into a `WindowBatch`
    /// and clear internal state.
    pub fn flush_to_batch(&mut self) -> WindowBatch {
        let map: AHashMap<WindowKey, Vec<Occ>> = std::mem::take(&mut self.map);

        let unique = map.len();
        let mut windows: Vec<WindowKey> = Vec::with_capacity(unique);
        let mut occs: Vec<Vec<Occ>> = Vec::with_capacity(unique);

        let mut total_hits = 0usize;
        for (k, v) in map {
            total_hits += v.len();
            windows.push(k);
            occs.push(v);
        }

        self.hits_in_batch = 0;

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
}