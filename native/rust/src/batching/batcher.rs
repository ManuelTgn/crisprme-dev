// src/batching/batcher.rs
//
// Drop-in replacement for your entire file.
//
// Notes:
// - Keeps your current behavior (does NOT clear on `feed_chunk` flush signal).
// - Scanner debug prints are behind `cfg!(debug_assertions)` so they won’t spam release builds.
// - Placeholder aligner is strand-aware and assumes window size == guide length (e.g., 23).
// - Mismatch counting excludes PAM bases (PAM already enforced by scanner).
// - Output is TSV: contig_id, pos, strand, mismatches, guide_aligned, target_aligned
//
// IMPORTANT: this file assumes:
// - `pam::ParsedPAM` exposes `.bytes` and `.revcomp` (as in your scanner).
// - `sequence_encoder` uses IUPAC bitmasks (0b1111 for N).
// - Reverse complement bitmasks are computed here via bit-level complement (no dependence on Iupac::complement()).

use crate::crispr::pam;
use crate::sequence::iupac::{sequence_encoder, Iupac};
use crate::sequence::scanner;

use ahash::AHashMap;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

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
fn unpack_occ(occ: Occ) -> (u32, u32, u8) {
    let contig_id = (occ >> 33) as u32;
    let pos = ((occ >> 1) & 0xFFFF_FFFF) as u32;
    let strand_bit = (occ & 1) as u8;
    (contig_id, pos, strand_bit)
}

/// Convert bitmask (0b0001/0010/0100/1000 sets) to its complement.
/// Works for ambiguity codes too by mapping each base bit independently.
#[inline(always)]
fn complement_mask(m: u8) -> u8 {
    // A(0001)->T(1000)
    // C(0010)->G(0100)
    // G(0100)->C(0010)
    // T(1000)->A(0001)
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

/// Determine the region of the window where we should count mismatches (exclude PAM).
///
/// Your scanner places PAM start as:
/// - fwd: pam_start_fwd = if right { 0 } else { size - plen }
/// so:
/// - right == false => PAM is on the RIGHT (end) of the window => compare [0 .. size-plen)
/// - right == true  => PAM is on the LEFT (start) of the window => compare [plen .. size)
#[inline(always)]
fn mismatch_range(size: usize, plen: usize, right: bool) -> (usize, usize) {
    debug_assert!(plen <= size);
    if right {
        (plen, size)
    } else {
        (0, size - plen)
    }
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
    ///
    /// Arguments:
    /// - contig_id: numeric contig ID managed by Python
    /// - chunk_start: start coordinate in contig (0-based)
    /// - chunk_seq: ASCII sequence for the chunk, including right overlap (size-1) if available
    /// - valid_len: length of the non-overlap "owned" region (usually 10Mb, shorter at contig end)
    ///
    /// Returns: FeedStatus(flushed, stats)
    pub fn feed_chunk(
        &mut self,
        contig_id: u32,
        chunk_start: u32,
        chunk_seq: &str,
        valid_len: usize,
    ) -> PyResult<FeedStatus> {
        // Convert chunk to IUPAC bitmask once (needed for window keys)
        let seq_bitmask: Vec<u8> = sequence_encoder(chunk_seq);

        // Run scanner on bitmask (no duplicate conversion)
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

        // Max start position (exclusive) such that window [p, p+size) fits.
        let max_start_excl = chunk_len - self.size + 1;
        let core_len = valid_len;

        // Overlap-aware acceptance interval in LOCAL coordinates.
        //
        // chunk_start == 0 => first chunk, no left-overlap conceptually:
        // accept only starts within [0, core_len)
        //
        // chunk_start != 0 => chunk has left overlap of overlap_left bases:
        // accept:
        //   A) core starts:     [overlap_left, overlap_left + core_len)
        //   B) recovery starts: [overlap_left - (size-1), overlap_left)
        //
        // Combined: [overlap_left-(size-1), overlap_left+core_len)
        let (accept_lo, mut accept_hi) = if chunk_start == 0 {
            (0usize, core_len)
        } else {
            let ov = self.overlap_left;
            let recovery = self.size.saturating_sub(1);
            let lo = ov.saturating_sub(recovery);
            let hi = ov + core_len;
            (lo, hi)
        };

        // Clamp hi to the range where windows fit
        if accept_hi > max_start_excl {
            accept_hi = max_start_excl;
        }

        // Empty acceptance range -> nothing to do
        if accept_hi <= accept_lo {
            let flushed = self.should_flush();
            // IMPORTANT: do NOT clear here; Python / flush_and_align will drain it.
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

            // Apply overlap-aware filter
            if p < accept_lo || p >= accept_hi {
                continue;
            }

            // Global contig coordinate: chunk_start corresponds to local position 0
            let pos_global = chunk_start as usize + p;
            if pos_global > (u32::MAX as usize) {
                return Err(PyErr::new::<PyValueError, _>("Position overflow"));
            }

            if cfg!(debug_assertions) {
                eprintln!(
                    "  [ACCEPTED] contig={} global_pos={} strand={}",
                    contig_id,
                    pos_global,
                    if strand[i] == 1 { '+' } else { '-' }
                );
            }

            let strand_bit = strand[i]; // 1=fwd (+), 0=rev (-)

            // Window key: own bytes once per unique window
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

        let flushed = self.should_flush();

        Ok(FeedStatus {
            flushed,
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
    /// Return all accepted (contig_id, pos, strand) triples currently stored in the batch.
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

    /// Placeholder alignment:
    /// - scanner works on full window size = spacer + PAM (e.g. 23)
    /// - BUT here we align only the spacer (e.g. 20)
    /// - '+' uses guide as-is; '-' uses revcomp(guide)
    /// - strand in TSV is the PAM strand (from scanner / Occ)
    ///
    /// Output TSV columns:
    /// contig_id  pos  strand  mismatches  guide_aligned  target_aligned
    pub fn flush_and_align_placeholder_tsv(
        &mut self,
        guide: &str,          // spacer ONLY (e.g. 20)
        max_mm: u16,
        out_path: PathBuf,
    ) -> PyResult<(usize, usize)> {
        // Encode spacer guide
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

        // Determine where the spacer lies inside the full window
        // right=false => PAM on RIGHT => spacer at [0..guide_len)
        // right=true  => PAM on LEFT  => spacer at [plen..plen+guide_len)
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

        // Precompute reverse-complement spacer guide bits (window is reference-forward)
        let guide_minus = revcomp_bits(&guide_plus);

        let guide_plus_str = bits_to_string(&guide_plus);
        let guide_minus_str = bits_to_string(&guide_minus);

        // Drain batch (consumes internal map)
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

            // Extract the spacer slice from the full candidate window
            let target_spacer_plus = &target_bits[spacer_lo_plus..spacer_hi_plus];
            let target_spacer_plus_str = bits_to_string(target_spacer_plus);

            let target_spacer_minus = &target_bits[spacer_lo_minus..spacer_hi_minus];
            let target_spacer_minus_str = bits_to_string(target_spacer_minus);

            // Mismatch counts (over spacer only)
            let mm_plus = count_mismatches_iupac_range(target_spacer_plus, &guide_plus, 0, guide_len);
            let mm_minus = count_mismatches_iupac_range(target_spacer_minus, &guide_minus, 0, guide_len);

            if mm_plus > max_mm && mm_minus > max_mm {
                continue;
            }

            let mut wrote_any = false;

            for occ in occs {
                let (contig_id, pos, strand_bit) = unpack_occ(occ);

                // PAM strand selects guide orientation; target stays reference-forward
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
    ///
    /// Intended to be called by `flush_and_align_*` methods.
    pub fn flush_to_batch(&mut self) -> WindowBatch {
        // Take ownership of the map without reallocating it entry-by-entry
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
        // self.map is already empty (mem::take)

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
