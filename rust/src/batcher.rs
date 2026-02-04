
use crate::pam;
use crate::iupac::Iupac;
use crate::scanner;


use ahash::AHashMap;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

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


#[pyclass]
#[derive(Clone)]
pub struct BatcherStats {
    #[pyo3(get)]
    pub hits_in_batch: usize,
    #[pyo3(get)]
    pub unique_windows: usize,
}


/// Simple status enum for Python.
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
        let seq_bitmask: Vec<u8> = chunk_seq
            .as_bytes()
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                Iupac::from_ascii(b)
                    .map(|iupac| iupac.0)
                    .map_err(|e| format!("Invalid base at chunk offset {i}: {e}"))
            })
            .collect::<Result<_, _>>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e))?;

        // Run scanner on bitmask (no duplicate conversion)
        let (pos_local, strand) = scanner::scan_targets_bitmask(
            &seq_bitmask,
            &self.pam,
            self.size,
            self.right,
            self.threads,
        ).map_err(|e| PyErr::new::<PyValueError, _>(e))?;

        // Overlap-aware filtering to avoid duplicates *and* avoid losing boundary windows.
        //
        // Python chunk() uses LEFT overlap only:
        // - chunk 0: ext_start=0, chunk_start=0, no left overlap
        // - chunk i>0: ext_start=core_start - overlap_left, chunk_start=ext_start
        //
        // Let core_len = valid_len (passed from Python).
        // In local chunk coordinates:
        // - overlap region is [0, overlap_left)
        // - core region begins at overlap_left
        //
        // We must accept:
        //   A) core starts:      p in [overlap_left, overlap_left + core_len)
        //   B) recovery starts:  p in [overlap_left - (size-1), overlap_left)
        //      (these are the windows that couldn't fit at the end of previous core chunk)
        //
        // For chunk 0 (chunk_start==0): accept only core-like [0, core_len)
        // because there is no previous chunk to recover.
        debug_assert_eq!(pos_local.len(), strand.len());

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
            if flushed {
                self.clear_batch();
            }
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

            let strand_bit = strand[i]; // 1=fwd, 0=rev

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
        if flushed {
            // Next step: real flush to aligner. For now, clear.
            self.clear_batch();
        }

        Ok(FeedStatus {
            flushed,
            stats: BatcherStats {
                hits_in_batch: self.hits_in_batch,
                unique_windows: self.map.len(),
            },
        })
    }

    /// Flush remaining data at end of genome. Returns stats of what was flushed.
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
                let contig_id = (occ >> 33) as u32;
                let pos = ((occ >> 1) & 0xFFFF_FFFF) as u32;
                let strand = (occ & 1) as u8;
                out.push((contig_id, pos, strand));
            }
        }

        Ok(out)
    }
}


impl TargetBatcher {
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



