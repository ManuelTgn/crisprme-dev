use ahash::AHashMap;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::crispr::pam;
use crate::sequence::scanner;
use crate::utils::thresholds::Thresholds;
use crate::crispr::guide::Guide;
use crate::alignment::alignment::Alignment;
use crate::alignment::cigarx::CigarxOp;
use crate::bindings;
use crate::memory::batch::{SequenceRingBatch, AlignmentRingBatch};
use crate::memory::ring::ring_buffer;
use crate::storage::reader::SequenceBatchDescr;
use crate::sequence::iupac::Iupac;


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
fn occ_contig_id(occ: Occ) -> u32 { (occ >> 33) as u32 }
#[inline(always)]
fn occ_pos(occ: Occ) -> u32 { ((occ >> 1) & 0xFFFF_FFFF) as u32 }
#[inline(always)]
fn occ_strand(occ: Occ) -> u8 { (occ & 1) as u8 }

#[derive(Debug)]
pub struct WindowBatch {
    /// Unique windows, each length == sequence_len (aka `size` used in scanning/aligning)
    pub windows: Vec<WindowKey>,
    /// Occurrences for each window (parallel to `windows`)
    pub occs: Vec<Vec<Occ>>,
}

impl WindowBatch {
    #[inline] pub fn len(&self) -> usize { self.windows.len() }
    #[inline] pub fn is_empty(&self) -> bool { self.windows.is_empty() }
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
    // scanner config
    size: usize,
    right: bool,
    threads: usize,
    overlap_left: usize,

    // flush policy
    batch_hits: usize,
    max_unique: usize,

    pam: pam::ParsedPAM,
    map: AHashMap<WindowKey, Vec<Occ>>,
    hits_in_batch: usize,

    // contig names (index == contig_id)
    contigs: Vec<String>,

    // output dir + per-contig writers (append)
    out_dir: Option<PathBuf>,
    writers: Vec<Option<BufWriter<File>>>,

    // cuda init
    cuda_initialized: bool,
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
        let pam = pam::ParsedPAM::new(pam_seq)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("Invalid PAM: {e}")))?;

        // To avoid losing kmers across chunk boundaries, overlap_left must be >= size-1
        if overlap_left < size.saturating_sub(1) {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "overlap_left={overlap_left} must be >= size-1={} to avoid losing kmers at chunk boundaries",
                size.saturating_sub(1)
            )));
        }

        Ok(Self {
            size, right, threads, overlap_left,
            batch_hits, max_unique,
            pam,
            map: AHashMap::new(),
            hits_in_batch: 0,
            contigs: Vec::new(),
            out_dir: None,
            writers: Vec::new(),
            cuda_initialized: false,
        })
    }

    /// Must be called once before flush_and_align().
    /// The order defines contig_id (0..n-1).
    pub fn set_contigs(&mut self, contigs: Vec<String>) {
        self.contigs = contigs;
        self.writers = vec![None; self.contigs.len()];
    }

    pub fn feed_chunk(
        &mut self,
        contig_id: u32,
        chunk_start: u32,
        chunk_seq: &str,
        valid_len: usize,
    ) -> PyResult<FeedStatus> {
        // Fast lossy conversion to IUPAC bitmasks (u8 masks)
        // NOTE: invalid chars become N (0b1111)
        let seq_bitmask: Vec<u8> = sequence_encoder(chunk_seq);

        // Run scanner on the already-encoded bitmask
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
            return Ok(FeedStatus {
                flushed: false,
                stats: BatcherStats {
                    hits_in_batch: self.hits_in_batch,
                    unique_windows: self.map.len(),
                },
            });
        }

        // Max local start position (exclusive) such that window [p, p+size) fits.
        let max_start_excl = chunk_len - self.size + 1;
        let core_len = valid_len;

        // Overlap-aware acceptance interval in LOCAL coordinates.
        let (accept_lo, mut accept_hi) = if chunk_start == 0 {
            (0usize, core_len)
        } else {
            let ov = self.overlap_left;
            let recovery = self.size.saturating_sub(1);
            let lo = ov.saturating_sub(recovery);
            let hi = ov + core_len;
            (lo, hi)
        };

        // Clamp hi to where windows still fit
        if accept_hi > max_start_excl {
            accept_hi = max_start_excl;
        }

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

        // Insert occurrences into the map
        for i in 0..pos_local.len() {
            let p = pos_local[i];

            // overlap-aware filter
            if p < accept_lo || p >= accept_hi {
                continue;
            }

            // global pos
            let pos_global = chunk_start as usize + p;
            if pos_global > (u32::MAX as usize) {
                return Err(PyErr::new::<PyValueError, _>("Position overflow"));
            }

            let strand_bit = strand[i]; // 1=fwd, 0=rev

            // Build window key (owned boxed slice)
            let start = p;
            let end = start + self.size;
            let key: WindowKey = seq_bitmask[start..end].to_vec().into_boxed_slice();

            let occ = pack_occ(contig_id, pos_global as u32, strand_bit);

            // single hash lookup
            self.map.entry(key).or_insert_with(Vec::new).push(occ);
            self.hits_in_batch += 1;
        }

        et flushed = self.should_flush();
        if flushed {
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

    /// Flush current batch and align, writing per-contig TSV files (append).
    ///
    /// Returns number of TSV rows written.
    pub fn flush_and_align(
        &mut self,
        output_dir: &str,
        guide: &str,
        mism: u32,
        qgap: u32,
        tgap: u32,
        device: u32,
    ) -> PyResult<usize> {
        if self.contigs.is_empty() {
            return Err(PyErr::new::<PyValueError, _>(
                "TargetBatcher.set_contigs(contigs) must be called before flush_and_align()",
            ));
        }

        let wb = self.flush_to_batch();
        if wb.is_empty() {
            return Ok(0);
        }

        // init output dir + ensure writers vector matches contigs
        self.ensure_output_dir(output_dir)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("output dir error: {e}")))?;

        // init cuda miner once
        if !self.cuda_initialized {
            bindings::miner::initialize(device);
            self.cuda_initialized = true;
        }

        let thresholds = Thresholds { mism, qgap, tgap };
        let guide_fwd: Guide = guide.into();
        let guide_rev = guide_fwd.reverse_complement();

        // Prepare GPU batches (single-slot ring buffers; minimal overhead)
        let n = wb.len();

        // Sequence slot bytes: n*(sequence_len + 4) matches your engine heuristic
        let seq_slot_bytes = n * (self.size + 4);
        let (seq_prod, _seq_cons) = ring_buffer::<SequenceRingBatch>(1, seq_slot_bytes, true);
        let mut seq_batch = seq_prod.acquire();
        seq_batch.descriptor = SequenceBatchDescr {
            sequence_count: n,
            sequence_len: self.size,
            global_offset: 0,
        };

        fill_sequence_batch_from_windows(&mut seq_batch, &wb.windows, self.size)?;

        // Copy sequences to GPU
        seq_batch.as_mut().sync_cpu_to_gpu(None);

        // Alignment slot: choose a big capacity to reduce miner loop iterations
        let alig_slot_bytes = 1_000_000usize * std::mem::size_of::<Alignment>();
        let (alig_prod, _alig_cons) = ring_buffer::<AlignmentRingBatch>(1, alig_slot_bytes, true);
        let mut alig_batch = alig_prod.acquire();

        // Mine positive (+)
        bindings::miner::pre_mine(&guide_fwd, self.size, &thresholds, b'+');
        let mut written = mine_expand_write_tsv(
            &self.contigs,
            &mut self.writers,
            &self.out_dir.as_ref().unwrap(),
            &seq_batch,
            &mut alig_batch,
            &wb,
            &guide_fwd,
        )?;
        bindings::miner::post_mine();

        // Mine negative (-)
        bindings::miner::pre_mine(&guide_rev, self.size, &thresholds, b'-');
        written += mine_expand_write_tsv(
            &self.contigs,
            &mut self.writers,
            &self.out_dir.as_ref().unwrap(),
            &seq_batch,
            &mut alig_batch,
            &wb,
            &guide_rev,
        )?;
        bindings::miner::post_mine();

        // Flush all open writers (optional but nice)
        for w in self.writers.iter_mut().flatten() {
            w.flush().map_err(|e| PyErr::new::<PyValueError, _>(format!("flush failed: {e}")))?;
        }

        Ok(written)
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
    fn ensure_output_dir(&mut self, output_dir: &str) -> std::io::Result<()> {
        let p = PathBuf::from(output_dir);
        std::fs::create_dir_all(&p)?;
        self.out_dir = Some(p);
        if self.writers.len() != self.contigs.len() {
            self.writers = vec![None; self.contigs.len()];
        }
        Ok(())
    }
}

impl TargetBatcher {
    // Convert the current batch (unique windows + occurrences) into a `WindowBatch`
    /// and clear internal state.
    ///
    /// Intended to be called by a future `flush_and_align()` method.
    pub fn flush_to_batch(&mut self) -> WindowBatch {
        // Take ownership of the map without reallocating it entry-by-entry
        let map: AHashMap<WindowKey, Vec<Occ>> = std::mem::take(&mut self.map);
        // avoid cloning keys/values in map

        let unique = map.len();  // unique windows 
        let mut windows: Vec<WindowKey> = Vec::with_capacity(unique);
        let mut occs: Vec<Vec<Occ>> = Vec::with_capacity(unique);

        let mut total_hits = 0usize;

        for (k, v) in map {
            total_hits += v.len();
            windows.push(k);
            occs.push(v);
        }

        // Reset counters/state
        self.hits_in_batch = 0;
        // self.map is already empty beacuse we used mem::take

        WindowBatch { windows, occs, total_hits }  // return window batch
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

fn fill_sequence_batch_from_windows(
    batch: &mut SequenceRingBatch,
    windows: &[WindowKey],
    slen: usize,
) -> PyResult<()> {
    let n = windows.len();
    if batch.descriptor.sequence_count != n || batch.descriptor.sequence_len != slen {
        return Err(PyErr::new::<PyValueError, _>("Sequence batch descriptor mismatch"));
    }

    // Write IUPAC bytes (Iupac is repr(u8) in your utils)
    let dst_iupac = batch.iupac_mut();
    let dst_bytes: &mut [u8] = unsafe {
        std::slice::from_raw_parts_mut(dst_iupac.as_mut_ptr() as *mut u8, n * slen)
    };

    for (i, w) in windows.iter().enumerate() {
        if w.len() != slen {
            return Err(PyErr::new::<PyValueError, _>("Window length mismatch"));
        }
        let off = i * slen;
        dst_bytes[off..off + slen].copy_from_slice(w.as_ref());
    }

    // IDs are indices into windows/occs
    let ids = batch.ids_mut();
    for i in 0..n {
        ids[i] = i as u32;
    }

    Ok(())
}

/// Open (or reuse) a per-contig TSV writer
fn ensure_contig_writer(
    contigs: &[String],
    writers: &mut [Option<BufWriter<File>>],
    out_dir: &Path,
    contig_id: usize,
) -> PyResult<&mut BufWriter<File>> {
    if contig_id >= contigs.len() || contig_id >= writers.len() {
        return Err(PyErr::new::<PyValueError, _>("contig_id out of range"));
    }
    if writers[contig_id].is_none() {
        let path = out_dir.join(format!("{}.tsv", contigs[contig_id]));
        let f = OpenOptions::new().create(true).append(true).open(&path)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("cannot open {:?}: {e}", path)))?;
        let mut w = BufWriter::new(f);

        // Optional header if file is new-ish (best effort: check metadata len)
        // If you want always header, remove this check.
        // NOTE: in append mode, checking file size requires re-open; skip for throughput.
        // We'll just write header once per process if you want; currently omitted.

        writers[contig_id] = Some(w);
    }
    Ok(writers[contig_id].as_mut().unwrap())
}

/// Mine (single strand pass), expand to occurrences, and append TSV rows.
/// Returns number of TSV rows written.
fn mine_expand_write_tsv(
    contigs: &[String],
    writers: &mut [Option<BufWriter<File>>],
    out_dir: &Path,
    seq_batch: &SequenceRingBatch,
    alig_batch: &mut AlignmentRingBatch,
    wb: &WindowBatch,
    guide_for_this_pass: &Guide,
) -> PyResult<usize> {
    let mut written = 0usize;

    loop {
        let complete = bindings::miner::mine(seq_batch, alig_batch);

        let acount = alig_batch.len();
        alig_batch.sync_gpu_to_cpu(Some(std::mem::size_of::<Alignment>() * acount));

        for al in alig_batch.alignments().iter() {
            let win_idx = al.id as usize;
            if win_idx >= wb.occs.len() {
                return Err(PyErr::new::<PyValueError, _>("alignment id out of range"));
            }

            let occ_list = &wb.occs[win_idx];
            let window_bytes = wb.windows[win_idx].as_ref();

            // Convert window bytes to IUPAC slice for rendering
            let window_iupac: &[Iupac] = unsafe {
                std::slice::from_raw_parts(window_bytes.as_ptr() as *const Iupac, window_bytes.len())
            };

            for &occ in occ_list.iter() {
                let cid = occ_contig_id(occ) as usize;
                let wstart = occ_pos(occ);
                let _cand_strand = occ_strand(occ); // available if you want it later

                let chrom = contigs.get(cid)
                    .ok_or_else(|| PyErr::new::<PyValueError, _>("contig id not in contigs list"))?;

                // POS = window_start + alignment.offset (as you did in CLI Alignments)
                let pos = wstart + (al.offset as u32);

                // Build aligned strings + counts
                let (g_aln, t_aln, mm, bdna, brna) = build_alignment_strings_and_counts(
                    guide_for_this_pass,
                    window_iupac,
                    al,
                );

                let strand_char = al.strand as char;

                let w = ensure_contig_writer(contigs, writers, out_dir, cid)?;
                // TSV row
                // chrom pos strand guide_aligned target_aligned mm bdna brna
                writeln!(
                    w,
                    "{chrom}\t{pos}\t{strand_char}\t{g_aln}\t{t_aln}\t{mm}\t{bdna}\t{brna}"
                ).map_err(|e| PyErr::new::<PyValueError, _>(format!("write failed: {e}")))?;

                written += 1;
            }
        }

        // reset batch length so next mine() overwrites cleanly
        alig_batch.set_len(0);

        if complete {
            break;
        }
    }

    Ok(written)
}

/// Build full-window aligned strings (guide padded with '-' on flanks) and compute mm/bdna/brna.
///
/// Counts mapping (consistent with your Results printing):
/// - Mismatch -> mm
/// - Insertion (gap in guide) -> bdna
/// - Deletion (gap in target) -> brna
fn build_alignment_strings_and_counts(
    guide: &Guide,
    target_window: &[Iupac],
    al: &Alignment,
) -> (String, String, u32, u32, u32) {
    let mut gline = String::new();
    let mut tline = String::new();

    let mut qidx: usize = 0;
    let mut tidx: usize = 0;

    let mut mm: u32 = 0;
    let mut bdna: u32 = 0;
    let mut brna: u32 = 0;

    // prefix before offset: keep target bases, pad guide with '-'
    for _ in 0..(al.offset as usize) {
        if tidx < target_window.len() {
            tline.push(target_window[tidx].to_utf8());
            gline.push('-');
            tidx += 1;
        }
    }

    for op in al.cigarx.operations() {
        match op {
            CigarxOp::Match => {
                tline.push(target_window[tidx].to_utf8());
                gline.push(guide[qidx].to_utf8());
                qidx += 1;
                tidx += 1;
            }
            CigarxOp::Mismatch => {
                mm += 1;
                tline.push(target_window[tidx].to_utf8());
                gline.push(guide[qidx].to_utf8());
                qidx += 1;
                tidx += 1;
            }
            CigarxOp::Deletion => {
                // gap in target -> RNA bulge (extra base in guide)
                brna += 1;
                tline.push('-');
                gline.push(guide[qidx].to_utf8());
                qidx += 1;
            }
            CigarxOp::Insertion => {
                // gap in guide -> DNA bulge (extra base in target)
                bdna += 1;
                tline.push(target_window[tidx].to_utf8());
                gline.push('-');
                tidx += 1;
            }
        }
    }

    // suffix: remaining target bases, pad guide with '-'
    while tidx < target_window.len() {
        tline.push(target_window[tidx].to_utf8());
        gline.push('-');
        tidx += 1;
    }

    (gline, tline, mm, bdna, brna)
}
