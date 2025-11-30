use std::path::PathBuf;

use crate::utils::{iupac_match, Thresholds, IUPAC};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use num_traits::{ops::bytes, PrimInt, Unsigned, Zero};
use rayon::prelude::*;

/// Metric collected by a miner during execution
pub struct Metrics {
    pub pruned: usize,
    pub mined: usize,
}

impl Metrics {
    /// Single valid result
    pub fn mined() -> Self {
        Self {
            pruned: 0,
            mined: 1,
        }
    }

    /// Single pruned result
    pub fn pruned() -> Self {
        Self {
            pruned: 1,
            mined: 0,
        }
    }

    /// Empty metrics
    pub fn empty() -> Self {
        Self {
            pruned: 0,
            mined: 0,
        }
    }

    /// Combine two metrics
    pub fn aggregate(&mut self, other: Metrics) {
        self.pruned += other.pruned;
        self.mined += other.mined;
    }
}

/// Miner DFS state
#[derive(Debug)]
struct State {
    qgaps: u32,
    tgaps: u32,
    misms: u32,
}

impl State {
    pub fn ed(&self) -> u32 {
        self.qgaps + self.tgaps + self.misms
    }

    pub fn invalid(&self, config: &Thresholds) -> bool {
        self.qgaps > config.qgap
            || self.tgaps > config.tgap
            || self.misms > config.mism
    }
}

/// Struct responsable for mining the alignments
#[derive(Debug, Clone)]
pub struct Miner<T> {
    /// Cache for beamlight search, allocated only once
    cache_min_ed: Vec<u32>,
    /// Current cigarx
    cigarx: CIGARX<T>,
    /// Current solutions
    solutions: Vec<Alignment<T>>,
    /// Thresholds for gaps and mismatches
    thresholds: Thresholds,
    /// Len of the guide
    guide_len: usize,
    /// Len of the target sequences
    seq_len: usize,
    /// Current sequence id
    seq_id: u32,
    /// Current starting offset
    start_offset: u8,
    /// Current guide
    q: Vec<IUPAC>,
    /// Current target
    t: Vec<IUPAC>
}

impl<T: PrimInt + Unsigned + Zero> Miner<T> {
    pub fn new(thresholds: &Thresholds, guide_len: usize, seq_len: usize) -> Self {
        Self {
            cache_min_ed: vec![0; (guide_len + 1) * (seq_len + 1)],
            guide_len,
            seq_len,
            cigarx: CIGARX(T::zero(), 0),
            solutions: vec![],
            thresholds: *thresholds,
            seq_id: 0,
            start_offset: 0,
            q: vec![],
            t: vec![],
        }
    }

    /// Recursive DFS to explore the alignments' state space
    fn mine_dfs(&mut self, qi: usize, ti: usize, state: State) -> Metrics {

        // Solution
        if qi == self.q.len() && !state.invalid(&self.thresholds) {
            self.solutions.push(Alignment { 
                cigarx: self.cigarx, 
                offset: self.start_offset, 
                seq_id: self.seq_id 
            });

            return Metrics::mined();
        }

        // Prune if minimum edit distance to solution is incompatible
        if state.ed() + self.min_ed(qi, ti) > self.max_ed() {
            return Metrics::pruned();
        }

        // Result metrics
        let mut metrics = Metrics {
            pruned: 0,
            mined: 0,
        };

        // Match/mismatch
        if qi < self.guide_len && ti < self.seq_len {
            let delta = if iupac_match(self.q[qi], self.t[ti]) {
                0
            } else {
                1
            };

            cigarx_push(&mut self.cigarx, if delta == 0 { b'=' } else { b'X' });        
            metrics.aggregate(self.mine_dfs(
                qi + 1,
                ti + 1,
                State {
                    qgaps: state.qgaps,
                    tgaps: state.tgaps,
                    misms: state.misms + delta,
                },
            ));
            cigarx_pop(&mut self.cigarx);
        }

        // Target gap
        if qi < self.guide_len {
            cigarx_push(&mut self.cigarx, b'I');
            metrics.aggregate(self.mine_dfs(
                qi + 1,
                ti,
                State {
                    qgaps: state.qgaps,
                    tgaps: state.tgaps + 1,
                    misms: state.misms,
                },
            ));
            cigarx_pop(&mut self.cigarx);
        }

        // Query gap
        if ti < self.seq_len {
            cigarx_push(&mut self.cigarx, b'D');
            metrics.aggregate(self.mine_dfs(
                qi,
                ti + 1,
                State {
                    qgaps: state.qgaps + 1,
                    tgaps: state.tgaps,
                    misms: state.misms,
                },
            ));
            cigarx_pop(&mut self.cigarx);
        }

        metrics

    }

    pub fn prepare(&mut self, q: &[IUPAC], t: &[IUPAC], id: u32) {
        self.cigarx = CIGARX(T::zero(), 0);
        self.solutions.clear();
        self.fill_min_ed(q, t);

        // NOTE: We should copy into our buffer instead of creating a new vector!
        self.q = Vec::from(q);
        self.t = Vec::from(t);

        self.start_offset = 0;
        self.seq_id = id;
    }

    /// Mine all valid alignments
    pub fn mine(&mut self, q: &[IUPAC], t: &[IUPAC], seq_id: u32) -> Metrics {
        
        self.prepare(q, t, seq_id);

        let mut metrics = Metrics::empty();
        for ti in 0..self.seq_len {
            self.start_offset = ti as u8;
            metrics.aggregate(self.mine_dfs(
                0,
                ti,
                State {
                    qgaps: 0,
                    tgaps: 0,
                    misms: 0,
                },
            ));
        }

        metrics
    }

    pub fn solutions(&self) -> &[Alignment<T>] {
        &self.solutions
    }

    /// Maximum pruning edit distance
    pub fn max_ed(&self) -> u32 {
        self.thresholds.ed()
    }

    /// Recompute the minimum edit distance cache
    fn fill_min_ed(&mut self, q: &[IUPAC], t: &[IUPAC]) {
        // Solution, all of the bottom row is a solution
        for ti in 0..=self.seq_len {
            *self.min_ed_mut(self.guide_len, ti) = 0;
        }

        // Insert remaining T as gaps
        for qi in 0..=self.guide_len {
            let value = (self.guide_len - qi) as u32;
            *self.min_ed_mut(qi, self.seq_len) = value;
        }

        for qi in (0..self.guide_len).rev() {
            for ti in (0..self.seq_len).rev() {
                // NOTE: This also handles wildcards
                let m = if iupac_match(q[qi], t[ti]) {
                    0
                } else {
                    1
                };

                let a = self.min_ed(qi + 1, ti + 1) + m;
                let b = self.min_ed(qi + 1, ti + 0) + 1;
                let c = self.min_ed(qi + 0, ti + 1) + 1;
                let value = a.min(b).min(c);

                *self.min_ed_mut(qi, ti) = value;
            }
        }

        /*
        // Print matrix
        for qi in 0..=self.config.query_len {
            println!();
            for ti in 0..=self.config.target_len {
                print!("{:>3}", self.min_ed(qi, ti));
            }
        }
        */

    }

    /// Return the minimum edit distance to a solution
    fn min_ed(&self, q: usize, t: usize) -> u32 {
        self.cache_min_ed[q * (self.seq_len + 1) + t]
    }

    fn min_ed_mut(&mut self, q: usize, t: usize) -> &mut u32 {
        &mut self.cache_min_ed[q * (self.seq_len + 1) + t]
    }
}



// CIGARX defined over a generic integer type
#[derive(Debug, Clone, Copy)]
pub struct CIGARX<T>(T, u8);

/// Encode a CIGARX value into an u8
pub fn cigarx_encode_single(value: u8) -> u8 {
    match value {

        b'=' => 0b00,
        b'I' => 0b01,
        b'D' => 0b10,
        b'X' => 0b11,

        _ => unimplemented!()
    }
}

/// Decode a single CIGARX value
pub fn cigarx_decode_single(value: u8) -> u8 {
    match value {

        0b00 => b'=',
        0b01 => b'I',
        0b10 => b'D',
        0b11 => b'X',

        _ => unimplemented!()
    }
}

pub fn cigarx_decode<T>(mut cigarx: CIGARX<T>) -> String 
where
    T: PrimInt + Unsigned + Zero 
{
    let mut result: Vec<u8> = vec![]; 
    for _ in 0..cigarx.1 {
        let value: T = cigarx.0 & T::from(0b11).unwrap();
        result.push(cigarx_decode_single(value.to_u8().unwrap()));
        cigarx.0 = cigarx.0 >> 2;
    }

    result.reverse();
    String::from(unsafe { str::from_utf8_unchecked(&result) })
}

/// Pack a CIGARX element into an integer type
pub fn cigarx_push<T>(current: &mut CIGARX<T>, value: u8)
where
    T: PrimInt + Unsigned + Zero 
{
    current.0 = current.0 << 2;
    current.0 = current.0 | T::from(cigarx_encode_single(value)).unwrap();
    current.1 += 1;
}

/// Pop last CIGARX element from the integer type
pub fn cigarx_pop<T>(current: &mut CIGARX<T>)
where 
    T: PrimInt + Unsigned + Zero 
{
    assert_ne!(current.1, 0);
    current.1 -= 1;
    current.0 = current.0 >> 2;
}  

#[derive(Debug, Clone, Copy)]
pub struct Alignment<T> {
    pub cigarx: CIGARX<T>,
    pub seq_id: u32,
    pub offset: u8,
}

#[derive(Debug)]
pub struct MinedGenome<T> {
    pub alignments: Vec<Alignment<T>>
}

const MINER_CPU_MULT: usize = 64;

/*
/// Implementation of a MinedGenome with a valid integer type
impl<T: PrimInt + Unsigned + Zero + Send> MinedGenome<T> {
    pub fn from_filtered_genome(genome: &FilteredGenome, guide: &[IUPAC], thresholds: &Thresholds) -> Self {

        // Split the mining between all available cores
        let chunk_size = genome.n / (num_cpus::get() * MINER_CPU_MULT);
     
        // Visualize chunks progress bars
        let bar = MultiProgress::new();
        let sty = ProgressStyle::with_template("{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}").unwrap()
            .progress_chars("##-");

        // Mine all aligments below the gap and mismatch threshold
        let indices: Vec<usize> = (0..genome.n).collect();
        let alignments:Vec<Alignment<T>> = indices 
            .par_chunks(chunk_size)
            .into_par_iter()
            .flat_map_iter(|items| {

                // Visualize this chunk progress
                let b = bar.add(ProgressBar::new(items.len() as u64));
                b.set_style(sty.clone());
                b.set_message("mining");

                let mut local_alignments: Vec<Alignment<T>> = vec![];
                let mut local_metrics = Metrics::empty();

                let mut local_miner: Miner<T> = Miner::new(thresholds, guide.len(), genome.seq_len);
                for &item_id in items {
                    
                    let seq_id = genome.ids[item_id];
                    let seq_beg = item_id * genome.seq_len;
                    let seq = &genome.sequences[seq_beg .. seq_beg + genome.seq_len];

                    local_metrics.aggregate(local_miner.mine(guide, seq, seq_id));
                    local_alignments.extend_from_slice(local_miner.solutions());
                    b.inc(1); 
                }

                // Visualize final statistics
                b.finish_with_message(format!("Complete! [{:<6}] (pruned: {})", local_alignments.len(), local_metrics.pruned));

                local_alignments

            }).collect();

        Self { alignments }
    }
}

use std::io::Write;
use std::fs::OpenOptions;
use std::sync::atomic::{AtomicU64, Ordering};
use std::os::unix::fs::FileExt;
use std::sync::Arc;

unsafe fn bytes_of<T>(val: &T) -> &[u8] {
    let size = core::mem::size_of::<T>();
    if size == 0 {
        return &[];
    }
    unsafe { core::slice::from_raw_parts(val as *const T as *const u8, size) }
}

pub fn parallel_mine_to_file<T>(genome: &FilteredGenome, guide: &[IUPAC], thresholds: &Thresholds, output: &PathBuf)
where
    T: PrimInt + Unsigned + Zero + Send
{

    // Create output file and atomic write offset
    let global_write_offset = Arc::new(AtomicU64::new(0));
    let file = Arc::new(OpenOptions::new()
        .create(true)
        .truncate(true)
        .write(true)
        .open(output)
        .expect("unable to create output file"));

    // CSV header
    //writeln!(file, "sequence_id,strand,offset,cigarx")
    //    .expect("unable to write output file"); 

    // Split the mining between all available cores
    let chunk_size = genome.n / (num_cpus::get() * MINER_CPU_MULT);
     
    // Visualize chunks progress bars
    let bar = MultiProgress::new();
    let sty = ProgressStyle::with_template("{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}").unwrap()
        .progress_chars("##-");

    // Mine all aligments below the gap and mismatch threshold
    let indices: Vec<usize> = (0..genome.n).collect();
    indices.par_chunks(chunk_size)
        .into_par_iter()
        .for_each(|items| {

            // Visualize this chunk progress
            let b = bar.add(ProgressBar::new(items.len() as u64));
            b.set_style(sty.clone());
            b.set_message("mining");

            let mut local_alignments: Vec<Alignment<T>> = vec![];
            let mut local_metrics = Metrics::empty();
            let mut skipped: usize = 0;

            let mut local_miner: Miner<T> = Miner::new(thresholds, guide.len(), genome.seq_len);
            for &item_id in items {
                    
                let seq_id = genome.ids[item_id];
                let seq_beg = item_id * genome.seq_len;
                let seq = &genome.sequences[seq_beg .. seq_beg + genome.seq_len];

                // Skip sequences with more than 3 Ns
                if seq.iter().filter(|&b| b.0 == 0b1111).count() > 4 {
                    b.set_message("Skipped sequence with more that 5 Ns!");
                    skipped += 1;
                    b.inc(1);
                    continue;
                }

                local_metrics.aggregate(local_miner.mine(guide, seq, seq_id));
                local_alignments.extend_from_slice(local_miner.solutions());
                b.inc(1); 
            }

            b.set_message("saving to file");

            let mut result = std::io::BufWriter::new(Vec::new());
            for alignment in &local_alignments {
                writeln!(result, "{},+,{},{}",
                    alignment.seq_id,
                    alignment.offset,
                    cigarx_decode(alignment.cigarx),
                ).expect("unable to write output file");
            }


            // Reserve file space atomically
            let bytes = result.into_inner().unwrap();
            let offset = global_write_offset.fetch_add(bytes.len() as u64, Ordering::SeqCst);

            // Write output bytes
            file.write_at(&bytes, offset)
                .expect("write failed");

            // Visualize final statistics
            b.finish_with_message(format!("Complete! [{:<6}] (pruned: {}, skipped: {})", 
                local_alignments.len(), local_metrics.pruned, skipped));
        });


}
*/
