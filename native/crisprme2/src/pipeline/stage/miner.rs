use bytemuck::Pod;
use columnar::{MemoryPool, pipeline::{Emit, Stage, PipelineError}, Share};
use crossbeam_channel::Receiver;
use itertools::izip;
use rand::Rng;

use crate::{bindings::{self, cuda}, model::{alignment::{SeqMinedBatch, SeqMinedFrame}, cigarx::{Cigarx, Cigarx64, CigarxOp}, input::{SeqBatch, SeqRowIdx}}, sequence::iupac::Iupac};

// ---------------------------------------------------------------------------
// Fake miner
// ---------------------------------------------------------------------------

/// Fake miner for now, it only checks for match/mismatch
pub struct Miner { pool: MemoryPool }

impl Miner {
    pub fn new(pool: &MemoryPool) -> Self {
        Self { pool: pool.clone() }
    }
}

impl Stage for Miner {

    type I = SeqBatch;
    type O = SeqMinedBatch;

    fn name() -> &'static str { "Miner" }

    #[tracing::instrument(name = "pipeline:miner", skip_all)]
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), PipelineError> {
        
        let mut sequences  = input.sequences.share();
        let mut occurences = input.occurences.share();

        input.sequences.with_cols(|cols| {

            let rows = cols.content.rows();
            tracing::info!("received {} rows to mine", rows);

            let mut mined = SeqMinedFrame::alloc(&self.pool, rows);
            mined.with_cols(|mut mined| {
                let zipped = izip!(
                    cols.content.iter(),
                    mined.seq_row_idx.iter_mut(),
                    mined.cigarx.iter_mut(),
                    mined.offset.iter_mut()
                );

                for (sequence, mined_seq_row, cigarx, offset) in zipped {
                     
                    *mined_seq_row = rand::rng().random_range(0..rows) as u32;

                    *cigarx = Cigarx64::default();
                    *offset = 2;

                    for j in 0..input.guide.len() {
                        if j % 3 == 0 {
                            cigarx.push(CigarxOp::Deletion);
                        } else if j % 5 == 0 {
                            cigarx.push(CigarxOp::Insertion);
                        } else if sequence[j].matches(input.guide[j]) {
                            cigarx.push(CigarxOp::Match);
                        } else {
                            cigarx.push(CigarxOp::Mismatch);
                        }
                    }
                }

            });

            emitter.emit(SeqMinedBatch {
                guide: input.guide.clone(),
                sequences: sequences.share(),
                occurences: occurences.share(),
                mined,
            })
        })
    }
}

// ---------------------------------------------------------------------------
// GPU miner
// ---------------------------------------------------------------------------

struct GpuBuffer<T> {
    /// Maximum number of `T` that can fit
    capacity: usize,
    /// Pointer to GPU memory
    dptr: *mut T
}

impl<T> GpuBuffer<T> {
    pub fn alloc(capacity: usize) -> Self {
        let dptr = bindings::cuda::malloc::<T>(capacity);
        Self { capacity, dptr }
    }
}

impl<T> Drop for GpuBuffer<T> {
    fn drop(&mut self) {
        bindings::cuda::free(self.dptr);
    }
}

// This should be fine
unsafe impl<T: Pod> Send for GpuBuffer<T> { }

/// Contains all GPU buffers used by a miner
struct GpuMinerBuffers {
    // Staging area of GPU memory for input sequences
    src_buffer_seq: GpuBuffer<Iupac>,
    // Staging area of GPU memory for output mined alignments
    dst_buffer_cigarx: GpuBuffer<Cigarx64>,
    dst_buffer_idx:    GpuBuffer<SeqRowIdx>,
    dst_buffer_offset: GpuBuffer<u8>,
}

pub struct GpuMiner {

    pool: MemoryPool,
    gpu: u32,

    /// Buffers on the GPU
    buffers: Option<GpuMinerBuffers>,

    /// Number of sequences that can fit in the src buffer
    src_seq_capacity: usize,
    /// Maximum size of each sequence
    src_seq_stride: usize,

    // Maximum number of alignments for each iteration of the miner
    dst_capacity: usize,

    // Total alignments mined
    total_mined: usize,
}

impl GpuMiner {

    pub fn new(pool: &MemoryPool, 
        src_seq_capacity: usize, src_seq_stride: usize, 
        dst_capacity: usize, gpu: u32) -> Self 
    {
        Self { 
            pool: pool.clone(),
            total_mined: 0, 
            buffers: None,
            src_seq_capacity,
            src_seq_stride,
            dst_capacity,
            gpu 
        }
    }
}

impl Stage for GpuMiner {

    type I = SeqBatch;
    type O = SeqMinedBatch;

    fn name() -> &'static str { "GpuMiner" }

    // Called by the asigned thread
    #[tracing::instrument(name = "pipeline:gpu_miner", skip_all)]
    fn initialize(&mut self) {
        
        // Initialize mining kernel for this gpu/thread
        bindings::miner::initialize(self.gpu);

        // Allocate GPU memory for the asigned GPU
        self.buffers = Some(GpuMinerBuffers { 
            src_buffer_seq:    GpuBuffer::alloc(self.src_seq_capacity * self.src_seq_stride), 
            dst_buffer_cigarx: GpuBuffer::alloc(self.dst_capacity), 
            dst_buffer_idx:    GpuBuffer::alloc(self.dst_capacity), 
            dst_buffer_offset: GpuBuffer::alloc(self.dst_capacity) 
        });
    }

    #[tracing::instrument(name = "pipeline:gpu_miner", skip_all)]
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), PipelineError> {

        // Allocated GPU buffers
        let buffers = self.buffers.as_mut().expect("Miner buffers not allocated!");

        // Prepare miner with context
        bindings::miner::prepare(&input.guide, input.seq_len, &input.thresholds);
        
        // Copy data from column region to staging buffer
        let mut src_seq_count = 0;
        input.sequences.with_cols(|cols| {
            src_seq_count = cols.content.rows();

            tracing::debug!("uploading columns chunks to GPU");

            // How many bytes for each sequence element
            let stride = cols.content.row_bytes();
            assert_eq!(stride, self.src_seq_stride);
            
            let mut offset = 0;
            for region in cols.content.chunks() {
                assert!(offset < self.src_seq_capacity, "offset is {offset}, but max capacity is {}", self.src_seq_capacity);
                //tracing::info!("moving {} bytes to src_seq_buffer", region.len() * stride);
                bindings::cuda::memcpy_to_gpu::<Iupac>(
                    region.as_ptr() as _,
                    unsafe { buffers.src_buffer_seq.dptr.add(offset) },
                    region.len() * stride
                );
                offset += region.len() * stride;
            }
        });

        // Mine alignments in a loop — kernel yields partial results when output buffer is full
        let mut sequences  = input.sequences.share();
        let mut occurences = input.occurences.share();

        loop {
            tracing::debug!("running columnar mine kernel");
            let (finish, mined_count) = bindings::miner::mine(
                buffers.src_buffer_seq.dptr as _,
                src_seq_count as u32,
                buffers.dst_buffer_cigarx.dptr as _,
                buffers.dst_buffer_idx.dptr as _,
                buffers.dst_buffer_offset.dptr as _,
                self.dst_capacity as u32
            );

            self.total_mined += mined_count;
            tracing::debug!("mined {} alignments (finish={})", mined_count, finish);
            if mined_count != 0 {

                // Download results from GPU and emit
                let mut mined = SeqMinedFrame::alloc(&self.pool, mined_count);
                mined.with_cols(|mut cols| {
                    let mut offset = 0;
                    for region in cols.cigarx.chunks_mut() {
                        bindings::cuda::memcpy_to_cpu::<Cigarx64>(
                            region.as_ptr() as _,
                            unsafe { buffers.dst_buffer_cigarx.dptr.add(offset) },
                            region.len()
                        );
                        offset += region.len();
                    }

                    offset = 0;
                    for region in cols.seq_row_idx.chunks_mut() {
                        bindings::cuda::memcpy_to_cpu::<SeqRowIdx>(
                            region.as_ptr() as _,
                            unsafe { buffers.dst_buffer_idx.dptr.add(offset) },
                            region.len()
                        );
                        offset += region.len();
                    }

                    offset = 0;
                    for region in cols.offset.chunks_mut() {
                        bindings::cuda::memcpy_to_cpu::<u8>(
                            region.as_ptr() as _,
                            unsafe { buffers.dst_buffer_offset.dptr.add(offset) },
                            region.len()
                        );
                        offset += region.len();
                    }
                });

                emitter.emit(SeqMinedBatch {
                    guide: input.guide.clone(),
                    sequences: sequences.share(),
                    occurences: occurences.share(),
                    mined
                })?;
            }
            if finish { break; }
        }

        bindings::miner::post_mine();
        Ok(())
    }

    fn shutdown(&mut self) -> Result<(), PipelineError> {
        bindings::miner::shutdown(self.gpu);
        Ok(())
    }
}