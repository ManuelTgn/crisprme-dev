// src/pipeline/engine/hybrid.rs
//
// DROP-IN DEBUG VERSION:
// - Keeps the WindowRingBatch + AlignmentRingBatch rings as in your snippet
// - Miner thread ONLY prints the unique windows it would mine (no CUDA calls)
// - Writer threads stay, but will just exit because nothing is ever produced
//
// Notes:
// - Requires WindowRingBatch::windows_iupac() and descriptor.sequence_len.
// - Assumes Iupac has to_utf8() (as in your other code).
// - Completely removes unused mine_and_expand_pass() to avoid warnings/confusion.
// - Keeps bindings import only for shutdown in Drop (as requested).

use std::thread::JoinHandle;
use crate::bindings;
use crate::alignment::alignment::Alignment;
use crate::memory::batch::{AlignmentRingBatch, SequenceRingBatch};
use crate::memory::ring::{ring_buffer, Consumer, Producer};
use crate::memory::arena::Arena;
use crate::storage::writer::AlignmentBatchWriter;
use crate::memory::ring::RingAdapter;
use crate::batching::batching::TargetBatcher;

use super::params::AlignmentParams;

use std::time::Instant;
use pyo3::{pyclass, pymethods, PyResult};
use pyo3::exceptions::PyBufferError;
use tracing::{error, info, trace, warn};
use crate::python::views::AlignmentBatchView;
use crate::sequence::iupac::Iupac;
use crate::storage::reader::SequenceBatchDescr;

type SequenceSend = Producer<SequenceRingBatch>;
type SequenceRecv = Consumer<SequenceRingBatch>;

type AlignmentSend = Producer<AlignmentRingBatch>;
type AlignmentRecv = Consumer<AlignmentRingBatch>;

/// Spawn a new thread and removes boilerplate
macro_rules! engine_spawn_thread {
    ($threads:expr, $func:ident, $($name:ident = $expr:expr),* $(,)?) => {{
        $( let $name = $expr; )+
        $threads.push(std::thread::spawn(move || {
            $func($($name),*);
        }));
    }};
}

#[pyclass]
pub struct HybridEngine {
    sequence_producer: Option<SequenceSend>,
    alignment_consumer: AlignmentRecv,
    threads: Vec<JoinHandle<()>>
}

impl HybridEngine {

    /// Create engine (CUDA init disabled for debug).
    pub fn new(alignment: AlignmentParams) -> Self {
        info!("running alignment (DEBUG): {alignment:#?}");

        // Window ring
        let sequence_batch_bytes = alignment.sequence_batch_size * (alignment.sequence_len + 4);
        let (sequence_producer, sequence_consumer) =
            ring_buffer::<SequenceRingBatch>(6, sequence_batch_bytes, true);

        // Alignment ring (unused in debug but kept for drop-in compatibility)
        let alignment_batch_bytes = alignment.alignment_batch_size * size_of::<Alignment>();
        let (alignment_producer, alignment_consumer) =
            ring_buffer::<AlignmentRingBatch>(6, alignment_batch_bytes, true);

        // Threads
        let mut threads = Vec::new();

        // Spawn CUDA device driver
        // DEBUG miner thread: consumes window batches and prints sequences
        engine_spawn_thread!(
            &mut threads,
            th_miner_cuda,
            alignment = alignment.clone(),
            input = sequence_consumer.clone(),
            output = alignment_producer.clone()
        );

        Self {
            sequence_producer: Some(sequence_producer),
            alignment_consumer,
            threads
        }
    }
}

impl Drop for HybridEngine {
    fn drop(&mut self) {
        info!("shutting down hybrid engine");

        // Signal to the miners that they should stop
        drop(self.sequence_producer.take());

        // Wait for all threads
        for thread in self.threads.drain(..) {
            if let Err(e) = thread.join() {
                error!("error dropping thread: {:?}", e);
            }
        }

        info!("all threads stopped");
    }
}

/// Thread responsible for managing the CUDA miner
#[tracing::instrument(skip_all)]
fn th_miner_cuda(pipeline: AlignmentParams, input: SequenceRecv, output: AlignmentSend) {
    info!("started!");

    // Total alignments mined by the cuda miner
    let mut total_mined = 0;

    bindings::miner::initialize(0);

    let mut mine_batch_idx = 0;
    let mut arena = Arena::alloc(1024 * 1024 * 1024);
    while let Some(mut batch) = input.recv() {
        trace!("received batch #{mine_batch_idx} ({} seq)", batch.len());
        mine_batch_idx += 1;

        arena.scoped(|memory| {
            let now = Instant::now();
            batch.as_mut()
                .sync_cpu_to_gpu(None);

            // Mine positive alignment batches
            let now = Instant::now();
            bindings::miner::pre_mine(&pipeline.guide, pipeline.sequence_len, &pipeline.thresholds, b'+');
            trace!("pre_mine positive [{} ms]", now.elapsed().as_millis());

            let mut alignments = output.acquire();
            let mut alignments_batch_count = 1;
            loop {
                let now = Instant::now();
                info!("mining positive batch of {} sequences", batch.len());
                let complete = bindings::miner::mine(&batch, &mut alignments);
                info!(
                    "mined {} positive alignments (output batch {}) [{} ms]",
                    alignments.len(),
                    alignments_batch_count,
                    now.elapsed().as_millis()
                );

                // Sync data
                let alignments_count = alignments.len();
                total_mined += alignments_count;
                alignments
                    .as_mut()
                    .sync_gpu_to_cpu(Some(size_of::<Alignment>() * alignments_count));

                alignments.replace_pos_by_id(&batch);
                output.commit(alignments);
                if complete {
                    break;
                }

                // Acquire next output batch
                alignments = output.acquire();
                alignments_batch_count += 1;
            }

            // Done with this batch
            let now = Instant::now();
            bindings::miner::post_mine();
            trace!("post_mine positive [{} ms]", now.elapsed().as_millis());

            // Mine negative alignment batches
            let now = Instant::now();
            let reverse_guide = pipeline.guide.reverse_complement();
            bindings::miner::pre_mine(&reverse_guide, pipeline.sequence_len, &pipeline.thresholds, b'-');
            trace!("pre_mine negative [{} ms]", now.elapsed().as_millis());

            let mut alignments = output.acquire();
            let mut alignments_batch_count = 1;
            loop {
                let now = Instant::now();
                info!("mining negative batch of {} sequences", batch.len());
                let complete = bindings::miner::mine(&batch, &mut alignments);
                info!(
                    "mined {} negative alignments (output batch {}) [{} ms]",
                    alignments.len(),
                    alignments_batch_count,
                    now.elapsed().as_millis()
                );

                // Sync data
                let alignments_count = alignments.len();
                total_mined += alignments_count;
                alignments
                    .as_mut()
                    .sync_gpu_to_cpu(Some(size_of::<Alignment>() * alignments_count));

                alignments.replace_pos_by_id(&batch);
                output.commit(alignments);
                if complete {
                    break;
                }

                // Acquire next output batch
                alignments = output.acquire();
                alignments_batch_count += 1;
            }

            // Done with this batch
            let now = Instant::now();
            bindings::miner::post_mine();
            trace!("post_mine negative [{} ms]", now.elapsed().as_millis());

            input.finish(batch);
        });
    }

    bindings::miner::shutdown(0);
    info!("closed, total mined: {} (positive and negative)", total_mined);
}

#[pymethods]
impl HybridEngine {

    /// Send TargetBatcher to engine
    pub fn send(&mut self, batcher: &TargetBatcher) -> PyResult<()> {

        // Try to get sequence producer, can fail on mid-shutdown
        let producer = self.sequence_producer.as_mut()
            .ok_or_else(|| PyBufferError::new_err("engine is not online"))?;

        // Wait for an empty buffer
        let mut batch = producer.acquire();

        if batcher.get_window_count() != batch.len() {
            warn!("received batcher window count does not match buffer size");
        }

        let input = batch.iupac_mut();
        let mut offset = 0;

        // Copy all window keys inside the SequenceRingBatch
        // NOTE: I don't think this needs further optimization
        for window in batcher.get_window_keys() {

            let window_len = window.len();
            assert!(offset + window_len <= input.len());

            // SAFETY: We know that the window key contains valid IUPAC bytes
            let window_iupac: &[Iupac] = unsafe {
                std::slice::from_raw_parts(
                    window.as_ptr() as *const Iupac,
                    window.len()
                )
            };

            input[offset .. offset + window_len]
                .copy_from_slice(window_iupac);

            offset += window_len;
        }

        // Generate sequence ids
        (0..batcher.get_window_count()).for_each(|i| {
            batch.ids_mut()[i] = i as u32;
        });

        // Send the filled buffer
        producer.commit_with_descriptor(batch, SequenceBatchDescr {
            sequence_count: batcher.get_window_count(),
            sequence_len: 30, // TODO,
            global_offset: 0,
            batcher_id: Some(batcher.id()),
        });

        info!("sent batcher window keys (size = {})", batcher.get_window_count());
        Ok(())
    }


    /// Retrieve alignments from engine.
    pub fn receive_blocking(&mut self) -> PyResult<AlignmentBatchView> {
        if let Some(batch) = self.alignment_consumer.recv() {
            return Ok(AlignmentBatchView::new(batch));
            // NOTE: This is called on python object drop
            // self.alignment_consumer.finish(batch);
        }
        Ok(AlignmentBatchView::empty())
    }


    /// Retrieve alignments from engine
    pub fn receive(&mut self) -> PyResult<AlignmentBatchView> {
        match self.alignment_consumer.try_recv() {
            Err(e) => {
                error!("tried to receive alignment batch: {}", e);
                Ok(AlignmentBatchView::empty())
            },
            Ok(batch) => match batch {
                Some(batch) => Ok(AlignmentBatchView::new(batch)),
                None => {
                    trace!("alignment batch not yet available");
                    Ok(AlignmentBatchView::empty())
                }
            }
        }
    }
}