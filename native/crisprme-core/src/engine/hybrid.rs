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
use tracing::{info, trace};

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
    sequence_producer: SequenceSend,
    sequence_consumer: SequenceRecv,
    alignment_producer: AlignmentSend,
    alignment_consumer: AlignmentRecv
}

impl HybridEngine {

    /// Create engine (CUDA init disabled for debug).
    pub fn new(alignment: AlignmentParams) -> Self {
        info!("running alignment (DEBUG): {alignment:#?}");

        // Window ring
        let sequence_batch_bytes = alignment.sequence_batch_size * (alignment.sequence_len + 4);
        let (sequence_producer, sequence_consumer) =
            ring_buffer::<SequenceRingBatch>(12, sequence_batch_bytes, true);

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

        // Drop internal producers. Miner never uses alignment_producer in debug mode.
        // drop(alignment_producer);

        Self {
            sequence_producer, 
            sequence_consumer,
            alignment_producer, 
            alignment_consumer
        }
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

    // TargetBatcher -> send -> SequenceRingBatch
    pub fn send(&mut self, batcher: &TargetBatcher) {
        let mut batch = self.sequence_producer.acquire();
        // Fill
        output.commit(batch);
    }

    pub fn receive(&mut self, batcher: &mut TargetBatcher) {
        if let Some(batch) = input.recv() {
            // TODO: Add alignments to TargetBatcher
            input.finish(batch);
        }
    }
}

/*
/// Thread responsible for DEBUG mining: just prints sequences that would be mined.
#[tracing::instrument(skip_all)]
fn th_miner_debug_print_windows(alignment: AlignmentParams, input: WindowRecv) {
    info!("DEBUG miner started (printing unique windows)");

    let mut batch_idx = 0usize;

    while let Some(wbatch) = input.recv() {
        batch_idx += 1;

        info!(
            "received window batch #{batch_idx} (windows={}, occs={})",
            wbatch.windows_len(),
            wbatch.occ_len()
        );

        let slen = wbatch.descriptor.sequence_len;
        let windows = wbatch.windows_iupac();

        // Optional: cap printing to avoid flooding logs for huge batches
        let max_print = alignment.debug_max_windows_to_print.unwrap_or(50);
        let n = wbatch.windows_len().min(max_print);

        for w in 0..n {
            let start = w * slen;
            let end = start + slen;
            let seq_slice = &windows[start..end];

            // Convert to string for debug
            let seq_string: String = seq_slice.iter().map(|b| b.to_utf8()).collect();

            println!("[BATCH {batch_idx}] window_id={w} seq={seq_string}");
        }

        if wbatch.windows_len() > n {
            println!(
                "[BATCH {batch_idx}] ... (printed {n}/{})",
                wbatch.windows_len()
            );
        }

        input.finish(wbatch);
    }

    info!("DEBUG miner closed");
}
*/