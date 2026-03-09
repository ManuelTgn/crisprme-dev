//! # Hybrid Alignment Engine
//!
//! [`HybridEngine`] orchestrates a multithreaded CUDA alignment pipeline.
//! It connects a Python-facing API to background threads via lock-free ring
//! buffers, keeping data movement off the GIL.
//!
//! ## Architecture
//!
//! ```text
//!  Python (send)          th_miner_cuda              Python (receive)
//!      │                       │                           │
//!      ▼                       ▼                           ▼
//! SequenceRingBatch ──────► [miner] ──► AlignmentRingBatch ──► AlignmentBatchView
//! (sequence_producer)    CUDA mine     (alignment_consumer)
//! ```
//!
//! ## Shutdown
//!
//! Dropping [`HybridEngine`] triggers an orderly shutdown:
//!
//! 1. `sequence_producer` is dropped, closing the input ring.
//! 2. The miner thread detects the closed ring (`recv()` returns `None`) and exits.
//! 3. The miner dropping `output` closes the alignment ring.
//! 4. [`Drop`] joins all threads before returning.

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

use std::collections::HashMap;
use std::thread::JoinHandle;
use crate::bindings;
use crate::alignment::alignment::Alignment;
use crate::memory::batch::{AlignmentRingBatch, SequenceRingBatch};
use crate::memory::ring::{ring_buffer, Consumer, Producer};
use crate::memory::arena::Arena;
use crate::storage::writer::{AlignmentBatchDescr, AlignmentBatchWriter};
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


/// Spawns a thread and pushes its [`JoinHandle`] onto `$threads`.
///
/// Each `$name = $expr` binding is moved into the thread closure, avoiding
/// the need for manual `let x = x.clone()` boilerplate at every call site.
///
/// # Example
/// ```rust,ignore
/// engine_spawn_thread!(
///     &mut threads,
///     th_miner_cuda,
///     params = params.clone(),
///     input  = sequence_consumer.clone(),
///     output = alignment_producer.clone(),
/// );
/// ```
macro_rules! engine_spawn_thread {
    ($threads:expr, $func:ident, $($name:ident = $expr:expr),* $(,)?) => {{
        $( let $name = $expr; )+
        $threads.push(std::thread::spawn(move || {
            $func($($name),*);
        }));
    }};
}

/// A CUDA-backed sequence alignment engine exposed to Python via PyO3.
///
/// Internally manages two ring buffers:
///
/// - **sequence ring** — carries [`SequenceRingBatch`] objects from Python
///   into the miner thread.
/// - **alignment ring** — carries [`AlignmentRingBatch`] results from the
///   miner thread back to Python.
///
/// Both rings are bounded and back-pressure: [`send`][HybridEngine::send]
/// will block if the sequence ring is full, and
/// [`receive_blocking`][HybridEngine::receive_blocking] will block until a
/// result batch is available.
///
/// # Python usage
/// ```python
/// engine = HybridEngine(params)
/// engine.send(batcher)
/// batch = engine.receive_blocking()
/// # batch is a zero-copy memoryview over the alignment ring slot
/// del engine  # triggers orderly shutdown
/// ```
#[pyclass]
pub struct HybridEngine {
    /// Alignment parameters used to configure the CUDA miner and to fill
    /// [`SequenceBatchDescr`] on every [`send`][HybridEngine::send] call.
    params: AlignmentParams,
    /// Write end of the sequence ring. Wrapped in `Option` so [`Drop`] can
    /// take and drop it to signal shutdown before joining threads.
    sequence_producer: Option<SequenceSend>,
    /// Read end of the alignment ring. Consumed by [`receive`][HybridEngine::receive] /
    /// [`receive_blocking`][HybridEngine::receive_blocking].
    alignment_consumer: AlignmentRecv,
    /// Handles for all background threads. Joined in [`Drop`] to guarantee
    /// threads have exited before the engine is considered destroyed.
    threads: Vec<JoinHandle<()>>
}

#[pymethods]
impl HybridEngine {

    /// Construct and start the engine.
    ///
    /// Allocates both ring buffers and spawns the CUDA miner thread.
    /// The engine is ready to accept [`send`][HybridEngine::send] calls
    /// immediately after this returns.
    ///
    /// # Ring buffer sizing
    /// - Sequence ring: 6 slots × `sequence_batch_size × (sequence_len + 4)` bytes.
    /// - Alignment ring: 6 slots × `alignment_batch_size × size_of::<Alignment>()` bytes.
    #[new]
    pub fn new(params: AlignmentParams) -> Self {
        
        info!("running alignment (DEBUG): {params:#?}");

        // Window ring
        let sequence_batch_bytes = params.sequence_batch_size * (params.sequence_len + 4);
        let (sequence_producer, sequence_consumer) =
            ring_buffer::<SequenceRingBatch>(6, sequence_batch_bytes, true);

        // Alignment ring (unused in debug but kept for drop-in compatibility)
        let alignment_batch_bytes = params.alignment_batch_size * size_of::<Alignment>();
        let (alignment_producer, alignment_consumer) =
            ring_buffer::<AlignmentRingBatch>(6, alignment_batch_bytes, true);

        // Threads
        let mut threads = Vec::new();

        // Spawn CUDA device driver
        // DEBUG miner thread: consumes window batches and prints sequences
        engine_spawn_thread!(
            &mut threads,
            th_miner_cuda,
            alignment = params.clone(),
            input = sequence_consumer.clone(),
            output = alignment_producer.clone()
        );

        // Spawn aggregator
        engine_spawn_thread!(
            &mut threads,
            th_aggregator,
            alignment = params.clone(),
            rx = alignment_consumer.clone()
        );

        Self {
            params,
            sequence_producer: Some(sequence_producer),
            alignment_consumer,
            threads
        }
    }

    /// Send a batch of sequences to the miner thread for alignment.
    pub fn send(&mut self, batcher: &mut TargetBatcher) -> PyResult<()> {

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

        // Create alignment stream
        let (reply_tx, reply_rx) = crossbeam::channel::unbounded();
        batcher.set_alignment_stream(reply_rx);

        // Send the filled buffer
        producer.commit_with_descriptor(batch, SequenceBatchDescr {
            sequence_count: batcher.get_window_count(),
            sequence_len: self.params.sequence_len,
            output_tx: Some(reply_tx),
            batcher_id: batcher.id(),
            global_offset: 0,
        });

        info!("sent batcher window keys (size = {})", batcher.get_window_count());
        Ok(())
    }


    /// Retrieve the next alignment batch, blocking until one is available.
    pub fn receive_blocking(&mut self, batcher: &mut TargetBatcher) -> PyResult<Vec<AlignmentBatchView>> {
        if let Some(rx) = batcher.extract_alignment_rx() {
            trace!("waiting for alignment batches from batcher {}", batcher.id());
            return rx.iter()
                .map(|batch| {
                    assert_eq!(batch.descriptor.batcher_id, batcher.id(), "received batch for wrong batcher");
                    Ok(AlignmentBatchView::new(batch))
                })
                .collect();
        }
        warn!("no alignment stream available for batcher {}, returning empty result", batcher.id());
        Ok(vec![])
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

/// Background thread: drives the CUDA miner for one GPU device.
///
/// Continuously reads [`SequenceRingBatch`] items from `input`, runs the CUDA
/// miner over them for both strands, and writes [`AlignmentRingBatch`] results
/// to `output`.  Exits when `input.recv()` returns `None`, which happens once
/// [`HybridEngine`] drops its `sequence_producer`.
///
/// # Per-batch flow
///
/// For each incoming sequence batch the thread performs two full mining passes:
///
/// 1. **Positive strand** (`+`) — [`bindings::miner::pre_mine`] configures the
///    guide, then [`bindings::miner::mine`] is called in a loop until it signals
///    completion, emitting one [`AlignmentRingBatch`] per loop iteration.
/// 2. **Negative strand** (`-`) — same flow using the reverse-complement guide.
///
/// After both passes, [`bindings::miner::post_mine`] cleans up GPU state and
/// `input.finish` returns the sequence slot to the ring for reuse.
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
                output.commit_with_descriptor(alignments, AlignmentBatchDescr {
                    batcher_id: batch.descriptor.batcher_id,
                    alignments_count,
                    // Only signal completion on the last batch of the negative strand to avoid triggering early finalization in the consumer
                    output_tx: None,
                });

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
                output.commit_with_descriptor(alignments, AlignmentBatchDescr {
                    batcher_id: batch.descriptor.batcher_id,
                    // NOTE: This can signal to the consumer that this is the last batch for this window batch, 
                    // so it can trigger any necessary finalization steps
                    output_tx: match complete {
                        true  => batch.descriptor.output_tx.take(),
                        false => None,
                    },
                    alignments_count,
                });

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

/// Aggregates alignment batches and forward them to the source TargetBatcher.
#[tracing::instrument(skip_all)]
fn th_aggregator(_pipeline: AlignmentParams, rx: AlignmentRecv) {

    let mut map: HashMap<usize, Vec<AlignmentRingBatch>> = HashMap::new();
    while let Some(mut batch) = rx.recv() {

        let tx = batch.descriptor.output_tx.take();
        let batcher_id = batch.descriptor.batcher_id;

        trace!("received batch to aggregate (target_batcher: {})", batcher_id);

        // Add batch to the corresponding batcher aggregate 
        map.entry(batcher_id)
            .or_default()
            .push(batch);

        // If there is a transmitter, this means that this is the last batch of the window batch, 
        // so we can trigger the finalization of the target batcher and submit all results
        if let Some(tx) = tx {
            let aggregates = map.remove(&batcher_id)
                .expect("alignment batch without aggregate");

            for inner in aggregates {
                if let Err(e) = tx.send(inner) {
                    error!("failed to send alignment batch to batcher: {}", e);
                }
            }
        }
    }
}