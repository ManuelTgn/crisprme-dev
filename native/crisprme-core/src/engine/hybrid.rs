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

        Self {
            params,
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

#[pymethods]
impl HybridEngine {

    /// Send a batch of query sequences to the alignment engine.
    ///
    /// Acquires the next free [`SequenceRingBatch`] slot (blocking if the ring
    /// is full), copies all window keys from `batcher` into it, assigns
    /// sequential IDs, then commits the slot for the miner thread to consume.
    ///
    /// # Errors
    ///
    /// Returns [`PyBufferError`] if the engine has already been shut down
    /// (i.e. [`Drop`] has been called and `sequence_producer` is `None`).
    ///
    /// # Panics
    ///
    /// Panics if any individual window key would write beyond the allocated
    /// IUPAC buffer (`offset + window_len > input.len()`).
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
            sequence_len: self.params.sequence_len,
            global_offset: 0,
            batcher_id: Some(batcher.id()),
        });

        info!("sent batcher window keys (size = {})", batcher.get_window_count());
        Ok(())
    }


    /// Retrieve the next alignment batch, blocking until one is available.
    ///
    /// Returns an [`AlignmentBatchView`] wrapping the ring buffer slot.
    /// The slot is **not** returned to the ring until the view is dropped, so
    /// callers should avoid holding onto it longer than necessary to prevent
    /// stalling the miner thread.
    ///
    /// Returns an empty view if the alignment ring has been closed (engine
    /// is shutting down and all results have been consumed).
    pub fn receive_blocking(&mut self) -> PyResult<AlignmentBatchView> {
        if let Some(batch) = self.alignment_consumer.recv() {
            return Ok(AlignmentBatchView::new(batch));
            // NOTE: This is called on python object drop
            // self.alignment_consumer.finish(batch);
        }
        Ok(AlignmentBatchView::empty())
    }


    /// Attempt to retrieve the next alignment batch without blocking.
    ///
    /// Returns an empty [`AlignmentBatchView`] immediately if no batch is
    /// currently available, rather than waiting.  Suitable for polling loops
    /// where the caller cannot afford to block the GIL.
    ///
    /// Errors from the ring (e.g. a poisoned state) are logged and converted
    /// to an empty view rather than propagated.
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