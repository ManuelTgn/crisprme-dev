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

// use crate::bindings;
use crate::alignment::alignment::Alignment;
use crate::memory::batch::{AlignmentRingBatch, WindowRingBatch};
use crate::memory::ring::{ring_buffer, Consumer, Producer};
use crate::storage::writer::AlignmentBatchWriter;

use super::params::AlignmentParams;

use std::time::Instant;
use tracing::info;

type WindowSend = Producer<WindowRingBatch>;
type WindowRecv = Consumer<WindowRingBatch>;

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

/// Handle returned to the caller so they can feed WindowRingBatch slots and then join threads.
pub struct EngineHandle {
    pub window_producer: WindowSend,
    threads: Vec<std::thread::JoinHandle<()>>,
}

impl EngineHandle {
    /// Drop the producer (or let it go out of scope) before calling join().
    pub fn join(self) {
        // When all producers are dropped, miner and writers will drain and exit.
        self.threads.into_iter().for_each(|h| h.join().unwrap());
    }
}

pub struct HybridEngine {
    device_id: i32,
}

impl HybridEngine {
    /// Create engine (CUDA init disabled for debug).
    pub fn new(device_id: i32) -> Self {
        // bindings::miner::initialize(device_id);
        Self { device_id }
    }

    /// Start an alignment that consumes WindowRingBatch.
    ///
    /// Returns an EngineHandle containing a Window producer; the caller fills/commits window batches.
    pub fn execute(&self, alignment: AlignmentParams) -> EngineHandle {
        info!("running alignment (DEBUG): {alignment:#?}");

        // Window ring
        let window_batch_bytes = alignment.window_batch_bytes;
        let (window_producer, window_consumer) =
            ring_buffer::<WindowRingBatch>(12, window_batch_bytes, true);

        // Alignment ring (unused in debug but kept for drop-in compatibility)
        let alignment_batch_bytes = alignment.alignment_batch_size * size_of::<Alignment>();
        let (alignment_producer, alignment_consumer) =
            ring_buffer::<AlignmentRingBatch>(6, alignment_batch_bytes, true);

        // Threads
        let mut threads = Vec::new();

        // DEBUG miner thread: consumes window batches and prints sequences
        engine_spawn_thread!(
            &mut threads,
            th_miner_debug_print_windows,
            alignment = alignment.clone(),
            input = window_consumer.clone()
        );

        // Writer threads (they will just idle/exit because no output is committed)
        let writer = AlignmentBatchWriter::open(&alignment.output_file);
        for _ in 0..6 {
            engine_spawn_thread!(
                &mut threads,
                th_writer,
                alignment = alignment.clone(),
                input = alignment_consumer.clone(),
                writer = writer.clone()
            );
        }

        // Drop internal producers. Miner never uses alignment_producer in debug mode.
        drop(alignment_producer);

        EngineHandle {
            window_producer,
            threads,
        }
    }
}

impl Drop for HybridEngine {
    fn drop(&mut self) {
        // In debug mode CUDA init is disabled, so shutdown must be disabled too.
        // bindings::miner::shutdown(self.device_id);
        let _ = self.device_id;
    }
}

/// Thread responsible for writing alignments batches
#[tracing::instrument(skip_all)]
fn th_writer(_alignment: AlignmentParams, input: AlignmentRecv, writer: AlignmentBatchWriter) {
    info!("writer started, {writer:#?}");
    let mut total_writes = 0usize;

    while let Some(batch) = input.recv() {
        total_writes += batch.len();

        let now = Instant::now();
        writer.write_from_memory(batch.alignments());
        info!(
            "wrote {} alignments to file [{} ms]",
            batch.len(),
            now.elapsed().as_millis()
        );

        input.finish(batch);
    }

    info!("writer closed, total writes: {total_writes}");
}

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
