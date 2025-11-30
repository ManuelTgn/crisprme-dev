use crate::common::sequence::Sequence;
use crate::{
    bindings,
    common::alignment::{Alignment, visualize},
    memory::{
        arena::Arena,
        batch::{AlignmentRingBatch, SequenceRingBatch},
        ring::{ring_buffer, Consumer, Producer, RingAdapter},
    },
    pipeline::PipelineDescriptor,
    storage::{reader::BinarySequenceBatchReader, writer::AlignmentBatchWriter},
};
use std::ops::DerefMut;
use std::time::Instant;
use tracing::{info, trace};

type BatchSend = Producer<SequenceRingBatch>;
type BatchRecv = Consumer<SequenceRingBatch>;

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

pub struct HybridEngine {}
impl HybridEngine {
    pub fn execute(&self, pipeline: PipelineDescriptor) {
        info!("running pipeline: {pipeline:#?}");

        // Ring buffer for sequence batches
        let sequence_batch_bytes = pipeline.sequence_batch_size * (pipeline.sequence_len + 4);
        let (batch_producer, batch_consumer) =
            ring_buffer::<SequenceRingBatch>(12, sequence_batch_bytes, true);

        // Ring buffer for alignment batches
        let alignment_batch_bytes = pipeline.alignment_batch_size * size_of::<Alignment>();
        let (alignment_producer, alignment_consumer) =
            ring_buffer::<AlignmentRingBatch>(6, alignment_batch_bytes, true);

        // Handles of the threads
        let mut threads = Vec::new();

        // Spawn I/O thread for reading input sequences
        engine_spawn_thread!(
            &mut threads,
            th_reader,
            pipeline = pipeline.clone(),
            producer = batch_producer.clone()
        );

        // Spawn thread for CUDA mining
        engine_spawn_thread!(
            &mut threads,
            th_miner_cuda,
            pipeline = pipeline.clone(),
            consumer = batch_consumer.clone(),
            producer = alignment_producer.clone()
        );

        /*
        // The remaining cores are used to mine on CPU
        let pipe = pipeline.clone();
        let consumer = batch_consumer.clone();
        let producer = alignment_producer.clone();
        threads.push(std::thread::spawn(move || {
            hybrid_cpu_miner_thread(pipe, consumer, producer);
        }));
        */

        // Spawn I/O thread for writing output alignment batches
        let writer = AlignmentBatchWriter::open(&pipeline.output_file);
        for _ in 0..6 {
            engine_spawn_thread!(
                &mut threads,
                th_writer,
                pipeline = pipeline.clone(),
                consumer = alignment_consumer.clone(),
                writer = writer.clone()
            );
        }

        // Drop producers to allow termination
        drop(alignment_producer);
        drop(batch_producer);

        // Wait for all threads to finish
        threads.into_iter().for_each(|h| h.join().unwrap());
    }
}

/// Thread responsible for reading the sequence batches
#[tracing::instrument(skip_all)]
fn th_reader(pipeline: PipelineDescriptor, output: BatchSend) {
    let reader = BinarySequenceBatchReader::open(
        &pipeline.sequence_file,
        pipeline.sequence_len,
        pipeline.sequence_batch_size,
    );

    info!("started, {reader:#?}");
    for (idx, real_size) in reader.batches() {

        // End at first batch
        //if idx > 0 { break; }

        let mut batch = output.acquire();

        let now = Instant::now();
        batch.descriptor = reader.describe(idx, real_size);
        reader.read_batch(idx, real_size, &mut batch);

        info!(
            "loaded sequence batch #{idx} ({real_size} seq) [{} ms]",
            now.elapsed().as_millis()
        );

        output.commit(batch);
    }

    info!("closed");
}

/// Thread responsible for writing alignments batches
#[tracing::instrument(skip_all)]
fn th_writer(pipeline: PipelineDescriptor, input: AlignmentRecv, writer: AlignmentBatchWriter) {
    info!("started, {writer:#?}");
    let mut total_writes = 0;
    let mut write_batch_idx = 0;
    while let Some(batch) = input.recv() {
        total_writes += batch.len();
        write_batch_idx += 1;

        let now = Instant::now();
        writer.write_from_memory(batch.alignments());
        info!(
            "wrote {} alignments to file [{} ms]",
            batch.len(),
            now.elapsed().as_millis()
        );

        input.finish(batch);
    }
    info!("closed, total writes: {total_writes}");
}

/// Thread responsible for managing the CUDA miner
#[tracing::instrument(skip_all)]
fn th_miner_cuda(pipeline: PipelineDescriptor, input: BatchRecv, output: AlignmentSend) {
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

            /*
            for seq in batch.sequences() {
                println!("seq: {seq}");
            }
            */

            /*
            let mut scores = memory.alloc_slice_fill(batch.len(), 255u8);
            batch.edit_distace_scores(&pipeline.guide, &mut scores);

            let mut mask = memory.alloc_slice_fill(batch.len(), false);
            for (i, seq) in batch.sequences().enumerate() {
                if seq.mutation_score() <= pipeline.mutation_max {
                    mask[i] = scores[i] <= pipeline.thresholds.ed() as u8;
                }
            }

            batch.apply_mask(mask.as_ref(), true);
            info!(
                "applyed mask for ed+mutation filter ({} sequences) [{} ms]",
                batch.len(),
                now.elapsed().as_millis()
            );
            */

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

                // SOLUTION
                /*
                let ids = batch.ids();
                let sequences: Vec<Sequence> = batch.sequences().collect();
                for alig in alignments.alignments() {

                    let g = pipeline.guide.to_string();
                    let s = sequences[alig.id as usize].to_string();
                    let i = ids[alig.id as usize]; 
                    let c = alig.cigarx.to_string();

                    println!("id: {}/{}, g: {}, s: {}, c: {}", i, alig.id, g, s, c);
                    visualize(
                        g.as_bytes(),
                        s.as_bytes(),
                        c.as_bytes(),
                        alig.offset as usize
                    );
                    println!();
                }
                */

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

                /*
                let ids = batch.ids();
                let sequences: Vec<Sequence> = batch.sequences().collect();
                for alig in alignments.alignments() {

                    let g = pipeline.guide.to_string();
                    let s = sequences[alig.id as usize].to_string();
                    let i = ids[alig.id as usize]; 
                    let c = alig.cigarx.to_string();

                    println!("id: {}/{}, g: {}, s: {}, c: {}", i, alig.id, g, s, c);
                    visualize(
                        g.as_bytes(),
                        s.as_bytes(),
                        c.as_bytes(),
                        alig.offset as usize
                    );
                    println!();
                }
                */

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

/// Thread responsible for managing a the CPU miner using all the remaining cores
#[tracing::instrument(skip_all)]
fn th_miner_cpu(pipeline: PipelineDescriptor, input: BatchRecv, output: AlignmentSend) {
    info!("started");

    let mut arena = Arena::alloc(1024 * 1024 * 1024);
    let mut mine_batch_idx = 0;
    while let Some(batch) = input.recv() {
        info!("mining batch #{mine_batch_idx}");
        mine_batch_idx += 1;

        arena.scoped(|memory| {
            /*
            for (id, seq) in batch.sequences_with_ids() {
                println!("seq: {:2}/{:?}", id, seq);
            }
            */

            /*
            let mut scores = memory.alloc_slice_fill(batch.len(), 255u8);
            batch.create_ed_mask(&pipeline.guide, &mut scores);

            let mut mask = memory.alloc_slice_fill(batch.len(), false);
            for (i, seq) in batch.sequences().enumerate() {
                if seq.mutation_score() <= pipeline.mutation_max {
                    mask[i] = scores[i] <= pipeline.thresholds.ed() as u8;
                }
            }

            batch.apply_mask(mask.as_ref(), true);

            let miner = SimpleMiner::new(&memory, &batch, &pipeline.guide, &pipeline.thresholds);

            // Mine all alignment in batches
            let mut alignments_memory = output.acquire();
            let mut alignments = AlignmentBatch::from_memory(&mut alignments_memory, false);
            for (_sidx, align) in miner {
                // Try to add alignment to current batch
                if !alignments.push(align) {
                    // Commit complete batch to writers
                    let alignment_count = alignments.len();
                    output.commit(alignments_memory, AlignmentBatchDescr { alignment_count });

                    // Acquire new alignment batch memory
                    alignments_memory = output.acquire();
                    alignments = AlignmentBatch::from_memory(&mut alignments_memory, false);
                    assert!(alignments.push(align));
                }
            }

            // Commit last probably not full batch
            let alignment_count = alignments.len();
            output.commit(alignments_memory, AlignmentBatchDescr { alignment_count });
            */

            let result = output.acquire();
            output.commit(result);

            // Signal that this sequence batch is now available
            input.finish(batch);
        });
    }

    info!("closed");
}
