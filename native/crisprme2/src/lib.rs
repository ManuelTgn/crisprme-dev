// modules used by the main function
mod alignment;
mod annotation;
mod batching;
mod bindings;
mod crispr;
mod engine;
mod error;
mod memory;
mod model;
mod pipeline;
pub mod python;
mod sequence;
mod storage;
mod utils;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::PyResult;

/// Finds all potential target candidates (CRISPR gRNAs) within a given sequence.
///
/// This function converts the input sequence and PAM into IUPAC bitmasks and performs a
/// parallelized scan to identify all positions where the target sequence and its
/// associated PAM match the defined criteria.
///
/// # Arguments
/// * `sequence` (str): The large DNA/RNA sequence to be scanned (e.g., a contig or chromosome).
/// * `contig` (str): The name/identifier of the sequence (e.g., "chr1").
/// * `pam_seq` (str): The Protospacer Adjacent Motif (PAM) sequence (e.g., "NGG").
/// * `k` (usize): The length of the target/protospacer sequence, excluding the PAM.
/// * `right` (bool): If `true`, the PAM is expected to be immediately *right* of the target sequence.
///                   If `false`, the PAM is expected to be immediately *left* of the target sequence.
/// * `threads` (usize): The number of threads to use for parallel scanning.
///
/// # Returns
/// A `list` of `target::Target` objects, where each object contains the position,
/// orientation, and bitmask sequence of a found target.
///
/// # Errors
/// Returns a `PyValueError` if input constraints are violated (e.g., invalid sizes or PAM sequence).
#[pyfunction]
pub fn extract_targets_rs(
    sequence: &str,
    pam_seq: &str,
    size: usize,
    right: bool,
    threads: usize,
) -> PyResult<(Vec<usize>, Vec<u8>)> {
    let pat = crispr::pam::ParsedPAM::new(pam_seq)
        .map_err(|e| PyErr::new::<PyValueError, _>(format!("Invalid PAM sequence: {e}")))?;

    // Execute the core parallel scanning logic and return the results
    sequence::scanner::scan_targets(sequence, &pat, size, right, threads)
        .map_err(|e| PyErr::new::<PyValueError, _>(e))
}

/// Defines the Python module structure and exposes Rust functions
#[pymodule]
pub mod _crisprme2_native {

    use std::{path::PathBuf, time::Instant};

    use columnar::{
        memory::CHUNK_SIZE,
        pipeline::{Driven, Pipeline, PipelineHandle},
        MemoryPool,
    };
    use itertools::izip;
    use pyo3::{
        Bound, Py, PyResult, Python, pyclass, pyfunction, pymethods, pymodule, types::{PyAnyMethods, PyList}
    };

    use crate::{
        bindings::cuda, model::{
            alignment::AlignmentFrame,
            input::{SEQ_MAX_LEN, SeqBatch, SeqFrame, SeqOccFrame}, occurence::Occurence,
        }, pipeline::{
            sink::{NullSink, writer::{CsvWriter, CsvWriterSink}}, source::reader::Reader, stage::{broadcast::Broadcast, miner::{GpuMiner, Miner}, resolve::Resolver, transform::PyTransform}
        }, sequence::{iupac::Iupac, sequence::Sequence}
    };

    #[pymodule_export]
    pub use columnar::python::PyBuffer;

    #[pymodule_export]
    pub use crate::batching::batching::TargetBatcher;

    #[pymodule_export]
    pub use crate::batching::batching::BatcherStats;

    #[pymodule_export]
    pub use crate::batching::batching::FeedStatus;

    #[pymodule_export]
    pub use crate::pipeline::stage::transform::PyAlignmentBatch;

    #[pymodule_export]
    pub use crate::extract_targets_rs;

    #[pymodule_export]
    pub use crate::crispr::guide::Guide;

    #[pymodule_export]
    pub use crate::alignment::thresholds::Thresholds;


    #[pyfunction]
    pub fn init_tracing() {
        tracing_subscriber::fmt()
            //.compact()
            .with_target(false)
            .with_file(false)
            .with_thread_ids(false)
            .with_max_level(tracing::Level::INFO)
            .init();
    }

    #[pyclass]
    struct PyPipeline {
        // Pipeline memory pool
        pool: MemoryPool,

        threshold: Thresholds,

        // Input sender (Option so we can drop it explicitly to signal EOF)
        input: Driven<SeqBatch>,
        handle: PipelineHandle,
    }

    #[pymethods]
    impl PyPipeline {

        fn send_debug_minable_data(&mut self, py: Python<'_>) -> PyResult<()> {

            const ROWS: usize = 500;

            let mut seqs = SeqFrame::alloc(&self.pool, ROWS);
            let mut occs = SeqOccFrame::alloc(&self.pool, ROWS * 3);

            // Create debug sequence
            let sequence = Sequence::from_utf8("GATTACAGATTACA");
            seqs.with_cols(|mut cols| {
                for content in cols.content.iter_mut() {
                    for j in 0..sequence.len() {
                        content[j] = sequence[j];
                    }
                }
            });

            // Create debug occurences
            occs.with_cols(|mut cols| {
                for (i, seq_idx) in cols.seq_row_idx.iter_mut().enumerate() {
                    *seq_idx = (i % 3) as u32;
                }
            });

            // Release GIL while sending so pipeline workers can acquire it
            py.detach(|| {
                self.input
                    .send(SeqBatch {
                        thresholds: Thresholds::new(1, 1, 2),
                        seq_len: sequence.len(),
                        guide: Guide::new("GATTACA"),
                        sequences: seqs,
                        occurences: occs,
                    })
                    .unwrap();
            });

            Ok(())
        }

        fn send_debug_data(&mut self, py: Python<'_>) -> PyResult<()> {

            const ROWS: usize = 10;

            let seq_len: usize = 24;
            let iupacs: [Iupac; 4] = [
                Iupac::from_utf8('A'),
                Iupac::from_utf8('C'),
                Iupac::from_utf8('T'),
                Iupac::from_utf8('G'),
            ];

            let mut seqs = SeqFrame::alloc(&self.pool, ROWS);
            let mut occs = SeqOccFrame::alloc(&self.pool, ROWS * 3);

            // Create debug sequences
            seqs.with_cols(|mut cols| {
                for (i, content) in cols.content.iter_mut().enumerate() {
                    for j in 0..seq_len {
                        content[j] = iupacs[(i + j) % 4];
                    }
                }
            });

            // Create debug occurences
            occs.with_cols(|mut cols| {
                for (i, seq_idx) in cols.seq_row_idx.iter_mut().enumerate() {
                    *seq_idx = (i % ROWS) as u32;
                }
            });

            // Release GIL while sending so pipeline workers can acquire it
            py.detach(|| {
                self.input
                    .send(SeqBatch {
                        thresholds: Thresholds::new(1, 1, 2),
                        seq_len,
                        guide: Guide::new("GATTACAGATTACA"),
                        sequences: seqs,
                        occurences: occs,
                    })
                    .unwrap();
            });

            Ok(())
        }

        /// Submit the content of a TargetBatcher
        pub fn submit(&mut self, py: Python<'_>, batcher: &mut TargetBatcher) -> PyResult<()> {

            assert!(batcher.get_sequence_len() <= SEQ_MAX_LEN,
                "window sequence should fit inside a SeqFrame");

            // Create compact representation
            let batch = batcher.flush_to_batch();

            // Copy sequences
            let mut seqs = SeqFrame::alloc(&self.pool, batch.len());
            seqs.with_cols(|mut cols| {
                for (i, content) in cols.content.iter_mut().enumerate() {
                    // Copy content to frame
                    let window = &batch.windows[i];
                    for j in 0..window.len() {
                        content[j] = Iupac::new(window[j]);
                    }
                }
            });

            // Copy occurences
            let total_occs = batch.occs.iter().map(|o| o.len()).sum();
            let mut occs = SeqOccFrame::alloc(&self.pool, total_occs);
            occs.with_cols(|mut cols| {
                // Each occurrence carries the index of the WINDOW (source sequence) it
                // belongs to, so seq_row_idx < source_seq_count (Broadcast/Reader contract).
                let iter = izip!(
                    cols.seq_row_idx.iter_mut(),
                    cols.occurence.iter_mut(),
                    batch.occs.iter().enumerate()
                        .flat_map(|(w, s)| s.iter().map(move |occ| (w as u32, *occ))),
                );
                for (dst_seq_id, dst_occ, (w, src_occ)) in iter {
                    *dst_seq_id = w;
                    *dst_occ = Occurence(src_occ);
                }
            });

            // Release GIL while sending so pipeline workers can acquire it
            py.detach(|| {
                self.input
                    .send(SeqBatch {
                        thresholds: self.threshold.clone(),
                        seq_len: batcher.get_sequence_len(),
                        guide: batcher.get_guide(),
                        sequences: seqs,
                        occurences: occs,
                    })
                    .unwrap();
            });

            Ok(())
        }

        /// Close the input and wait for all pipeline workers to finish.
        /// Must be called explicitly: dropping PyPipeline from Python will deadlock
        /// because worker threads need the GIL to call Python transforms.
        fn close(&mut self, py: Python<'_>) {
            self.input.close();
            py.detach(|| {
                // Release GIL so worker threads can finish their Python calls
                self.handle.join();
            });
        }
    }

    /// Similar to PyPipeline but with a source stage that reads batches of sequences from disk.
    #[pyclass]
    struct PySourcedPipeline {
        started_at: Instant,
        handle: PipelineHandle,
        pool: MemoryPool,
    }

    impl Drop for PySourcedPipeline {
        fn drop(&mut self) {
            tracing::info!("pipeline took {:.2} s", 
                self.started_at.elapsed().as_secs_f32());
        }
    }

    /// Create a driven pipeline with transforms
    #[pyfunction]
    fn pipeline<'py>(
        chunks: usize, 
        threshold: Thresholds,
        transforms: Bound<'py, PyList>
    ) -> PyResult<PyPipeline> {
        
        // Create memory pool and pin all chunks for DMA from GPU
        let pool = MemoryPool::new(CHUNK_SIZE * chunks, |ptr, bytes| {
            tracing::trace!("pinning chunk (ptr = {:?}, bytes = {})", ptr, bytes);
            cuda::pin(ptr, bytes);
        });
        
        tracing::info!("building pipeline...");
        let (input, pipeline) = Pipeline::driven(10);

        let mut pipeline = pipeline
            .stage(1, |pool, _| GpuMiner::new(pool, 100_000, 32, 100_000, 0))
            .stage(2, |pool, _| Resolver::new(pool))
            .stage(2, |pool, _| Broadcast::new(pool));

        // Add all transform stages
        tracing::info!("adding transform stages: ");
        for elem in transforms {
            tracing::info!("\t{:?}", elem.get_type().getattr("__name__").unwrap());

            let transform = elem.unbind();
            pipeline = pipeline.stage_once(|_| PyTransform::new(transform))
        }

        // Add sink stage
        //let pipeline = pipeline.sink(2, |_, _| NullSink::<AlignmentFrame>::new());
        let csv_writer = CsvWriter::open("results.csv".into());
        let pipeline = pipeline.sink(2, {
            let csv_writer_clone = csv_writer.clone();
            move |_, _| { 
                CsvWriterSink::new(&csv_writer_clone)
            }
        });
        
        tracing::info!("pipeline ready!");
        let handle = pipeline.execute(&pool, 3);
        Ok(PyPipeline {
            threshold,
            handle,
            input,
            pool,
        })
    }


    /// Create a dataset pipeline that reads batches of sequences from disk, applies transforms, and writes results to disk.
    #[pyfunction]
    fn dataset_pipeline<'py>(
        chunks: usize, 
        transforms: Bound<'py, PyList>, 
        folder: PathBuf,
        batch_size: usize,
        guide: Guide,
        thresholds: Thresholds,
        sequence_len: usize,
    ) -> PyResult<PySourcedPipeline> {

        // Create memory pool and pin all chunks for DMA from GPU
        let pool = MemoryPool::new(CHUNK_SIZE * chunks, |ptr, bytes| {
            tracing::trace!("pinning chunk (ptr = {:?}, bytes = {})", ptr, bytes);
            cuda::pin(ptr, bytes);
        });

        let seq_path = folder.join("sequences.bin");
        let pos_path = folder.join("positions.bin");

        tracing::info!("building pipeline...");
        let pipeline = Pipeline::source(1, move |pool, _| {
                Reader::open(seq_path.clone(), pos_path.clone(), sequence_len, batch_size, guide.clone(), thresholds, pool.clone())
                    .expect("unable to create reader")
        });

        let mut pipeline = pipeline
            .stage(1, |pool, _| GpuMiner::new(pool, 100_000, 32, 100_000, 0))
            .stage(2, |pool, _| Resolver::new(pool))
            .stage(2, |pool, _| Broadcast::new(pool));

        // Add all transform stages
        tracing::info!("adding transform stages: ");
        for elem in transforms {
            tracing::info!("\t{:?}", elem.get_type().getattr("__name__").unwrap());

            let transform = elem.unbind();
            pipeline = pipeline.stage_once(|_| PyTransform::new(transform))
        }

        // Add sink stage
        let csv_writer = CsvWriter::open("results.csv".into());
        let pipeline = pipeline.sink(2, {
            let csv_writer_clone = csv_writer.clone();
            move |_, _| { 
                CsvWriterSink::new(&csv_writer_clone)
            }
        });

        tracing::info!("pipeline ready!");
        let handle = pipeline.execute(&pool, 3);
        Ok(PySourcedPipeline {
            started_at: Instant::now(),
            handle,
            pool,
        })
    }

    /// Install the Rust -> Python logging bridge.
    ///
    /// Call this **once**, early, passing the `CrisprmeLoggers` bundle.
    /// It composes a compact stderr layer (dev console) with the
    /// [`PyLoggerLayer`], so every `tracing` event in the native core is
    /// mirrored into `basic.log` / `verbose.log` / `errors.log`.
    ///
    /// `TRACE` is filtered out to keep hot-path `trace!` events off the GIL.
    /// Uses `try_init`, so a second call (or a prior `init_tracing`) is a
    /// silent no-op rather than a panic — but you should call only one of
    /// the two initialisers
    #[pyfunction]
    fn init_logging(loggers: &Bound<'_, PyAny>) -> PyResult<()> {
        use tracing_subscriber::prelude::*;
        use tracing_subscriber::{filter::LevelFilter, fmt};

        let py_layer = crate::python::pylog::PyLoggerLayer::from_bundle(loggers)?;

        tracing_subscriber::registry()
            .with(LevelFilter::DEBUG)
            .with(fmt::layer().compact().with_target(false))
            .with(py_layer)
            .try_init()
            .ok();

        Ok(())
    }
}
