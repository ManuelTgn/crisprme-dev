// modules used by the main function
mod bindings;
mod crispr;
mod utils;
mod memory;
mod alignment;
mod sequence;
mod batching;
mod storage;
mod engine;
mod error;
mod model;
mod pipeline;
mod annotation;
pub mod python;

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
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
    use std::array;

    use columnar::{MemoryPool, memory::CHUNK_SIZE, pipeline::{Driven, Pipeline, PipelineHandle}};
    use pyo3::{Bound, PyResult, Python, pyclass, pyfunction, pymethods, types::{PyAnyMethods, PyList}};

    use crate::{crispr::guide::Guide, model::{alignment::AlignmentFrame, input::{SeqBatch, SeqFrame, SeqOccFrame}}, pipeline::{merge::MergeJoinSorted, miner::Miner, resolve::Resolver, sink::NullSink, transform::PyTransform}, sequence::iupac::Iupac};

    /*
    use columnar::ext::pyo3::PyColumnView;

    use crate::{model::alignment::aligned::PyAlignmentBatch, pipeline::{PyPipeline, *}};
    use super::*;

    #[pymodule_init]
    fn _crisprme2_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
        
        // add the top-level function to the Python module
        // m.add_function(wrap_pyfunction!(extract_targets_rs, m)?)?;

        //m.add_function(wrap_pyfunction!(initialize_engine_logger, m)?)?;

        // Allows python to create a new pipeline
        m.add_function(wrap_pyfunction!(create_pipeline, m)?)?;

        m.add_class::<PyAlignmentBatch>()?;
        m.add_class::<PyColumnView>()?;
        m.add_class::<PyPipeline>()?;

        /*
        m.add_class::<TargetBatcher>()?;
        m.add_class::<FeedStatus>()?;
        m.add_class::<BatcherStats>()?;
        m.add_class::<HybridEngine>()?;
        m.add_class::<AlignmentParams>()?;
        m.add_class::<Thresholds>()?;
        m.add_class::<Guide>()?;
        m.add_class::<AlignmentBatchView>()?;
        */

        Ok(())
    }
    */

    pub use columnar::python::PyBuffer;
    pub use crate::pipeline::transform::PyAlignment;

    #[pyfunction]
    pub fn initialize_engine_logger() {
        tracing_subscriber::fmt()
                .compact()
                .with_target(false)
                .with_thread_ids(true)
                .with_max_level(tracing::Level::TRACE)
                .init();
    }

    #[pyclass]
    struct PyPipeline {

        // Pipeline memory pool
        pool: MemoryPool,

        // Input sender (Option so we can drop it explicitly to signal EOF)
        input: Driven<SeqBatch>,
        handle: PipelineHandle
    }

    #[pymethods]
    impl PyPipeline {

        fn send_debug_data(&mut self, py: Python<'_>) -> PyResult<()> {

            const ROWS: usize = 50_000;

            let iupacs: [Iupac; 4] = [
                Iupac::from_utf8('A'),
                Iupac::from_utf8('C'),
                Iupac::from_utf8('T'),
                Iupac::from_utf8('G')
            ];

            let mut seqs = SeqFrame::alloc(&self.pool, ROWS);
            let mut occs = SeqOccFrame::alloc(&self.pool, ROWS * 3);

            // Create debug sequences
            seqs.with_cols(|mut cols| {
                for (i, (id, content)) in 
                    cols.id.iter_mut().zip(cols.content.iter_mut()).enumerate() 
                {
                    *id = i as u32;
                    *content = array::from_fn(|j| 
                        iupacs[(i + j) % 4])
                }
            });

            // Create debug occurences
            occs.with_cols(|mut cols| {
                for (i, seq_id) in cols.seq_id.iter_mut().enumerate() {
                    *seq_id = (i / 3) as u32;
                }
            });

            // Release GIL while sending so pipeline workers can acquire it
            py.detach(|| {
                self.input.send(
                    SeqBatch { 
                        seq_len: 24, 
                        guide: Guide::new("GATTACA"), 
                        sequences: seqs, 
                        occurences: occs 
                    }
                ).unwrap();
            });

            Ok(())
        }

        /// Close the input and wait for all pipeline workers to finish.
        /// Must be called explicitly: dropping PyPipeline from Python will deadlock
        /// because worker threads need the GIL to call Python transforms.
        fn close(&mut self, py: Python<'_>) {
            self.input.close();
            py.detach(|| { // Release GIL so worker threads can finish their Python calls
                self.handle.join();
            });
        }
    }

    /// Create a driven pipeline with transforms
    #[pyfunction]
    fn pipeline<'py>(transforms: Bound<'py, PyList>) -> PyResult<PyPipeline> {

        let pool = MemoryPool::new(CHUNK_SIZE * 1000, |_, _| { });
        let (input, pipeline) = Pipeline::driven(10);

        let mut pipeline = pipeline
            .stage(2, |pool, _| Miner::new(pool))
            .stage(2, |pool, _| Resolver::new(pool))
            .stage(2, |pool, _| MergeJoinSorted::new(pool));

        // Add all transform stages
        println!("adding transform stages: ");
        for elem in transforms {
            println!("\t{:?}", elem.get_type()
                .getattr("__name__")
                .unwrap());

            let transform = elem.unbind();
            pipeline = pipeline.stage_once(|_| PyTransform::new(transform))
        }

        // Add sink stage
        let pipeline = pipeline.sink(2, |_, _| NullSink::<AlignmentFrame>::new());
        println!("pipeline ready!");

        let handle = pipeline.execute(&pool, 10);
        Ok(PyPipeline {
            handle,
            input,
            pool
        })
    }
}
