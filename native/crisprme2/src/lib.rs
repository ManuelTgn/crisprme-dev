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
mod partition;
mod pipeline;
pub mod python;
mod sequence;
mod storage;
mod utils;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::PyResult;

/// Defines the Python module structure and exposes Rust functions
#[pymodule]
pub mod _crisprme2_native {

    use crate::python::pylog::PyLoggerLayer;
    use tracing_subscriber::prelude::*;

    use std::{path::PathBuf, time::Instant};

    use columnar::{
        memory::CHUNK_SIZE,
        pipeline::{Driven, Pipeline, PipelineHandle},
        MemoryPool,
    };
    use itertools::izip;
    use pyo3::{
        exceptions::{PyOSError, PyValueError},
        pyclass, pyfunction, pymethods, pymodule,
        types::{PyAny, PyAnyMethods, PyList},
        Bound, Py, PyResult, Python,
    };

    use crate::{
        annotation::features::FeatureRegistry,
        bindings::cuda,
        crispr::pam::PAM,
        model::{
            alignment::AlignmentFrame,
            input::{SeqBatch, SeqFrame, SeqOccFrame, SEQ_MAX_LEN},
            occurence::Occurence,
        },
        partition::PrimaryCriteria,
        pipeline::{
            sink::{
                writer::{ContigLabels, PamContext, PamPlacement, TsvWriter, TsvWriterSink},
                NullSink,
            },
            source::reader::Reader,
            stage::{
                broadcast::Broadcast,
                miner::{GpuMiner, Miner},
                resolve::Resolver,
                transform::PyTransform,
            },
        },
        sequence::{iupac::Iupac, sequence::Sequence},
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
    pub use crate::crispr::guide::Guide;

    #[pymodule_export]
    pub use crate::alignment::thresholds::Thresholds;

    #[pymodule_export]
    pub use crate::python::pyannotation::PyRegistry;

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
                        pam_len: 0, // debug batches carry no PAM
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
                        pam_len: 0, // debug batches carry no PAM
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
            assert!(
                batcher.get_sequence_len() <= SEQ_MAX_LEN,
                "window sequence should fit inside a SeqFrame"
            );

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
                    batch
                        .occs
                        .iter()
                        .enumerate()
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
                        pam_len: batcher.get_pam_len(),
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
            tracing::info!(
                "pipeline took {:.2} s",
                self.started_at.elapsed().as_secs_f32()
            );
        }
    }

    /// Create a driven pipeline with transforms
    #[pyfunction]
    fn pipeline<'py>(
        chunks: usize,
        threshold: Thresholds,
        transforms: Bound<'py, PyList>,
        pam: &str,
        upstream: bool,
        outpath: PathBuf,
        contigs: Vec<String>,
        annotation_beds: Vec<String>,
        annotation_names: Vec<String>,
    ) -> PyResult<PyPipeline> {
        // Validate the PAM before allocating a multi-GB pool.
        let parsed_pam = PAM::new(pam)
            .map_err(|e| PyValueError::new_err(format!("invalid PAM {pam:?}: {e}")))?;
        let pam_ctx = PamContext::new(&parsed_pam, PamPlacement::from_upstream(upstream));
        tracing::info!(
            "guide column layout: {}",
            if upstream {
                "<PAM><guide>"
            } else {
                "<guide><PAM>"
            }
        );

        let contigs = ContigLabels::from_names(contigs)?;

        // Sink-side registries: built independently from the same BED paths the
        // annotation transform used. Deterministic bit assignment guarantees
        // decode here matches the transform's encode
        let annotation_registries: std::sync::Arc<[FeatureRegistry]> = annotation_beds
            .iter()
            .map(FeatureRegistry::from_bed)
            .collect::<Result<Vec<_>, _>>()? // AnnotationError -> PyErr
            .into();
        if annotation_beds.len() != annotation_names.len() {
            return Err(PyValueError::new_err(format!(
                "Mismatching number of annotation BEDs ({}) and annotation column names ({})",
                annotation_beds.len(),
                annotation_names.len()
            )));
        }

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
        let tsv_writer = TsvWriter::open(
            &outpath,
            pam_ctx,
            contigs,
            annotation_registries,
            &annotation_names,
        )
        .map_err(|e| {
            PyOSError::new_err(format!("cannot open TSV report {}: {e}", outpath.display()))
        })?;

        let pipeline = pipeline.sink(2, {
            let tsv_writer_clone = tsv_writer.clone();
            move |_, _| TsvWriterSink::new(&tsv_writer_clone)
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

    /// Split the intermediate report into primary + alternative reports.
    ///
    /// Clusters hits within `n` bases (strand-aware, keyed on contig/strand/start)
    /// and selects one primary per cluster under `criteria` — an ordered list of
    /// `(field, direction)` pairs, or `None` for the default (edit distance, then
    /// dna/rna bulges, then mismatches; all ascending). Returns
    /// `(clusters, primary, alternative)` counts.
    #[pyfunction]
    #[pyo3(signature = (intermediate, primary_out, alternative_out, n, criteria=None))]
    fn partition_report(
        intermediate: PathBuf,
        primary_out: PathBuf,
        alternative_out: PathBuf,
        n: u32,
        criteria: Option<Vec<String>>,
    ) -> PyResult<(usize, usize, usize)> {
        let criteria = match criteria {
            Some(spec) => PrimaryCriteria::from_spec(&spec).map_err(PyValueError::new_err)?,
            None => PrimaryCriteria::default(),
        };
        let stats = crate::partition::partition_report(
            &intermediate,
            &primary_out,
            &alternative_out,
            n,
            &criteria,
        )?;
        Ok((stats.clusters, stats.primary, stats.alternative))
    }

    /// Install the Rust -> Python logging bridge.
    ///
    /// Call this **once**, early, passing the `CrisprmeLoggers` bundle.
    /// It composes a compact stderr layer (dev console) with the
    /// [`PyLoggerLayer`], so every `tracing` event in the native core is
    /// mirrored into `basic.log` / `verbose.log` / `errors.log`.
    ///
    /// `TRACE` is filtered out to keep hot-path `trace!` events off the GIL.
    #[pyfunction]
    fn init_logging(loggers: &Bound<'_, PyAny>) -> PyResult<bool> {
        use tracing_subscriber::prelude::*;
        use tracing_subscriber::{filter::LevelFilter, fmt};

        let py_layer = crate::python::pylog::PyLoggerLayer::from_bundle(loggers)?;

        let installed = tracing_subscriber::registry()
            .with(LevelFilter::DEBUG)
            .with(py_layer)
            .try_init()
            .is_ok();

        Ok(installed) // report install status instead of hiding it
    }

    /// Concrete PAM variants in `pam_id` order: entry i is the ASCII PAM whose
    /// `Occurence::pam_id() == i`. Lets Python build a weight table aligned to the
    /// `pam_id` column.
    #[pyfunction]
    fn pam_variants_ascii(motif: &str) -> PyResult<Vec<String>> {
        let pam =
            PAM::new(motif).map_err(|e| PyValueError::new_err(format!("Invalid PAM: {e}")))?;
        (0..pam.variant_count() as u16)
            .map(|id| {
                pam.pam_variant_ascii(id)
                    .map_err(|e| PyValueError::new_err(format!("variant {id}: {e}")))
            })
            .collect()
    }
}
