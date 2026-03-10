use std::sync::Arc;

use columnar::{buffer::Schema, pipeline::Pipeline, pool::{BatchMut, BatchRef, ConnectorRx, ConnectorTx, Pool, connector_mut, connector_ref}};
use pyo3::{Py, PyAny, PyResult, Python, pyclass, pyfunction, pymethods};

use crate::{crispr::guide::Guide, model::{alignment::{AlignmentSchema, MinedBatchMetadata, MinedSchema, ResolvedBatchMetadata, ResolvedSchema, aligned::PyAlignmentBatch}, input::{SeqBatchMetadata, SeqOccSchema, SeqSchema, occurence}}, pipeline::{broadcast::AlignmentBroadcast, resolve::AlignmentSimpleResolve, scanner::MineScanner, transform::AlignmentPythonTransform}, sequence::sequence::Sequence};


pub mod broadcast;
pub mod resolve;
pub mod transform;
pub mod scanner;

#[pyclass]
pub struct PyPipeline {

    pub inseq_tx: ConnectorTx<BatchRef<SeqSchema, SeqBatchMetadata>>,
    pub mutat_rx: ConnectorRx<BatchMut<AlignmentSchema, ()>>,

    pub sequences:  Arc<Pool<SeqSchema>>,
    pub occurences: Arc<Pool<SeqOccSchema>>,

    // NOTE: Must be declared at the end for drop order
    pub pipeline: Pipeline<()>
}

#[pymethods]
impl PyPipeline {

    // Example submit
    fn submit(&mut self, py: Python<'_>) {
        println!("Submitting new sequence batch...");

        use crate::model::input::sequences::schema  as ss;
        use crate::model::input::occurences::schema as os;

        let mut seq = self.sequences.acquire()
            .unwrap();

        let sequence = Sequence::from_ascii_lossy("GATTACA");
        let guide = Guide::from_ascii_bytes_lossy(b"GAT");

        seq.set_len(10);
        seq.mutate(
            (ss::id,),
            |(ids,)| {
                for i in 0..10 {
                    ids[i] = i as u32;
                }
            }
        );

        // Create 4 occurence batches
        let occurences: Vec<_> = (1..5).into_iter()
            .map(|i| {

                let mut occ = self.occurences.acquire()
                    .unwrap();

                occ.set_len(i * 2);
                occ.mutate(
                    (os::seq_id, os::occurence), 
                    |(seq_ids, occurences)| {
                        for j in 0..i {
                            occurences[i] = occurence(i as u32, i as u32, i as u8);
                            seq_ids[j] = i as u32;
                        } 
                    }
                );

                occ.freeze()
            })
            .collect();

        let seq = seq.with_metadata(SeqBatchMetadata { 
            seq_len: sequence.len(),
            occurences, 
            guide
        }).freeze();

        // Release GIL so stage threads can acquire it
        py.detach(|| {
            self.inseq_tx.send(seq)
                .unwrap();
        });

        println!("Submission done!")
    }

    // Example receive
    fn receive(&mut self, py: Python<'_>) -> PyResult<PyAlignmentBatch> {
        println!("Receiving alignment batch...");

        let result = py.detach(|| {
            self.mutat_rx.recv()
                .unwrap()
        });

        Ok(PyAlignmentBatch { 
            batch: Some(result) 
        })
    }

    fn wait(&mut self) {
        self.pipeline.wait();
    }
}

fn make_pool<S: Schema>(slots: usize, elements: usize) -> Arc<Pool<S>> {
    Arc::new(Pool::new(slots, elements))
}

#[pyfunction]
pub fn create_pipeline<'py>(transform: Py<PyAny>) -> PyResult<PyPipeline> {

    Python::initialize();

    // ------ Pools -----------------------------------------------------------

    let sequences  = make_pool::<SeqSchema>(2, 16);
    let occurences = make_pool::<SeqOccSchema>(8, 16);
    let mined      = make_pool::<MinedSchema>(2, 16);
    let resolved   = make_pool::<ResolvedSchema>(2, 16);
    let aligned    = make_pool::<AlignmentSchema>(2, 16);

    // ------ Connectors ------------------------------------------------------

    let (inseq_tx, inseq_rx) = connector_ref::<SeqSchema, SeqBatchMetadata>(2);
    let (mined_tx, mined_rx) = connector_mut::<MinedSchema, MinedBatchMetadata>(2);
    let (rslvd_tx, rslvs_rx) = connector_mut::<ResolvedSchema, ResolvedBatchMetadata>(2);
    let (align_tx, align_rx) = connector_mut::<AlignmentSchema, ()>(2);
    let (mutat_tx, mutat_rx) = connector_mut::<AlignmentSchema, ()>(2);

    // ------ Stages ----------------------------------------------------------

    let mut pipeline = Pipeline::new(());

    pipeline.stage("mine", 1, inseq_rx, mined_tx, move |_ctx| {
        MineScanner::new(mined.clone())
    });

    pipeline.stage("resolve", 1, mined_rx, rslvd_tx, move |_ctx| {
        AlignmentSimpleResolve::new(resolved.clone(), 1024 * 1024 * 10)
    });

    pipeline.stage("broadcast", 1, rslvs_rx, align_tx, move |_ctx| {
        AlignmentBroadcast::new(aligned.clone(), 1024 *1024 * 20)
    });

    pipeline.stage("transform", 1, align_rx, mutat_tx, move |_ctx| {
        let t = Python::try_attach(|py| transform.clone_ref(py));
        AlignmentPythonTransform::new(t.unwrap())
    });

    Ok(PyPipeline { 
        inseq_tx, mutat_rx, sequences, occurences, pipeline 
    })
}