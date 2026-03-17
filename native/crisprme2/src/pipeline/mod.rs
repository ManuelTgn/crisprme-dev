use std::sync::Arc;

pub mod merge;
pub mod resolve;
pub mod transform;
pub mod miner;

pub mod sink;

/*
#[pyclass]
pub struct PyPipeline {

    pub inseq_tx: Option<ConnectorTx<BatchRef<SeqSchema, SeqBatchMetadata>>>,
    pub mutat_rx: ConnectorRx<BatchMut<AlignmentSchema, ()>>,

    pub sequences:  Arc<Pool<SeqSchema>>,
    pub occurences: Arc<Pool<SeqOccSchema>>,

    // NOTE: Must be declared at the end for drop order
    pub pipeline: Pipeline<()>
}

#[pymethods]
impl PyPipeline {

    fn close(&mut self) {
        self.inseq_tx.take();
    }

    // Example submit
    fn submit_example(&mut self, py: Python<'_>) {
        println!("Submitting new sequence batch...");

        use crate::model::input::sequences::schema  as ss;
        use crate::model::input::occurences::schema as os;

        let mut seq = self.sequences.acquire()
            .unwrap();

        let sequence = Sequence::from_ascii_lossy("GATTACA");
        let guide = Guide::from_ascii_bytes_lossy(b"GAT");

        println!("creating input batch with 10 rows:");
        seq.set_len(10);
        seq.mutate(
            (ss::id,),
            |(ids,)| {
                for i in 0..10 {
                    println!("\tid[{i}] = {i}");
                    ids[i] = i as u32;
                }
            }
        );

        // Create 4 occurence batches
        println!("Creating occurence batches:");
        let occurences: Vec<_> = (1..5).into_iter()
            .map(|i| {

                let mut occ = self.occurences.acquire()
                    .unwrap();

                let rows = i * 2;
                occ.set_len(rows);

                println!("\tbatch {i} with {rows} rows");
                occ.mutate(
                    (os::seq_id, os::occurence), 
                    |(seq_ids, occurences)| {
                        for j in 0..rows {
                            
                            let this_occ = occurence(i as u32, j as u32, i as u8);
                            println!("\t\tseq_id[{j}] = {i}, occurences[{j}] = {this_occ}");

                            occurences[j] = this_occ;
                            seq_ids[j] = i as u32;
                        } 
                    }
                );

                occ.freeze()
            })
            .collect();

        println!("composing batches, occurences: {}", occurences.len());
        let seq = seq.with_metadata(SeqBatchMetadata { 
            seq_len: sequence.len(),
            occurences, 
            guide
        }).freeze();

        // Release GIL so stage threads can acquire it
        py.detach(|| {
            self.inseq_tx.as_ref().unwrap().send(seq)
                .unwrap();
        });

        println!("Submission done!")
    }

    // Example receive
    fn receive(&mut self, py: Python<'_>) -> PyResult<(bool, PyAlignmentBatch)> {
        println!("Receiving alignment batch...");

        let result = py.detach(|| {
            self.mutat_rx.recv()
        });

        match result {
            Ok(result) => Ok((false, PyAlignmentBatch { batch: Some(result) })),
            Err(e) => {
                eprint!("error: {e}");
                Ok((true, PyAlignmentBatch { batch: None } ))
            }
        }
    }

    fn wait(&mut self) {
        self.pipeline.wait();
    }
}

fn make_pool<S: Schema>(slots: usize, elements: usize) -> Arc<Pool<S>> {
    Arc::new(Pool::new(slots, elements))
}

#[pyfunction]
pub fn create_pipeline(py: Python<'_>, transform: Py<PyAny>) -> PyResult<PyPipeline> {

    // ------ Pools -----------------------------------------------------------

    let sequences  = make_pool::<SeqSchema>(8, 16);
    let occurences = make_pool::<SeqOccSchema>(8, 16);
    let mined      = make_pool::<MinedSchema>(8, 16);
    let resolved   = make_pool::<ResolvedSchema>(8, 16);
    let aligned    = make_pool::<AlignmentSchema>(8, 16);

    // ------ Connectors ------------------------------------------------------

    let (inseq_tx, inseq_rx) = connector_ref::<SeqSchema, SeqBatchMetadata>(12);
    let (mined_tx, mined_rx) = connector_mut::<MinedSchema, MinedBatchMetadata>(12);
    let (rslvd_tx, rslvs_rx) = connector_mut::<ResolvedSchema, ResolvedBatchMetadata>(12);
    let (align_tx, align_rx) = connector_mut::<AlignmentSchema, ()>(12);
    let (mutat_tx, mutat_rx) = connector_mut::<AlignmentSchema, ()>(12);

    // ------ Stages ----------------------------------------------------------

    let mut pipeline = Pipeline::new(());
    py.detach(|| {

        pipeline.stage("mine", 4, inseq_rx, mined_tx, move |_ctx| {
            MineScanner::new(mined.clone())
        });

        pipeline.stage("resolve", 4, mined_rx, rslvd_tx, move |_ctx| {
            AlignmentSimpleResolve::new(resolved.clone(), 1024 * 1024 * 10)
        });

        pipeline.stage("broadcast", 4, rslvs_rx, align_tx, move |_ctx| {
            AlignmentBroadcast::new(aligned.clone(), 1024 *1024 * 20)
        });
    });
    
    pipeline.stage("transform", 4, align_rx, mutat_tx, move |_ctx| {
        let t = Python::try_attach(|py| transform.clone_ref(py));
        AlignmentPythonTransform::new(t.unwrap())
    });

    Ok(PyPipeline {
        inseq_tx: Some(inseq_tx),
        mutat_rx,
        sequences,
        occurences,
        pipeline
    })
}
 */