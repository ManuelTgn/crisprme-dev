use std::time::Instant;

use columnar::{pipeline::{Emit, Stage, StageError}, pool::BatchMut};
use pyo3::{Py, PyAny, Python};

use crate::model::alignment::{AlignmentSchema, aligned::PyAlignmentBatch};

/// Modifies an alignment batch using python
pub struct AlignmentPythonTransform {
    transform: Py<PyAny>,
}

impl AlignmentPythonTransform {
    pub fn new(transform: Py<PyAny>) -> Self {
        Self { transform }
    }
}

impl Stage for AlignmentPythonTransform {

    type Input  = BatchMut<AlignmentSchema, ()>;
    type Output = BatchMut<AlignmentSchema, ()>;

    fn process<E>(&mut self, input: Self::Input, emitter: &mut E) -> Result<(), StageError>
    where
        E: Emit<Self::Output> 
    {
        let _span = tracing::debug_span!("alignment-python-transform")
            .entered();

        let py_batch = PyAlignmentBatch { 
            batch: Some(input) 
        };

        let start = Instant::now();
        let result = Python::attach(|py| {
                let py_batch = Py::new(py, py_batch)
                    .expect("unable to attach buffer to python");

                self.transform.call1(py, (&py_batch,))
                    .expect("unable to call transform function on batch");

                // Take the batch back after callback returns
                let mut inner = py_batch.borrow_mut(py);
                inner.batch.take().unwrap()
            });

        tracing::debug!("python took {:.2} (s) to process {} rows",
            start.elapsed().as_secs_f32(), result.len());

        emitter.emit(result)
    }
}