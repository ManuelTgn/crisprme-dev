use std::array;

use columnar::{memory::region_min_size, pipeline::{Emit, Stage, StageError}, python::PyBuffer};
use pyo3::{Py, PyAny, PyResult, Python, pyclass, pymethods};

use crate::model::alignment::AlignmentFrame;

#[pyclass(unsendable)]
pub struct PyAlignmentBatch {

    seq_id: PyBuffer,
    offset: PyBuffer,
    rguide: PyBuffer,
    rseq:   PyBuffer,

    features: [PyBuffer; 10],
    scores:   [PyBuffer; 4],
}

#[pymethods]
impl PyAlignmentBatch {

    fn seq_id(&self) -> PyResult<PyBuffer> { Ok(self.seq_id) }
    fn offset(&self) -> PyResult<PyBuffer> { Ok(self.offset) }
    fn rguide(&self) -> PyResult<PyBuffer> { Ok(self.rguide) }
    fn rseq(&self)   -> PyResult<PyBuffer> { Ok(self.rseq)   }

    fn feature(&self, idx: usize) -> PyResult<PyBuffer> { Ok(self.features[idx]) }
    fn score(&self, idx: usize)   -> PyResult<PyBuffer> { Ok(self.scores[idx])   }

}

/// Applies a transformation using a python callable
pub struct PyTransform(Py<PyAny>);
impl PyTransform {
    pub fn new(transform: Py<PyAny>) -> Self {
        Self(transform)
    }
}

impl Stage for PyTransform {
    type I = AlignmentFrame;
    type O = AlignmentFrame;

    fn name() -> &'static str { "PyTransform" }

    #[tracing::instrument(name = "pipeline:py_transform", skip_all)]
    fn process(&mut self, mut input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError> {
        input.with_cols(|mut cols| {

            // Size of smallest continous memory region
            let stride = region_min_size(&[
                cols.rguide.row_bytes(),
                cols.rseq.row_bytes(),
                cols.features.row_bytes(),
                cols.offset.row_bytes(),
                cols.scores.row_bytes()
            ]);

            let mut features: [_; 10] = cols.features.split();
            let mut scores:   [_;  4] = cols.scores.split();

            let mut row = 0;
            let total_rows = cols.occurence.rows();
            while row < total_rows {
                let len = stride.min(total_rows - row);

                // Get continous regions of memory, read-only
                
                let slice_seq_id = cols.seq_row_idx.slice(row, len);
                let slice_offset = cols.offset.slice(row, len);
                let slice_rguide = cols.rguide.slice(row, len);
                let slice_rseq   = cols.rseq.slice(row, len);

                // Get continous regions of memory, mutable

                let mut slice_features: Vec<_> = features.iter_mut().map(|s| s.slice_mut(row, len) ).collect();
                let mut slice_scores: Vec<_> = scores.iter_mut().map(|s| s.slice_mut(row, len) ).collect();

                // Create PyBuffer for all regions

                let seq_id = unsafe { PyBuffer::from_slice(slice_seq_id) };
                let offset = unsafe { PyBuffer::from_slice(slice_offset) };
                let rguide = unsafe { PyBuffer::from_array(slice_rguide) };
                let rseq   = unsafe { PyBuffer::from_array(slice_rseq)   };

                let features = array::from_fn(|i| unsafe { PyBuffer::from_slice_mut(slice_features[i]) });
                let scores   = array::from_fn(|i| unsafe { PyBuffer::from_slice_mut(slice_scores[i])   });

                tracing::debug!("running python transform");
                Python::attach(|py| {

                    let input = PyAlignmentBatch {
                        seq_id,
                        offset,
                        rguide,
                        rseq,
                        features,
                        scores,
                    };

                    self.0.call1(py, (input,))
                        .expect("python transform object caused an error");
                });

                row += len;
            }
        });

        // Forward
        emitter.emit(input)
    }
}
