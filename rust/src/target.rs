use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]

pub struct Target {
    #[pyo3(get)]
    pub contig: String,
    #[pyo3(get)]
    pub position: usize,
    #[pyo3(get)]
    pub orientation: bool,  // 0 or 1
    #[pyo3(get)]
    pub target: String,
}

impl Target {
    pub fn new(contig: &str, pos: usize, ori: bool, target_seq: &str) -> Self {
        Self {
            contig: contig.to_string(),
            position: pos,
            orientation: ori,
            target: target_seq.to_string(),
        }
    }
}