use pyo3::{pyclass, pymethods};

/// Threholds for filtering and alignment processes
#[derive(Debug, Clone, Copy)]
#[pyclass]
pub struct Thresholds {
    /// Max allowed gaps in guide
    #[pyo3(get, set)]
    pub qgap: u32,
    /// Max allowed gaps in target
    #[pyo3(get, set)]
    pub tgap: u32,
    /// Max allowed mismatches
    #[pyo3(get, set)]
    pub mism: u32,
}

impl Thresholds {
    /// Calculate the max edit distance based on the thresholds
    pub fn ed(&self) -> u32 {
        self.qgap + self.tgap + self.mism
    }
}

#[pymethods]
impl Thresholds {
    #[new]
    pub fn new(qgap: u32, tgap: u32, mism: u32) -> Self {
        Self { qgap, tgap, mism }
    }
}
