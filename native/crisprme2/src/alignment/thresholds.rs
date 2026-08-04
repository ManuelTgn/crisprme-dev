use pyo3::{pyclass, pymethods};

/// Threholds for filtering and alignment processes
#[derive(Debug, Clone, Copy)]
#[pyclass]
pub struct Thresholds {
    /// Maximum **RNA bulges** ... Python: `brna`.
    #[pyo3(get, set)]
    pub qgap: u32,
    /// Maximum **DNA bulges** ... Python: `bdna`.
    #[pyo3(get, set)]
    pub tgap: u32,
    /// Maximum mismatches.
    #[pyo3(get, set)]
    pub mism: u32,
    /// Effective (normalized) edit-distance cap: `mism + qgap + tgap` when the
    /// user passes 0, otherwise clamped to that sum
    #[pyo3(get)]
    pub max_ed: u32,
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
    pub fn new(qgap: u32, tgap: u32, mism: u32, max_ed: u32) -> Self {
        let sum = qgap + tgap + mism;
        // 0 == "disabled" -> use the sum; otherwise a value above the sum can
        // never bind, so clamp it
        let max_ed = if max_ed == 0 { sum } else { max_ed.min(sum) };
        Self {
            qgap,
            tgap,
            mism,
            max_ed,
        }
    }
}
