
/// Threholds for the filtering and mining processes
#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    /// Max allowed gaps in query
    pub qgap: u32,
    /// Max allowed gaps in target
    pub tgap: u32,
    /// Max allowed mismatches
    pub mism: u32,
}

impl Thresholds {
    /// Calculate the max edit distance based on the thresholds
    pub fn ed(&self) -> u32 {
        self.qgap + self.tgap + self.mism
    }
}
