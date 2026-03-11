use crate::annotation::annotate::annotate_batch;
use crate::annotation::features::FeatureRegistry;

use pyo3::prelude::*;
use pyo3::types::PyBytes;


/*
/// Python wrapper around FeatureRegistry
/// 
/// Exposes only safe operations needed for annotation
#[pyclass]
pub struct PyRegistry {
    inner: FeatureRegistry,
}

#[pymethods]
impl PyRegistry {

    /// Create a registry from a BED file
    /// 
    /// Args:
    ///     paths (str): Path to BED file
    /// 
    /// Raises:
    ///     IOError: If file cannot be read
    ///     ValueError: If BED file is malformed
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        // Construct feature registry from input BED file
        let registry = FeatureRegistry::from_bed(path)?;
        Ok(Self { inner: registry })
    }

    /// Return total number of unique features
    fn num_features(&self) -> usize {
        self.inner.num_features()
    }

    /// Return feature ID for a given feature name
    /// Returns None if feature is not present
    fn get_feature_id(&self, feature: &str) -> Option<usize> {
        self.inner.get_id(feature)
    }

    /// Annotate a batch of targets
    /// Args:
    ///     hits (List[List[int]]): Each inner list contains feature IDs 
    ///         overlapping a target
    /// 
    /// Returns:
    ///     List[bytes]: Compact byte arrays encoding feature overlaps
    /// 
    /// Raises:
    ///     ValueError: If feature IDs are invalid
    ///     ValueError: If input is empty
    fn annotate_batch<'py>(
        &self,
        py: Python<'py>,
        hits: Vec<Vec<usize>>
    ) -> PyResult<Vec<&'py PyBytes>> {
        let results = annotate_batch(
            hits, 
            self.inner.num_features()
        )?;

        Ok(results
            .into_iter()
            .map(|bits| PyBytes::new(py, bits.as_slice()))
            .collect())
    }
}
*/
