use crate::error::crisprme_errors::AnnotationError;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;




/// Immutable registry mapping feature names to stable numeric IDs.
/// 
/// This structure is built once from a BED file and then shared across threads
/// for fast annotation of millions of CRISPR off-targets.
/// 
/// The registry guarantees:
/// - Deterministic feature indexing
/// - O(1) lookup time
/// - Thread safety (no mutability)
#[derive(Debug)]
pub struct FeatureRegistry {
    features: Vec<String>,
    feature_to_id: HashMap<String, usize>,
}


impl FeatureRegistry {

    /// Build a registry form a BED file.
    /// 
    /// The 4th column is interpreted as the feature name
    /// 
    /// Returns:
    /// - `Ok(Self)` on success
    /// - `AnnotationError` on I/O or malformed file 
    pub fn from_bed<P: AsRef<Path>>(path: P) -> Result<Self, AnnotationError> {
        // Open input BED file and initialize stdin reader
        let file = File::open(path)
            .map_err(|e| AnnotationError::IoError(e.to_string()))?;
        let reader = BufReader::new(file);

        // Initialize features list and features to id map
        let mut features = Vec::new();
        let mut features_to_id = HashMap::new();

        // Read each line in BED file to retrieve features set
        for (line_number, line) in reader.lines().enumerate() {
            let line = line
                .map_err(|e| AnnotationError::IoError(e.to_string()))?;

            if line.starts_with('#') { 
                continue;  // Skip comments
            }

            // Parse fields in BED lines
            let fields: Vec<&str> = line.split('\t').collect();  
            if fields.len() < 4 {
                return Err(AnnotationError::MalformedBed { line: line_number + 1 });
            }

            // Add feature to features set for current BED file
            let feature = fields[3];
            if !features_to_id.contains_key(feature) {
                let id = features.len();
                features.push(feature.to_string());
                features_to_id.insert(feature.to_string(), id);
            }
        }

        Ok( Self { features, feature_to_id })
    }

    #[inline(always)]
    pub fn num_features(&self) -> usize {
        self.features.len()
    }

    /// Return the numeric ID associated with a feature name
    /// 
    /// Returns:
    /// - `Some(id)` if the feature exists
    /// - `None` if the feature is not present in the registry
    /// 
    /// Notes:
    /// - This method performs O(1) average-time lookup
    /// - No string allocation is performed 
    /// - Safe for concurrent reads (registry is immutable after construction)
    #[inline(always)]
    pub fn get_id(&self, feature: &str) -> Option<usize> {
        self.feature_to_id.get(feature).copied()
    }
}


