use crate::annotation::bitset::AnnotationBits;
use crate::error::crisprme_errors::AnnotationError;

use rayon::{prelude::*, result};

/// Annotate a batch of targets in parallel.
/// 
/// Each element of `hits_per_target` represents one target.
/// - `hits_per_target[i]` contains the feature IDs overlapping target i
/// - An empty inner vector means no overlap
/// 
/// Returns:
/// - Vector of AnnotationBits, one per target
/// 
/// Errors:
/// - EmptyInput if the outer vector is empty
/// - Invalid FeatureId if any feature id is out of range
pub fn annotate_batch(
    hits_per_target: Vec<Vec<usize>>, 
    num_features: usize
) -> Result<Vec<AnnotationBits>, AnnotationError> {
    if hits_per_target.is_empty() {
        return Err(AnnotationError::EmptyInput);
    }

    // Proceed with bits flipping
    hits_per_target
        .into_par_iter()
        .map(|feature_ids| {

            // Always allocate bitset
            let mut bits = AnnotationBits::new(num_features);

            // If empty -> return zero bitset
            if feature_ids.is_empty() {
                return Ok(bits);
            }

            for fid in feature_ids {
                bits.set(fid)?;  // will validate range
            }

            Ok(bits)
        })
        .collect()
}

