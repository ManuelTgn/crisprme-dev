use std::collections::HashMap;
use pyo3::prelude::*;

use crate::target::Target;  // raw output from scanner

/// Represents a collection of target sites grouped by their unique sequence.
/// 
/// The key is the target sequence (IUPAC bitmasks), and the value is a vector
/// of all occurrences (positions and orientations) of that specific sequence.
/// This structure efficiently collapses redundant target sequences found during the scan.
/// 
/// This struct is exposed to Python via PyO3.
#[pyclass]
#[derive(Debug, Clone)]
pub struct HashedTargets {
    /// The map where keys are the unique target bitmasks (`Vec<u8>`) and values are
    /// vectors of occurrence data (`(contig, position, orientation)`).
    /// This directly replaces the Python `Dict[bytes, Target]` structure
    #[pyo3(get)]
    pub targets: HashMap<Vec<u8>, Vec<(String, usize, bool)>>,
}

#[pymethods]
impl HashedTargets {
    /// internal constructor to initialize the HashedTargets struct
    pub fn new() -> Self {
        HashedTargets {
            targets: HashMap::new(),
        }
    }
}

/// Performs the core logic of grouping raw scan results by sequence.
/// 
/// This function iterates over the `raw_targets` (the complete list of found sites)
/// and consolidates all occurrences that share the exact same sequence bitmask 
/// into a single entry in a HashMap. This replaces the slow Python dictionary creation loop.
///
/// # Arguments
/// * `raw_targets` - A vector of all individual `Target` matches found by the scanner.
/// 
/// # Returns
/// A fully constructed `HashedTargets` object containing the collapsed, unique targets
pub fn hash_and_group_targets(raw_targets: Vec<Target>) -> HashedTargets {
    let mut targets_map: HashMap<Vec<u8>, Vec<(String, usize, bool)>> = HashMap::new();

    for target in raw_targets {
        // use HashMap's entry API for efficient lookup and insertion
        let entry = targets_map
            // the unique key is the target sequence (Vec<u8>)
            .entry(target.target)
            // if the key is new, insert an empty vector
            .or_insert_with(Vec::new);

        // append the target data (contig, position, strand) to the correct sequence group
        entry.push((target.contig, target.position, target.orientation));
    }

    HashedTargets { targets: targets_map}
}


