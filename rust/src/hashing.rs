use std::collections::HashMap;
use ahash::AHashMap;
use pyo3::prelude::*;
use std::fs::OpenOptions;
use std::io::{Write, BufWriter, Result as IOResult, Seek};
use serde::Serialize;
use rmp_serde;  // handle MessagePack

use crate::target::Target;  // raw output from scanner

// Type alias for the complex value data, just for cleaner code
type OccurrenceData = Vec<(String, usize, bool)>;

// Define the fixed-width record for the Index File
#[derive(Debug, Serialize)]  // Serialize for writing to the index file
struct TargetIndexRecord {
    // Offset in the sequences file (Targets.bin)
    seq_offset: u64,
    // Length of the sequence (k-mer)
    seq_len: u8,
    // Offset in the occurrences file (Occurrences.bin)
    data_offset: u64,
    // Number of occurrence records (count of (contig, pos, strand) tuples)
    data_count: u32,
}

/// Represents a collection of target sites grouped by their unique sequence.
/// 
/// The key is the target sequence (IUPAC bitmasks), and the value is a vector
/// of all occurrences (positions and orientations) of that specific sequence.
/// This structure efficiently collapses redundant target sequences found during the scan.
/// 
/// This struct is exposed to Python via PyO3.
#[pyclass]
#[derive(Debug, Clone, Serialize)]
pub struct HashedTargets {
    /// The map where keys are the unique target bitmasks (`Vec<u8>`) and values are
    /// vectors of occurrence data (`(contig, position, orientation)`).
    /// This directly replaces the Python `Dict[bytes, Target]` structure
    #[pyo3(get)]
    pub targets: HashMap<Vec<u8>, OccurrenceData>,
}

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
    let mut targets_map: AHashMap<Vec<u8>, OccurrenceData> = AHashMap::with_capacity(raw_targets.len());

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

    HashedTargets { targets: targets_map.into_iter().collect() }
}

impl HashedTargets {
    /// Saves the target data into a three-file system, optimized for maximum write speed.
    ///
    /// This implementation uses zero-copy serialization (`bincode::serialize_into`) 
    /// for both binary files, eliminating unnecessary heap allocations.
    ///
    /// # Files Generated:
    /// 1. Targets.bin: Contiguous target sequences (raw binary bytes).
    /// 2. Occurrences.bin: Densely packed serialized occurrence records (bincode).
    /// 3. Index.bin: Fixed-width records pointing to the data (bincode).
    ///
    /// # Arguments
    /// * `path_prefix` - The base path for the files (e.g., "/data/scan_results").
    pub fn save_indexed_binary(&self, path_prefix: &str, is_first_chunk: bool) -> IOResult<()> {
        let seq_path = format!("{}_targets.bin", path_prefix);
        let data_path = format!("{}_occurrences.bin", path_prefix);
        let index_path = format!("{}_index.bin", path_prefix);

        // --- Determine Open options (Truncate vs. Append) ---
        // Open or create the files, and set the cursor to the end for appending
        let mut binding = OpenOptions::new();
        let open_options = binding.create(true).write(true); 
        
        if is_first_chunk { 
            // Truncate mode: Clear existing files and start fresh
            open_options.truncate(true);
        } else {
            // Append mode: Keep existing content and continue writing from the end
            open_options.append(true);
        };

        // --- Writers for the three files ---
        let seq_file = open_options.open(seq_path)?;
        let data_file = open_options.open(data_path)?;
        let index_file = open_options.open(index_path)?;
        
        let mut seq_writer = BufWriter::new(seq_file);
        let mut data_writer = BufWriter::new(data_file);
        let mut index_writer = BufWriter::new(index_file);

        // Initialize Offsets based on existing file size
        // Read the current position (end of file) to determine the starting offset for this chunk.
        let mut current_seq_offset: u64 = seq_writer.stream_position()?;
        let mut current_data_offset: u64 = data_writer.stream_position()?;
        // Index offset is implicitly tracked by consecutive fixed-width writes, 
        // but reading the initial position is good practice.
        let mut _current_index_offset: u64 = index_writer.stream_position()?; 

        // Iterate over the HashMap to populate all three files simultaneously
        for (sequence, occurrences) in self.targets.iter() {
            let data_count = occurrences.len() as u32;  // count target occurrences

            // 1. Write Occurrence Data (Occurrences.bin) - ZERO COPY
            // The occurrence data offset is the current end of file position.
            let data_offset = current_data_offset;

            // Calculate size before writing to update the offset for the next record
            let encoded_size = bincode::serialized_size(occurrences)
                .map_err(|e| {
                    eprintln!("Error serializing occurrence data for sequence: {:?}", sequence);
                    std::io::Error::new(std::io::ErrorKind::Other, format!("Bincode serialization failed: {}", e))
                })?;
            
            // Write directly into the buffer (no intermediate Vec<u8> allocation) (Occurrences.bin)
            bincode::serialize_into(&mut data_writer, occurrences)
                .map_err(|e| {
                    eprintln!("Error serializing occurrence data for sequence {:?}", sequence);
                    std::io::Error::new(std::io::ErrorKind::Other, format!("Bincode serialization failed: {}", e))
                })?;

            // Update data offset for the next record
            current_data_offset += encoded_size;

            // 2. Write Target Sequence (Targets.bin) - Zero-copy write
            seq_writer.write_all(sequence.as_slice())?;

            let total_seq_bytes_written = sequence.len() as u64 + 1; // +1 for the newline
            let seq_len = sequence.len() as u8; // Original sequence length (assuming target < 255)

            // 3. Write Index Record (Index.bin) - Zero-copy
            let record = TargetIndexRecord {
                seq_offset: current_seq_offset,
                seq_len: seq_len, 
                data_offset: data_offset,
                data_count: data_count,
            };

            bincode::serialize_into(&mut index_writer, &record)
                .map_err(|e| {
                    eprintln!("Error writing index record for sequence: {:?}", sequence);
                    std::io::Error::new(std::io::ErrorKind::Other, format!("Bincode index write failed: {}", e))
                })?;

            // Update offsets for the next record
            current_seq_offset += total_seq_bytes_written; // CRITICAL: Update by total bytes written
        }

        // --- CRITICAL FLUSHING STEPS ---
        // Explicitly drop writers to force flushing and closing files
        std::mem::drop(seq_writer);
        std::mem::drop(data_writer);
        std::mem::drop(index_writer);
        
        // If the code reaches here, all writers have been successfully flushed and closed.
        Ok(())
    }
}