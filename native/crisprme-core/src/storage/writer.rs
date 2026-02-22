use std::fs::{File, OpenOptions};
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;

use crate::alignment::alignment::Alignment;

#[derive(Clone, Debug, Default)]
pub struct AlignmentBatchDescr {
    pub alignments_count: usize,
    pub batcher_id: Option<usize>,
}

/// Multithread alignment writer
#[derive(Clone, Debug)]
pub struct AlignmentBatchWriter {
    /// File handler
    file: Arc<File>,
    /// Global write offset
    offset: Arc<AtomicUsize>,
}

impl AlignmentBatchWriter {
    /// Create a new writer
    pub fn open<P: AsRef<Path>>(path: P) -> Self {
        let file = OpenOptions::new()
            .truncate(true)
            .create(true)
            .write(true)
            .open(path)
            .expect("unable to open write file");

        Self {
            offset: Arc::new(AtomicUsize::new(0)),
            file: Arc::new(file),
        }
    }

    /// Write a slice of alignments
    pub fn write_from_memory(&self, alignments: &[Alignment]) {
        let bytes = std::mem::size_of_val(alignments);
        let offset = self
            .offset
            .fetch_add(bytes, std::sync::atomic::Ordering::SeqCst);

        // SAFETY: Alignment is repr(C)
        // Reinterpret the memory of the batch as bytes
        let memory = unsafe { std::slice::from_raw_parts(alignments.as_ptr() as *const u8, bytes) };

        let mut byte_written = 0;
        while byte_written < bytes {
            let written = self
                .file
                .write_at(&memory[byte_written..], offset as u64)
                .expect("unable to write alignments");

            byte_written += written;
        }
    }
}
