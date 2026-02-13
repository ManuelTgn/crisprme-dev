use std::io::Read;
use std::os::unix::fs::FileExt;
use std::path::PathBuf;
use std::sync::Arc;
use std::fs::File;

use crate::sequence::iupac::Iupac;
use crate::memory::arena::{ArenaBox, Memory};
use crate::memory::batch::SequenceRingBatch;

/// Iterate over all ids arrays present in the binary file
pub struct BinaryPositionReader<'mem, R: Read> {
    mem: &'mem Memory<'mem>,
    reader: R,
}

impl<'mem, R: Read> BinaryPositionReader<'mem, R> {
    pub fn new(mem: &'mem Memory, reader: R) -> Self {
        Self { mem, reader }
    }
}

impl<'mem, R: Read> Iterator for BinaryPositionReader<'mem, R> {
    type Item = ArenaBox<'mem, [u32]>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut len_buf = [0u8; 4];
        if let Err(e) = self.reader.read_exact(&mut len_buf) {
            return if e.kind() == std::io::ErrorKind::UnexpectedEof {
                None
            } else {
                panic!("unable to read array len");
            };
        }

        let len = u32::from_le_bytes(len_buf) as usize;
        let mut result = self.mem.alloc_slice_fill(len, 0);
        for slot in &mut *result {
            let mut buffer = [0; 4];
            self.reader
                .read_exact(&mut buffer)
                .expect("unable to read array bytes");

            *slot = u32::from_le_bytes(buffer);
        }

        Some(result)
    }
}

/// Iterate over all sequences present in the binary file
pub struct BinarySequenceReader<'mem, R: Read> {
    buffer: ArenaBox<'mem, [u8]>,
    sequence_len: usize,
    reader: R,
}

impl<'mem, R: Read> BinarySequenceReader<'mem, R> {
    pub fn new(mem: &'mem Memory, reader: R, sequence_len: usize) -> Self {
        let buffer = mem.alloc_slice_fill(sequence_len, 0u8);
        Self {
            buffer,
            sequence_len,
            reader,
        }
    }
}

impl<'mem, R: Read> Iterator for BinarySequenceReader<'mem, R> {
    type Item = &'mem [Iupac];
    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_exact(&mut self.buffer) {
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => None,
            Err(_) => panic!("unable to read sequence"),
            Ok(()) => Some(unsafe {
                std::slice::from_raw_parts(self.buffer.as_ptr() as *const Iupac, self.sequence_len)
            }),
        }
    }
}

/// Describes a `SequenceBatch` content
#[derive(Debug, Clone, Default)]
pub struct SequenceBatchDescr {
    pub sequence_count: usize,
    pub sequence_len: usize,
    pub global_offset: usize,
}

/// Iterate over batches of sequences in a file
#[derive(Debug)]
pub struct BinarySequenceBatchReader {
    sequence_len: usize,
    sequence_count: usize,
    batch_size: usize,
    file: Arc<File>,
    cuda: bool,
}

impl BinarySequenceBatchReader {
    pub fn open(input: &PathBuf, sequence_len: usize, batch_size: usize) -> Self {
        let file = File::open(input).expect("unable to open file");
        let file_len = file.metadata().unwrap().len();
        let sequence_count = file_len as usize / sequence_len;
        Self {
            file: Arc::new(file),
            sequence_count,
            sequence_len,
            batch_size,
            cuda: false,
        }
    }

    pub fn new(file: File, sequence_len: usize, batch_size: usize, cuda: bool) -> Self {
        let file_len = file.metadata().unwrap().len();
        let sequence_count = file_len / sequence_len as u64;
        Self {
            file: Arc::new(file),
            sequence_count: sequence_count as usize,
            sequence_len,
            batch_size,
            cuda,
        }
    }

    /// Return standard batch size
    pub fn batch_size(&self) -> usize {
        self.batch_size
    }

    /// Return an iterator over all batches indices and sizes
    pub fn batches<'s>(&'s self) -> impl Iterator<Item = (usize, usize)> + 's {
        let num_batches = self.sequence_count.div_ceil(self.batch_size);
        (0..num_batches).map(move |idx| {
            let running_elem_count = (idx + 1) * self.batch_size;
            let size = if running_elem_count > self.sequence_count {
                self.sequence_count % self.batch_size
            } else {
                self.batch_size
            };
            (idx, size)
        })
    }

    fn read_at_iupac(&self, buffer: &mut [Iupac], offset: usize) {
        // SAFETY: Iupac is repr(C)
        // Reinterpret the input buffer as a slice of bytes
        let buffer =
            unsafe { std::slice::from_raw_parts_mut(buffer.as_mut_ptr() as *mut u8, buffer.len()) };

        let mut bytes_read = 0;
        while bytes_read < buffer.len() {
            let read = self
                .file
                .read_at(&mut buffer[bytes_read..], (offset + bytes_read) as u64)
                .expect("unable to read bytes");

            bytes_read += read;
        }
    }

    pub fn describe(&self, idx: usize, real_count: usize) -> SequenceBatchDescr {
        SequenceBatchDescr {
            global_offset: idx * self.batch_size,
            sequence_count: real_count,
            sequence_len: self.sequence_len,
        }
    }

    pub fn read_batch(&self, idx: usize, real_count: usize, batch: &mut SequenceRingBatch) {
        assert!(batch.iupac().len() >= real_count * self.sequence_len);
        assert!(batch.ids().len() >= real_count);

        // Calculate offset of the inner sequence indices relative to the entire file
        let global_offset = idx * self.batch_size;
        self.read_at_iupac(
            &mut batch.iupac_mut()[..real_count * self.sequence_len],
            global_offset * self.sequence_len,
        );

        // Generate sequence ids
        (0..real_count).for_each(|i| {
            batch.ids_mut()[i] = (global_offset + i) as u32;
        });
    }
}
