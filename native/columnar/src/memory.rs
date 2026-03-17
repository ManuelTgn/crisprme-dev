//! Chunk-based memory pool with typed access via `bytemuck`.
//!
//! Memory is pre-allocated as fixed-size [`CHUNK_SIZE`] byte blocks managed by
//! a lock-free [`MemoryPool`]. Individual chunks are returned automatically on
//! drop. [`ChunkArray`] groups one or more chunks into a contiguous logical
//! array with typed, per-chunk and region-based iterators.

use bytemuck::Pod;
use crossbeam::queue::ArrayQueue;
use tracing::{Level, error, span};
use std::{mem::ManuallyDrop, sync::{Arc, Weak}, time::Instant};

/// Size of a single memory chunk (64 KB)
pub const CHUNK_SIZE: usize = 65536;

/// Maximum alowed seconds before `try_acquire_spinning` fails
pub const QUEUE_ACQUIRE_MAX_SECS: f32 = 1.0;

#[derive(Debug)]
pub enum MemoryError {
    MaximumTriesReached,
    OutOfMemory
}

/// Raw backing storage for one chunk, cache-line aligned.
#[repr(align(64))]
#[derive(Debug)]
pub struct ChunkBytes([u8; CHUNK_SIZE]);
impl ChunkBytes {
    pub fn boxed() -> Box<Self> {
        Box::new(Self([0u8; CHUNK_SIZE]))
    }
}

/// A single chunk handle. Returns its memory to the pool on drop.
pub struct Chunk {
    memory: ManuallyDrop<Box<ChunkBytes>>,
    queue: Weak<MemoryPoolQueue>,
}

impl Drop for Chunk {
    fn drop(&mut self) {
        let memory = unsafe { ManuallyDrop::take(&mut self.memory) };
        if let Some(queue) = self.queue.upgrade() {
            let _ = queue.push(memory);
        }
    }
}

impl Chunk {

    /// How many elements of type T can a chunk hold
    pub const fn capacity_of<T: Pod>() -> usize {
        CHUNK_SIZE / std::mem::size_of::<T>()
    }

    /// Typed read-only view of the chunk bytes
    pub fn as_slice<T: Pod>(&self) -> &[T] {
        bytemuck::cast_slice(&self.memory.0)
    }

    /// Typed mutable view of the chunk bytes
    pub fn as_slice_mut<T: Pod>(&mut self) -> &mut [T] {
        bytemuck::cast_slice_mut(&mut self.memory.0)
    }
}

/// A logical array spanning one or more [`Chunk`]s.
///
/// Provides typed element access, per-chunk iteration, and sub-range
/// region iterators. The last chunk may be partially used.
#[derive(Default)]
pub struct ChunkArray {
    chunks: Vec<Chunk>,
    len: usize
}

impl ChunkArray {

    /// Allocate enough chunks from `pool` to hold `count` of type `T`
    pub fn alloc<T: Pod>(pool: &MemoryPool, count: usize) -> Result<Self, MemoryError> {
        let used_bytes = count.checked_mul(std::mem::size_of::<T>())
            .expect("overflow");

        let chunks_count = used_bytes.div_ceil(CHUNK_SIZE);
        let mut chunks = Vec::with_capacity(chunks_count);
        for _ in 0..chunks_count {
            // Get chunk without considering backpressure
            chunks.push(pool.acquire()
                .ok_or(MemoryError::OutOfMemory)?);
        }

        tracing::debug!("ChunkArray with {} chunks ({} bytes not used)",
            chunks_count, chunks_count * CHUNK_SIZE - used_bytes);

        Ok(Self { chunks, len: used_bytes })
    }

    pub fn len(&self) -> usize { self.len }
    pub fn chunk_count(&self) -> usize { self.chunks.len() }
    pub fn capacity(&self) -> usize { self.chunk_count() * CHUNK_SIZE }

    /// Total number of actual `T` elements across all chunks
    pub fn len_of<T: Pod>(&self) -> usize {
        self.len() / std::mem::size_of::<T>()
    }

    /// Total number of storable 'T' elements across all chunks
    pub fn capacity_of<T: Pod>(&self) -> usize {
        self.capacity() / std::mem::size_of::<T>()
    }

    /// Map an index for an element of type 'T' to a `(chunk_index, offset_within_chunk)` pair
    pub const fn map_index<T: Pod>(&self, idx: usize) -> (usize, usize) {
        let capacity = Chunk::capacity_of::<T>();
        (idx / capacity, idx % capacity)
    }

    /// Reference to element `idx` of type 'T' (computes chunk + offset)
    pub fn get<T: Pod>(&self, idx: usize) -> &T {
        let (cidx, offset) = self.map_index::<T>(idx);
        self.get_fast(cidx, offset)
    }

    /// Mutable reference to element `idx` of type 'T'
    pub fn get_mut<T: Pod>(&mut self, idx: usize) -> &mut T {
        let (cidx, offset) = self.map_index::<T>(idx);
        self.get_fast_mut(cidx, offset)
    }

    /// Reference of type 'T' by pre-computed `(chunk_index, offset)`, skips index math
    pub fn get_fast<T: Pod>(&self, cidx: usize, offset: usize) -> &T {
        &self.chunks[cidx].as_slice()[offset]
    }

    /// Mutable reference of type 'T' by pre-computed `(chunk_index, offset)`
    pub fn get_fast_mut<T: Pod>(&mut self, cidx: usize, offset: usize) -> &mut T {
        &mut self.chunks[cidx].as_slice_mut()[offset]
    }

    /// Per-chunk slices of `T` (last slice may be shorter than a full chunk)
    pub fn chunks<T: Pod>(&self) -> impl Iterator<Item = &[T]> {

        let total_elements = self.len_of::<T>();
        let chunk_capacity = Chunk::capacity_of::<T>();

        self.chunks.iter().enumerate().map(move |(i, chunk)| {

            let start = i * chunk_capacity;
            let remaining = total_elements.saturating_sub(start);
            let len = remaining.min(chunk_capacity);

            &chunk.as_slice::<T>()[..len]
        })
    }

    /// Per-chunk mutable slices of `T`
    pub fn chunks_mut<T: Pod>(&mut self) -> impl Iterator<Item = &mut [T]> {

        let total_elements = self.len_of::<T>();
        let chunk_capacity = Chunk::capacity_of::<T>();

        self.chunks.iter_mut().enumerate().map(move |(i, chunk)| {

            let start = i * chunk_capacity;
            let remaining = total_elements.saturating_sub(start);
            let len = remaining.min(chunk_capacity);

            &mut chunk.as_slice_mut::<T>()[..len]
        })
    }

    /// Slice into an index range within a single chunk
    pub fn subchunk<T: Pod>(&self, idx: usize, len: usize) -> &[T] {
        let (chunk, offset) = self.map_index::<T>(idx);
        &self.chunks[chunk]
            .as_slice::<T>()[offset .. offset + len]
    }

    /// Mutable slice into an index range within a single chunk
    pub fn subchunk_mut<T: Pod>(&mut self, idx: usize, len: usize) -> &mut [T] {
        let (chunk, offset) = self.map_index::<T>(idx);
        &mut self.chunks[chunk]
            .as_slice_mut::<T>()[offset .. offset + len]
    }
}

/// Minimum elements-per-chunk across a set of element sizes
pub fn region_min_size(element_sizes: &[usize]) -> usize {
    element_sizes.iter()
        .map(|s| CHUNK_SIZE / *s)
        .min().unwrap()
}

type MemoryPoolQueue = ArrayQueue<Box<ChunkBytes>>;

/// Lock-free pool of pre-allocated [`Chunk`]s.
///
/// Chunks are acquired with [`acquire`](Self::acquire) or
/// [`try_acquire_spinning`](Self::try_acquire_spinning) and returned
/// automatically when dropped.
#[derive(Debug, Clone)]
pub struct MemoryPool {
    chunks: Arc<MemoryPoolQueue>,
}

impl MemoryPool {

    /// Allocate `bytes / CHUNK_SIZE` chunks, calling `visitor` for each one
    pub fn new<V>(bytes: usize, visitor: V) -> Self
    where
        V: Fn(*const u8, usize)
    {
        let chunks_count = bytes.div_ceil(CHUNK_SIZE);
        let chunks = ArrayQueue::new(chunks_count);
        for _ in 0..chunks_count {
            let boxed = ChunkBytes::boxed();
            visitor(boxed.as_ref().0.as_ptr(), CHUNK_SIZE);
            chunks.push(boxed)
                .expect("Unable to create memory pool");
        }

        Self { chunks: Arc::new(chunks) }
    }

    /// Try to pop a chunk from the pool, returns `None` if empty
    pub fn acquire(&self) -> Option<Chunk> {
        self.chunks.pop().map(|memory| Chunk {
            memory: ManuallyDrop::new(memory),
            queue: Arc::downgrade(&self.chunks),
        })
    }

    /// Spin-wait for a chunk
    pub fn acquire_spinning(&self) -> Chunk {
        loop {
            if let Some(chunk) = self.acquire() {
                return chunk;
            }
            // No reason for busy loop, the pipeline has backpressure
            std::thread::yield_now();
        }
    }

    /// Spin-wait for a chunk, fail after [`QUEUE_ACQUIRE_MAX_SECS`] seconds
    pub fn try_acquire_spinning(&self) -> Result<Chunk, MemoryError> {

        // Busy loop for some tries
        for _ in 0..100 {
            if let Some(chunk) = self.acquire() {
                return Ok(chunk);
            }
        }

        let instant = Instant::now();
        loop {
            // No reason for busy loop, the pipeline has backpressure
            std::thread::yield_now();

            // We have a chunk
            if let Some(chunk) = self.acquire() {
                return Ok(chunk);
            }

            // Max seconds reached, probably a deadlock or inadequate pipeline resources
            if instant.elapsed().as_secs_f32() > QUEUE_ACQUIRE_MAX_SECS {
                error!("memory allocation failed: MaximumTriesReached");
                return Err(MemoryError::MaximumTriesReached);
            }
        }
    }

    /// Total number of chunk slots in the pool
    pub fn capacity(&self) -> usize {
        return self.chunks.capacity()
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn iteration_mut() {
        let pool = MemoryPool::new(1_000_000, |_, _| { });
        let mut span = ChunkArray::alloc::<[u32; 4]>(&pool, 256)
            .unwrap();

        for chunk in span.chunks_mut::<[u32; 4]>() {
            for v in chunk {
                v[0] = 0; v[1] = 1; v[2] = 2; v[3] = 3;
            }
        }
    }

    #[test]
    fn chunk_returns_to_pool_on_drop() {
        let pool = MemoryPool::new(CHUNK_SIZE, |_, _| { });
        assert_eq!(pool.capacity(), 1);

        let chunk = pool.acquire().unwrap();
        assert!(pool.acquire().is_none());

        drop(chunk);
        assert!(pool.acquire().is_some());
    }

    #[test]
    fn chunk_span_correct_chunk_count() {
        let pool = MemoryPool::new(CHUNK_SIZE * 4, |_, _| { });

        // Exactly one chunk worth of u32s
        let elems_per_chunk = CHUNK_SIZE / std::mem::size_of::<u32>();
        let span = ChunkArray::alloc::<u32>(&pool, elems_per_chunk).unwrap();
        assert_eq!(span.chunk_count(), 1);
        assert_eq!(span.len(), CHUNK_SIZE);
        drop(span);

        // One element over -> needs two chunks
        let span = ChunkArray::alloc::<u32>(&pool, elems_per_chunk + 1).unwrap();
        assert_eq!(span.chunk_count(), 2);
        drop(span);

        // Zero elements -> zero chunks
        let span = ChunkArray::alloc::<u32>(&pool, 0).unwrap();
        assert_eq!(span.chunk_count(), 0);
        assert_eq!(span.len(), 0);
    }

    #[test]
    fn write_then_read_back() {
        let pool = MemoryPool::new(CHUNK_SIZE * 4, |_, _| { });
        let mut span = ChunkArray::alloc::<u64>(&pool, 1024).unwrap();

        // Write sequential values
        let mut val = 0u64;
        for chunk in span.chunks_mut::<u64>() {
            for elem in chunk {
                *elem = val;
                val += 1;
            }
        }

        // Read them back
        let mut expected = 0u64;
        for chunk in span.chunks::<u64>() {
            for &elem in chunk {
                assert_eq!(elem, expected);
                expected += 1;
            }
        }
        assert_eq!(expected, 1024);
    }

    #[test]
    fn pool_exhaustion() {
        let pool = MemoryPool::new(CHUNK_SIZE * 2, |_, _| { });

        let _c1 = pool.acquire().unwrap();
        let _c2 = pool.acquire().unwrap();
        assert!(pool.acquire().is_none());

        let result = pool.try_acquire_spinning();
        assert!(matches!(result, Err(MemoryError::MaximumTriesReached)));
    }

    #[test]
    fn visitor_called_for_each_chunk() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let count = AtomicUsize::new(0);
        let total_bytes = AtomicUsize::new(0);

        let _pool = MemoryPool::new(CHUNK_SIZE * 3, |_ptr, size| {
            count.fetch_add(1, Ordering::Relaxed);
            total_bytes.fetch_add(size, Ordering::Relaxed);
        });

        assert_eq!(count.load(Ordering::Relaxed), 3);
        assert_eq!(total_bytes.load(Ordering::Relaxed), CHUNK_SIZE * 3);
    }

    #[test]
    fn chunk_span_capacity_vs_used() {
        let pool = MemoryPool::new(CHUNK_SIZE * 4, |_, _| { });

        // Allocate slightly more than one chunk
        let elems = CHUNK_SIZE / std::mem::size_of::<u32>() + 1;
        let span = ChunkArray::alloc::<u32>(&pool, elems).unwrap();

        assert_eq!(span.len(), elems * std::mem::size_of::<u32>());
        assert_eq!(span.capacity(), CHUNK_SIZE * 2);
        assert!(span.capacity() >= span.len());
    }

    #[test]
    fn last_chunk_slice_is_truncated() {
        let pool = MemoryPool::new(CHUNK_SIZE * 4, |_, _| { });
        let elems_per_chunk = CHUNK_SIZE / std::mem::size_of::<u32>();

        // 1.5 chunks worth of elements
        let total = elems_per_chunk + elems_per_chunk / 2;
        let span = ChunkArray::alloc::<u32>(&pool, total).unwrap();

        let slices: Vec<_> = span.chunks::<u32>().collect();
        assert_eq!(slices.len(), 2);
        assert_eq!(slices[0].len(), elems_per_chunk);
        assert_eq!(slices[1].len(), elems_per_chunk / 2);
    }
}