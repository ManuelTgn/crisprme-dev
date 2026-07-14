//! Typed column access over untyped `DynColumn` storage.
//!
//! A [`Shape`] ZST encodes the memory layout so that [`Column`] can provide
//! safe, typed iterators and accessors without per-element dispatch.
//!
//! Two shapes are provided:
//! - [`Scalar`] -- one `T` per row, layout `(rows,)`
//! - [`Array<N>`] -- one `[T; N]` per row, layout `(rows, N)`, row-major
//!
//! [`ColumnGroup`] provides N independent scalar sub-columns backed by N
//! separate [`DynColumn`] slots.

use crate::{
    frame::DynColumn,
    memory::{ChunkArray, MemoryPool},
    shared::Share,
};
use bytemuck::Pod;
use std::marker::PhantomData;

/// Typed, borrowed view over a [`DynColumn`].
///
/// Created via [`ColumnScalar::new`] from a `&mut DynColumn` reference.
pub struct Column<'frame, T: Pod> {
    inner: &'frame mut DynColumn,
    _type: PhantomData<T>,
}

impl<'frame, T: Pod> Column<'frame, T> {
    /// Create a new typed column from a DynColumn reference
    pub fn new(inner: &'frame mut DynColumn) -> Self {
        Self {
            _type: PhantomData,
            inner,
        }
    }

    /// Allocate this column type from the pool
    pub fn alloc(&mut self, pool: &MemoryPool, rows: usize) {
        *self.inner =
            DynColumn::new(ChunkArray::alloc::<T>(pool, rows).expect("unable to allocate column"));
    }

    /// Allocate this column with the same amount of rows of another
    pub fn alloc_like<K: Pod>(&mut self, pool: &MemoryPool, other: &Column<'_, K>) {
        self.alloc(pool, other.rows())
    }

    /// Number of rows in this column
    pub fn rows(&self) -> usize {
        self.inner.as_ref().len_of::<T>()
    }

    /// Size in bytes of the rows
    pub const fn row_bytes(&self) -> usize {
        std::mem::size_of::<T>()
    }

    /// Returns the index of an element in memory
    pub fn index(&self, row: usize) -> (usize, usize) {
        self.inner.as_ref().map_index::<T>(row)
    }

    /// Share: freeze `other` and point `self` at the same Arc'd data
    pub fn shared(&mut self, other: &mut Column<'_, T>) {
        *self.inner = other.inner.share();
    }

    /// Take: move ownership from `other` into `self`, leaving `other` empty
    pub fn taken(&mut self, other: &mut Column<'_, T>) {
        *self.inner = other.inner.take();
    }

    /// Reference to row `row`
    pub fn get(&self, row: usize) -> &T {
        self.inner.as_ref().get::<T>(row)
    }

    /// Mutable reference to row `row`
    pub fn get_mut(&mut self, row: usize) -> &mut T {
        self.inner.as_mut().get_mut::<T>(row)
    }

    /// Reference to row at `chunk` and `offset`
    pub fn get_fast(&self, chunk: usize, offset: usize) -> &T {
        self.inner.as_ref().get_fast::<T>(chunk, offset)
    }

    /// Mutable reference to row at `chunk` and `offset`
    pub fn get_fast_mut(&mut self, chunk: usize, offset: usize) -> &mut T {
        self.inner.as_mut().get_fast_mut::<T>(chunk, offset)
    }

    /// Per-chunk slices (aligned to `CHUNK_SIZE`)
    pub fn chunks(&self) -> impl Iterator<Item = &[T]> {
        self.inner.as_ref().chunks::<T>()
    }

    /// Per-chunk mutable slices
    pub fn chunks_mut(&mut self) -> impl Iterator<Item = &mut [T]> {
        self.inner.as_mut().chunks_mut::<T>()
    }

    pub fn slice(&mut self, idx: usize, len: usize) -> &[T] {
        self.inner.as_ref().subchunk::<T>(idx, len)
    }

    pub fn slice_mut(&mut self, idx: usize, len: usize) -> &mut [T] {
        self.inner.as_mut().subchunk_mut::<T>(idx, len)
    }

    /// Flat iterator over all rows
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.chunks().flat_map(|c| c.iter())
    }

    /// Flat mutable iterator over all rows
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.chunks_mut().flat_map(|c| c.iter_mut())
    }
}

/// N independent scalar sub-columns, each backed by its own [`DynColumn`] slot.
///
/// Created by the derive macro for `[T; N]` fields annotated with `#[columnar(group)]`.
/// Access individual sub-columns via [`col`](ColumnGroup::col).
pub struct ColumnGroup<'frame, T: Pod, const N: usize> {
    columns: [&'frame mut DynColumn; N],
    _type: PhantomData<T>,
}

impl<'frame, T: Pod, const N: usize> ColumnGroup<'frame, T, N> {
    /// Create a group from N mutable slot references
    pub fn new(columns: [&'frame mut DynColumn; N]) -> Self {
        Self {
            _type: PhantomData,
            columns,
        }
    }

    /// Number of rows (same across all sub-columns; reads from the first)
    pub fn rows(&self) -> usize {
        self.columns[0].as_ref().len_of::<T>()
    }

    /// Size in bytes of the rows
    pub const fn row_bytes(&self) -> usize {
        std::mem::size_of::<T>()
    }

    /// Share: freeze `other` and point `self` at the same Arc'd data
    pub fn shared(&mut self, other: &mut ColumnGroup<'_, T, N>) {
        for i in 0..N {
            *self.columns[i] = other.columns[i].share();
        }
    }

    /// Take: move ownership from `other` into `self`, leaving `other` empty
    pub fn taken(&mut self, other: &mut ColumnGroup<'_, T, N>) {
        for i in 0..N {
            *self.columns[i] = other.columns[i].take();
        }
    }

    /// Access sub-column `idx`
    pub fn col(&mut self, idx: usize) -> Column<'_, T> {
        debug_assert!(idx < N, "idx >= {}, in group", N);
        Column::new(self.columns[idx])
    }

    /// Split the group into all of its columns
    pub fn split(self) -> [Column<'frame, T>; N] {
        self.columns.map(|c| Column::new(c))
    }

    /// Allocate all N sub-columns with the given row count
    pub fn alloc(&mut self, pool: &MemoryPool, rows: usize) {
        for slot in self.columns.iter_mut() {
            **slot = DynColumn::new(
                ChunkArray::alloc::<T>(pool, rows).expect("unable to allocate group sub-column"),
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod test {
    use super::*;
    use crate::memory::{CHUNK_SIZE, MemoryPool};

    fn make_pool() -> MemoryPool {
        MemoryPool::new(CHUNK_SIZE * 8, |_, _| {})
    }

    // -- Scalar tests --

    #[test]
    fn scalar_alloc_and_rows() {
        let pool = make_pool();
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<u32>::new(&mut dyn_col);
        col.alloc(&pool, 100);
        assert_eq!(col.rows(), 100);
    }

    #[test]
    fn scalar_write_then_read() {
        let pool = make_pool();
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<u32>::new(&mut dyn_col);
        col.alloc(&pool, 64);

        for (i, v) in col.iter_mut().enumerate() {
            *v = i as u32;
        }

        for (i, v) in col.iter().enumerate() {
            assert_eq!(*v, i as u32);
        }
    }

    #[test]
    fn scalar_get_and_get_mut() {
        let pool = make_pool();
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<u32>::new(&mut dyn_col);
        col.alloc(&pool, 16);

        *col.get_mut(5) = 42;
        assert_eq!(*col.get(5), 42);
    }

    #[test]
    fn scalar_chunks_yield_all_elements() {
        let pool = make_pool();
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<u32>::new(&mut dyn_col);
        col.alloc(&pool, 200);

        let total: usize = col.chunks().map(|c| c.len()).sum();
        assert_eq!(total, 200);
    }

    #[test]
    fn scalar_cross_chunk_boundary() {
        let pool = make_pool();
        let elems_per_chunk = CHUNK_SIZE / std::mem::size_of::<u32>();
        // Allocate across 2 chunks
        let count = elems_per_chunk + 100;
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<u32>::new(&mut dyn_col);
        col.alloc(&pool, count);

        for (i, v) in col.iter_mut().enumerate() {
            *v = i as u32;
        }

        // Verify continuity across chunk boundary
        for (i, v) in col.iter().enumerate() {
            assert_eq!(*v, i as u32, "mismatch at index {i}");
        }
        assert_eq!(col.rows(), count);
    }

    // -- Array tests --

    #[test]
    fn array_alloc_and_rows() {
        let pool = make_pool();
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<[u32; 4]>::new(&mut dyn_col);
        col.alloc(&pool, 50);
        assert_eq!(col.rows(), 50);
    }

    #[test]
    fn array_elements_are_contiguous() {
        let pool = make_pool();
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<[u32; 4]>::new(&mut dyn_col);
        col.alloc(&pool, 32);

        for (i, arr) in col.iter_mut().enumerate() {
            *arr = [
                i as u32 * 4,
                i as u32 * 4 + 1,
                i as u32 * 4 + 2,
                i as u32 * 4 + 3,
            ];
        }

        for (i, arr) in col.iter().enumerate() {
            assert_eq!(arr[0], i as u32 * 4);
            assert_eq!(arr[1], i as u32 * 4 + 1);
            assert_eq!(arr[2], i as u32 * 4 + 2);
            assert_eq!(arr[3], i as u32 * 4 + 3);
        }
    }

    #[test]
    fn array_get_and_get_mut() {
        let pool = make_pool();
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<[u32; 4]>::new(&mut dyn_col);
        col.alloc(&pool, 10);

        *col.get_mut(7) = [10, 20, 30, 40];
        assert_eq!(*col.get(7), [10, 20, 30, 40]);
    }

    #[test]
    fn array_cross_chunk_boundary() {
        let pool = make_pool();
        // [u32; 4] is 16 bytes, chunk fits 65536/16 = 4096 elements
        let count = 4096 + 100;
        let mut dyn_col = DynColumn::Empty;
        let mut col = Column::<[u32; 4]>::new(&mut dyn_col);
        col.alloc(&pool, count);

        for (i, arr) in col.iter_mut().enumerate() {
            *arr = [i as u32; 4];
        }

        for (i, arr) in col.iter().enumerate() {
            assert_eq!(*arr, [i as u32; 4], "mismatch at row {i}");
        }
        assert_eq!(col.rows(), count);
    }

    // -- ColumnGroup tests --

    #[test]
    fn group_alloc_and_rows() {
        let pool = make_pool();
        let mut slots = [DynColumn::Empty, DynColumn::Empty, DynColumn::Empty];
        let [s0, s1, s2] = &mut slots;
        let mut grp = ColumnGroup::<u32, 3>::new([s0, s1, s2]);
        grp.alloc(&pool, 100);
        assert_eq!(grp.rows(), 100);
    }

    #[test]
    fn group_columns_are_independent() {
        let pool = make_pool();
        let mut slots = [DynColumn::Empty, DynColumn::Empty, DynColumn::Empty];
        let [s0, s1, s2] = &mut slots;
        let mut grp = ColumnGroup::<u32, 3>::new([s0, s1, s2]);
        grp.alloc(&pool, 16);

        // Write different values to each sub-column
        for (i, v) in grp.col(0).iter_mut().enumerate() {
            *v = i as u32;
        }
        for (i, v) in grp.col(1).iter_mut().enumerate() {
            *v = (i as u32) * 10;
        }
        for (i, v) in grp.col(2).iter_mut().enumerate() {
            *v = (i as u32) * 100;
        }

        // Read back and verify independence
        let c0: Vec<u32> = grp.col(0).iter().copied().collect();
        let c1: Vec<u32> = grp.col(1).iter().copied().collect();
        let c2: Vec<u32> = grp.col(2).iter().copied().collect();

        assert_eq!(c0, (0..16u32).collect::<Vec<_>>());
        assert_eq!(c1, (0..16u32).map(|i| i * 10).collect::<Vec<_>>());
        assert_eq!(c2, (0..16u32).map(|i| i * 100).collect::<Vec<_>>());
    }

    #[test]
    fn group_sub_column_correct_length() {
        let pool = make_pool();
        let mut slots = [
            DynColumn::Empty,
            DynColumn::Empty,
            DynColumn::Empty,
            DynColumn::Empty,
        ];
        let [s0, s1, s2, s3] = &mut slots;
        let mut grp = ColumnGroup::<u32, 4>::new([s0, s1, s2, s3]);
        grp.alloc(&pool, 50);

        for idx in 0..4 {
            let count: usize = grp.col(idx).chunks().map(|c| c.len()).sum();
            assert_eq!(count, 50, "sub-column {idx} should have 50 elements");
        }
    }

    // -- shared_from / taken_from tests --

    #[test]
    fn shared_from_freezes_source() {
        let pool = make_pool();
        let mut src_dyn = DynColumn::Empty;
        let mut dst_dyn = DynColumn::Empty;

        {
            let mut src = Column::<u32>::new(&mut src_dyn);
            src.alloc(&pool, 16);
            for (i, v) in src.iter_mut().enumerate() {
                *v = i as u32;
            }
        }

        {
            let mut src = Column::<u32>::new(&mut src_dyn);
            let mut dst = Column::<u32>::new(&mut dst_dyn);
            dst.shared(&mut src);
        }

        assert!(src_dyn.is_frozen());
        assert!(dst_dyn.is_frozen());

        // Both see the same data
        let src_col = Column::<u32>::new(&mut src_dyn);
        let src_vals: Vec<u32> = src_col.iter().copied().collect();

        let dst_col = Column::<u32>::new(&mut dst_dyn);
        let dst_vals: Vec<u32> = dst_col.iter().copied().collect();

        assert_eq!(src_vals, dst_vals);
    }

    #[test]
    fn taken_from_empties_source() {
        let pool = make_pool();
        let mut src_dyn = DynColumn::Empty;
        let mut dst_dyn = DynColumn::Empty;

        {
            let mut src = Column::<u32>::new(&mut src_dyn);
            src.alloc(&pool, 8);
            for (i, v) in src.iter_mut().enumerate() {
                *v = i as u32;
            }
        }

        {
            let mut src = Column::<u32>::new(&mut src_dyn);
            let mut dst = Column::<u32>::new(&mut dst_dyn);
            dst.taken(&mut src);
        }

        assert!(src_dyn.is_empty());
        assert!(dst_dyn.is_owned());

        let dst_col = Column::<u32>::new(&mut dst_dyn);
        let vals: Vec<u32> = dst_col.iter().copied().collect();
        assert_eq!(vals, (0..8u32).collect::<Vec<_>>());
    }

    // -- zip across same-type columns --

    #[test]
    fn zip_same_type_columns() {
        let pool = make_pool();
        let mut a_dyn = DynColumn::Empty;
        let mut b_dyn = DynColumn::Empty;

        {
            let mut a = Column::<u32>::new(&mut a_dyn);
            a.alloc(&pool, 64);
            for (i, v) in a.iter_mut().enumerate() {
                *v = i as u32;
            }
        }

        {
            let mut b = Column::<u32>::new(&mut b_dyn);
            b.alloc(&pool, 64);
        }

        // Chunk-level zip: a.chunks() and b.chunks_mut() are aligned
        let a = Column::<u32>::new(&mut a_dyn);
        let mut b = Column::<u32>::new(&mut b_dyn);

        for (a_chunk, b_chunk) in a.chunks().zip(b.chunks_mut()) {
            for (a_val, b_val) in a_chunk.iter().zip(b_chunk.iter_mut()) {
                *b_val = *a_val * 2;
            }
        }

        let b_ref = Column::<u32>::new(&mut b_dyn);
        for (i, v) in b_ref.iter().enumerate() {
            assert_eq!(*v, i as u32 * 2);
        }
    }
}
