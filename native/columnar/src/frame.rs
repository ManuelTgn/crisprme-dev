//! Untyped, fixed-width frame of up to [`FRAME_MAX_SLOTS`] columns.
//!
//! A [`DynFrame`] holds [`DynColumn`] slots (each a `ShareCell<ChunkArray>`)
//! and a row count. Typed access is layered on top via `shape::Column`.

use crate::{
    memory::ChunkArray,
    shared::{Share, ShareCell},
};

/// An untyped column: a shareable chunk array
pub type DynColumn = ShareCell<ChunkArray>;

/// Fixed-width collection of [`DynColumn`] slots.
pub struct DynFrame {
    slots: Vec<DynColumn>,
}

impl DynFrame {
    /// Create a frame with all slots empty and zero rows
    pub fn empty(slots: usize) -> Self {
        let mut v = Vec::with_capacity(slots);
        v.resize_with(slots, DynColumn::default);
        Self { slots: v }
    }

    /// Returns a raw pointer to the slot storage.
    ///
    /// # Safety
    /// The caller must ensure exclusive access to each slot index and that
    /// `idx < self.slots.len()` when dereferencing.
    pub unsafe fn slots_ptr(&mut self) -> *mut DynColumn {
        self.slots.as_mut_ptr()
    }
}

impl Share for DynFrame {
    /// Freeze every slot and return a shared copy of the frame
    fn share(&mut self) -> Self {
        let mut result = DynFrame::empty(self.slots.len());
        for (src, dst) in self.slots.iter_mut().zip(result.slots.iter_mut()) {
            *dst = src.share();
        }
        result
    }
}
