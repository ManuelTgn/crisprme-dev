//! Schema-tagged wrapper around [`DynFrame`].
//!
//! [`Columnar<S>`] pairs a [`DynFrame`] with a compile-time [`Schema`] marker
//! so that a derive macro can generate typed column accessors.

use std::marker::PhantomData;

use crate::{frame::DynFrame, memory::MemoryPool, shared::Share};

/// Marker trait for types whose fields map to [`DynFrame`] slots.
/// Intended to be implemented by a derive macro.
pub unsafe trait Schema {
    /// Total number of [`DynColumn`](crate::frame::DynColumn) slots this schema requires.
    const SLOTS: usize;
}

/// A [`DynFrame`] tagged with a compile-time schema `S`.
pub struct TypedFrame<S: Schema> {
    pub _schema: PhantomData<S>,
    pub frame: DynFrame,
}

impl<S: Schema> TypedFrame<S> {
    /// Create an empty frame with this schema
    pub fn empty() -> Self {
        Self {
            _schema: PhantomData,
            frame: DynFrame::empty(S::SLOTS),
        }
    }

    /// Mutable access to the inner frame (used by derive macro)
    pub fn frame_mut(&mut self) -> &mut DynFrame {
        &mut self.frame
    }

    /// Attach a schema to an existing frame.
    /// # Safety
    /// The caller must ensure the frame's layout matches `S`.
    pub unsafe fn attach(frame: DynFrame) -> Self {
        Self {
            _schema: PhantomData,
            frame,
        }
    }

    /// Consume the wrapper and return the inner frame
    pub fn detach(self) -> DynFrame {
        self.frame
    }
}

impl<S: Schema> Share for TypedFrame<S> {
    fn share(&mut self) -> Self {
        unsafe { Self::attach(self.frame.share()) }
    }
}

/// Access columns from multiple frames without nested closures.
///
/// ```ignore
/// cols!(mut s = sequences, mut p = positions, mut m = merged => {
///     m.seq_id.shared_from(&mut s.id);
/// });
/// ```
#[macro_export]
macro_rules! cols {
    // Base case: single binding
    (mut $name:ident = $frame:expr => $body:expr) => {
        $frame.with_cols(|mut $name| $body)
    };
    ($name:ident = $frame:expr => $body:expr) => {
        $frame.with_cols(|$name| $body)
    };
    // Recursive case: peel off first binding, nest the rest
    (mut $name:ident = $frame:expr, $($rest:tt)+) => {
        $frame.with_cols(|mut $name| cols!($($rest)+))
    };
    ($name:ident = $frame:expr, $($rest:tt)+) => {
        $frame.with_cols(|$name| cols!($($rest)+))
    };
}
