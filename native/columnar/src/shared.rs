//! Copy-on-write cell with three states: `Owned`, `Frozen` (shared via `Arc`), or `Empty`.
//!
//! [`ShareCell`] is the building block for zero-copy column sharing.
//! Freezing an owned cell wraps it in an `Arc`; subsequent [`Share::share`]
//! calls clone the `Arc` without copying the data.

use std::sync::Arc;

/// A cell that is either exclusively owned, shared-immutable, or vacant.
///
/// - `Owned(T)` -- mutable access allowed
/// - `Frozen(Arc<T>)` -- read-only, cheaply cloneable
/// - `Empty` -- no data
pub enum ShareCell<T> {
    Frozen(Arc<T>),
    Owned(T),
    Empty
}

impl<T> ShareCell<T> {

    /// Create a new owned cell
    pub fn new(value: T) -> Self {
        Self::Owned(value)
    }

    pub fn is_owned(&self)  -> bool { matches!(self, ShareCell::Owned(_))  }
    pub fn is_frozen(&self) -> bool { matches!(self, ShareCell::Frozen(_)) }
    pub fn is_empty(&self)  -> bool { matches!(self, ShareCell::Empty)     }

    /// Move the inner value out, leaving `self` empty
    pub fn take(&mut self) -> ShareCell<T> {
        let old = std::mem::replace(self, ShareCell::Empty);
        return old;
    }

    /// Transition to `Frozen`: wraps an owned value in `Arc`, no-ops if already frozen
    pub fn freeze(&mut self) {
        let old = std::mem::replace(self, ShareCell::Empty);
        *self = match old {
            ShareCell::Frozen(arc) => ShareCell::Frozen(arc),
            ShareCell::Owned(v) => ShareCell::Frozen(Arc::new(v)),
            ShareCell::Empty => ShareCell::Empty,
        };
    }

    /// Drop the inner value, transition to `Empty`
    pub fn clear(&mut self) {
        *self = ShareCell::Empty
    }

    /// Try to get a reference (`Owned` or `Frozen`), `None` if empty
    pub fn get_ref(&self) -> Option<&T> {
        match self {
            ShareCell::Frozen(arc) => Some(arc),
            ShareCell::Owned(v) => Some(v),
            _ => None,
        }
    }

    /// Try to get a mutable reference (`Owned` only), `None` otherwise
    pub fn get_mut(&mut self) -> Option<&mut T> {
        match self {
            ShareCell::Owned(v) => Some(v),
            _ => None,
        }
    }

    /// Reference to inner data. Panics if empty.
    pub fn as_ref(&self) -> &T {
        self.get_ref()
            .expect("cell should not be empty")
    }

    /// Mutable reference to inner data. Panics if not owned.
    pub fn as_mut(&mut self) -> &mut T {
        self.get_mut()
            .expect("cell should not empty or frozen")
    }

}

/// Produce a shared copy of `self`, freezing the source in the process.
pub trait Share {
    fn share(&mut self) -> Self;
}

impl<T> Share for ShareCell<T> {
    /// Freeze, then clone the `Arc` -- both source and result are `Frozen`.
    fn share(&mut self) -> Self {
        self.freeze();
        match self {
            ShareCell::Frozen(arc) => ShareCell::Frozen(arc.clone()),
            ShareCell::Owned(_) => unreachable!(),
            ShareCell::Empty => ShareCell::Empty,
        }
    }
}

/// A default ShareCell is empty
impl<T> Default for ShareCell<T> {
    fn default() -> Self {
        Self::Empty
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn share_owned_freezes_source() {
        let mut cell = ShareCell::new(4);
        assert!(cell.is_owned());
        let shared = cell.share();
        assert!(shared.is_frozen());
        assert!(cell.is_frozen());
    }

    #[test]
    fn share_empty_stays_empty() {
        let mut cell = ShareCell::<u32>::Empty;
        let shared = cell.share();
        assert!(cell.is_empty());
        assert!(shared.is_empty());
    }
}