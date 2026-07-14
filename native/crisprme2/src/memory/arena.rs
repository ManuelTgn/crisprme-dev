//! Bump-allocated arena memory utilities.
//!
//! This module provides a thin wrapper around [`bump_scope::Bump`] to enable
//! fast, short-lived allocations in performance-critical paths
//! (e.g., batch preparation before CUDA alignment).
//!
//! # Why a bump arena?
//!
//! Bump allocation is ideal when:
//! - many small allocations are needed,
//! - all allocations share the same lifetime,
//! - deallocation can happen all at once.
//!
//! In CRISPR alignment workflows, batches of targets, buffers, or
//! intermediate structures often share a common lifetime (e.g., one alignment
//! round). A bump allocator:
//!
//! - avoids per-allocation `malloc/free` overhead,
//! - reduces allocator contention,
//! - improves cache locality,
//! - enables deterministic memory reclamation.
//!
//! Memory is reclaimed wholesale when the arena is dropped or when a scope
//! ends.
//!
//! # Safety model
//!
//! Allocations from a bump arena are tied to the arena lifetime (`'arena`).
//! All objects allocated within an arena must not outlive the arena.
//!
//! The `transmute_box` helper provided here is `unsafe` internally and must
//! only be used when layout compatibility between types is guaranteed.
//!
//! # Typical usage
//!
//! ```ignore
//! let mut arena = Arena::alloc(1 << 20); // 1 MB
//! let scope = arena.scope();
//! let value = scope.alloc(42);
//! ```
//!
//! All memory allocated inside `scope` is invalidated when the scope ends.
//!

use bump_scope::{Bump, BumpBox, BumpScope};

/// Arena allocator backed by [`bump_scope::Bump`].
///
/// This is a thin newtype wrapper that allows you to control how bump
/// allocation is exposed across the codebase while keeping the underlying
/// allocator accessible via `Deref`.
///
/// The arena owns a contiguous memory region from which allocations are
/// served linearly (bump-pointer style).
///
/// # Lifetime
/// Objects allocated from this arena are valid for the lifetime of the arena
/// (or a derived [`BumpScope`]).
///
/// Dropping the arena frees *all* allocations at once.
pub struct Arena(Bump);

impl Arena {
    /// Create a new arena with a preallocated buffer of `size` bytes.
    ///
    /// # Arguments
    /// * `size` — Initial capacity in bytes.
    ///
    /// # Notes
    /// - This does **not** limit total allocation size; the bump allocator
    ///   may grow internally if required (depending on `bump_scope` behavior).
    /// - Choosing a reasonable size upfront reduces reallocations.
    #[inline]
    pub fn alloc(size: usize) -> Self {
        Self(Bump::with_size(size))
    }
}

/// Allocation type within an [`Arena`].
///
/// `ArenaBox<'arena, T>` behaves similarly to `Box<T>` but:
/// - does **not** allocate via the global allocator,
/// - is freed automatically when the arena (or scope) is dropped,
/// - cannot outlive the arena lifetime `'arena`.
///
/// This is ideal for short-lived, high-frequency allocations.
pub type ArenaBox<'arena, T> = BumpBox<'arena, T>;

/// A scoped lifetime view into the arena.
///
/// Scopes allow structured allocation lifetimes inside a longer-lived arena.
/// When a scope ends, all allocations inside that scope are invalidated.
///
/// Use this to isolate temporary allocations within larger workflows.
pub type Memory<'arena> = BumpScope<'arena>;

/// Transmute an arena allocation from type `A` to type `B`.
///
/// This is a zero-cost pointer cast of the underlying allocation.
///
/// # Safety
///
/// This function is **extremely unsafe** unless all of the following hold:
///
/// - `A` and `B` have identical memory layout (`size_of::<A>() == size_of::<B>()`)
/// - Alignment requirements of `B` are compatible with `A`
/// - The underlying bytes are a valid representation of `B`
/// - No drop semantics are violated
///
/// Misuse can lead to:
/// - Undefined behavior
/// - Memory corruption
/// - Violations of Rust aliasing rules
///
/// # When is this acceptable?
///
/// - Reinterpreting raw byte buffers as structured types
/// - Casting between repr(C) or repr(transparent) types
/// - Performance-critical code where layout guarantees are externally enforced
///
/// # Recommendation
///
/// Prefer safe conversions unless this is strictly required for performance
/// or FFI integration.
///
/// # Example
///
/// ```ignore
/// let a: ArenaBox<u32> = scope.alloc(42);
/// let b: ArenaBox<[u8; 4]> = transmute_box(a);
/// ```
#[inline]
pub unsafe fn transmute_box<'arena, A, B>(src: ArenaBox<'arena, A>) -> ArenaBox<'arena, B> {
    // Catches mistakes in debug builds
    debug_assert_eq!(std::mem::size_of::<A>(), std::mem::size_of::<B>());

    let ptr = src.into_raw();
    ArenaBox::from_raw(ptr.cast::<B>())
}

/// Deref to underlying [`Bump`] allocator.
///
/// This allows calling allocator methods directly:
///
/// ```ignore
/// let mut arena = Arena::alloc(1024);
/// let scope = arena.scope();
/// ```
impl std::ops::Deref for Arena {
    type Target = Bump;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Mutable deref to underlying allocator.
///
/// Allows mutable operations such as creating scopes or adjusting allocation state.
impl std::ops::DerefMut for Arena {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
