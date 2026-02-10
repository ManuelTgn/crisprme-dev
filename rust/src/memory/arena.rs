use bump_scope::{Bump, BumpBox, BumpScope};

/// Arena allocator based on bump_scope
pub struct Arena(Bump);

impl Arena {
    pub fn alloc(size: usize) -> Self {
        Self(Bump::with_size(size))
    }
}

/// Allocation on the arena
pub type ArenaBox<'arena, T> = BumpBox<'arena, T>;

/// Transmute a box from a type to another
pub fn transmute_box<'arena, A, B>(src: ArenaBox<'arena, A>) -> ArenaBox<'arena, B> {
    unsafe {
        let ptr = src.into_raw();
        ArenaBox::from_raw(ptr.cast::<B>())
    }
}

/// Scope in the arena
pub type Memory<'arena> = BumpScope<'arena>;

impl std::ops::Deref for Arena {
    type Target = Bump;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Arena {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
