use bump_scope::{Bump, BumpScope};

pub struct Arena(Bump);
impl Arena {

    /// Create new arena with `bytes` capacity
    pub fn with_capacity(bytes: usize) -> Self {
        Self(Bump::with_size(bytes))
    }

    /// Create a temporary memory region to work with
    pub fn scoped<R>(&mut self, f: impl FnOnce(&mut Memory) -> R) -> R {
        self.0.scoped(f)
    }
}

pub type Memory<'arena> = BumpScope<'arena>;