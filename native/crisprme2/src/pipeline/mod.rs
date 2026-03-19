
pub mod sink;
pub mod stage;

#[cfg(test)]
pub mod test {
    use columnar::{MemoryPool, memory::CHUNK_SIZE, pipeline::{Emit, PipelineError}};
    use std::cell::RefCell;

    /// Creates a memory pool sized for use in unit tests.
    pub fn make_pool() -> MemoryPool {
        MemoryPool::new(CHUNK_SIZE * 100, |_, _| {})
    }

    /// Collects items emitted by a stage for inspection in tests.
    pub struct Collector<T>(pub RefCell<Vec<T>>);

    impl<T> Collector<T> {
        pub fn new() -> Self {
            Self(RefCell::new(vec![]))
        }

        pub fn into_inner(self) -> Vec<T> {
            self.0.into_inner()
        }
    }

    impl<T> Emit<T> for Collector<T> {
        fn emit(&self, item: T) -> Result<(), PipelineError> {
            self.0.borrow_mut().push(item);
            Ok(())
        }
    }
}