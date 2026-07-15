use columnar::pipeline::{PipelineError, Sink};
use std::marker::PhantomData;

pub mod writer;

/// A Sink that does nothing, it just consumes output
pub struct NullSink<T>(PhantomData<T>);

impl<T> NullSink<T> {
    pub fn new() -> Self {
        Self(PhantomData)
    }
}

impl<T> Sink for NullSink<T>
where
    T: Send + 'static,
{
    type I = T;

    fn name() -> &'static str {
        "NullSink"
    }
    fn consume(&mut self, _item: Self::I) -> Result<(), PipelineError> {
        Ok(())
    }
}
