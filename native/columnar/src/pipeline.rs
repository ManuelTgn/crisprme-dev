
//! Multi-threaded data processing pipeline.
//!
//! A [`Pipeline`] is a chain of [`Stage`]s connected by bounded crossbeam channels.
//! Each stage runs in one or more worker threads. Data flows left-to-right; backpressure
//! is applied automatically when a downstream channel is full.
//!
//! Dropping a [`PipelineHandle`] closes the input and joins all worker threads.

use std::thread::JoinHandle;
use crossbeam::channel::{Receiver, Sender, TrySendError, bounded};
use metrics::{counter, histogram};
use crate::MemoryPool;

/// Error returned when a stage fails or a channel is disconnected.
#[derive(Debug)]
pub struct StageError;

/// Sink for items produced by a [`Stage`]. Wraps crossbeam send errors.
pub trait Emit<T> {
    fn emit(&self, item: T) -> Result<(), StageError>;
}

impl<T> Emit<T> for Sender<T> {

    fn emit(&self, item: T) -> Result<(), StageError> {
        if let Err(e) = self.try_send(item) {
            match e {
                TrySendError::Disconnected(_) => return Err(StageError),
                TrySendError::Full(item) => {
                    counter!("pipeline.backpressure").increment(1);
                    self.send(item)
                        .map_err(|_| StageError)?;
                },
            }
        }
        Ok(())
    }
}

// =============================================================================
// Stage trait
// =============================================================================

/// A processing stage that transforms items of type `I` into items of type `O`.
///
/// Implement [`process`](Stage::process) to define the transformation.
/// [`run`](Stage::run) is the thread entry-point and drives the recv/process loop.
pub trait Stage: Send + 'static {

    type I: Send + 'static;
    type O: Send + 'static;

    // Get name of the stage
    fn name() -> &'static str;

    /// Transform one input item, emitting zero or more outputs via `emitter`.
    fn process(&mut self, input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError>;

    /// Thread entry-point: receives items and calls [`process`](Stage::process) until the channel closes.
    fn run(&mut self, src: Receiver<Self::I>, dst: impl Emit<Self::O>) -> Result<(), StageError> {
        while let Ok(item) = src.recv() {
            counter!("pipeline.items", "stage" => Self::name()).increment(1);

            let t = std::time::Instant::now();
            self.process(item, &dst)?;
            histogram!("pipeline.elapsed_ns", "stage" => Self::name())
                .record(t.elapsed().as_nanos() as f64);
        }
        self.shutdown()
    }

    /// Called once after the input channel closes. Override to flush or release resources.
    fn shutdown(&mut self) -> Result<(), StageError> { Ok(()) }

}

/// A stage that generates items autonomously — no external input required.
///
/// Used as the entry point of a [`pipeline_with_source`] pipeline.
pub trait Source: Send + 'static {

    type O: Send + 'static;

    // Get name of the stage
    fn name() -> &'static str;

    /// Produce the next item, or `None` to signal exhaustion.
    fn next(&mut self) -> Result<Option<Self::O>, StageError>;

    /// Thread entry-point: calls [`next`](Source::next) until it returns `None` or an error.
    fn run(&mut self, dst: impl Emit<Self::O>) -> Result<(), StageError> {
        while let Ok(Some(item)) = self.next() {
            counter!("pipeline.items", "stage" => Self::name()).increment(1);

            let t = std::time::Instant::now();
            dst.emit(item)?;
            histogram!("pipeline.elapsed_ns", "stage" => Self::name())
                .record(t.elapsed().as_nanos() as f64);
        }
        self.shutdown()
    }

    /// Called once after the last item has been emitted.
    fn shutdown(&mut self) -> Result<(), StageError> { Ok(()) }

}

/// A terminal stage that consumes items without producing output.
///
/// Used as the endpoint of a [`Pipeline::sink`] pipeline.
pub trait Sink: Send + 'static {

    type I: Send + 'static;

    // Get name of the stage
    fn name() -> &'static str;

    /// Thread entry-point: receives items and calls [`consume`](Sink::consume) until the channel closes.
    fn run(&mut self, src: Receiver<Self::I>) -> Result<(), StageError> {
        while let Ok(item) = src.recv() {
            counter!("pipeline.items", "stage" => Self::name()).increment(1);

            let t = std::time::Instant::now();
            self.consume(item)?;
            histogram!("pipeline.elapsed_ns", "stage" => Self::name())
                .record(t.elapsed().as_nanos() as f64);
            
        }
        self.shutdown()
    }

    /// Process one item from the pipeline's output.
    fn consume(&mut self, item: Self::I) -> Result<(), StageError>;

    /// Called once after the input channel closes.
    fn shutdown(&mut self) -> Result<(), StageError> { Ok(()) }
}

// =============================================================================
// Pipeline
// =============================================================================

// A dynamic recursive function to construct the pipeline
type PipelineFn<T> = Box<dyn FnOnce(&MemoryPool, usize) -> (Vec<JoinHandle<()>>, Receiver<T>)>;

/// A linear pipeline builder.
///
/// `T` is the output type of the last node added.
pub struct Pipeline<T: Send + 'static> {
    f: PipelineFn<T>,
}

impl Pipeline<()> {

    /// Create a driven pipeline entry-point
    pub fn driven<T: Send + 'static>(input_cap: usize) -> (Driven<T>, Pipeline<T>) {
        
        let (tx, rx) = bounded::<T>(input_cap);
        let driven = Driven(Some(tx));

        let pipeline = Pipeline {
            f: Box::new(move |_pool, _cap| {
                let handles = Vec::new();
                (handles, rx)
            }),
        };

        (driven, pipeline)
    }

    /// Create a source for the pipeline
    pub fn source<F, S>(workers: usize, f: F) -> Pipeline<S::O>
    where
        S: Source,
        F: Fn(&MemoryPool, usize) -> S + Send + 'static,
    {
        Pipeline {
            f: Box::new(move |pool, cap| {
                
                // Create channel connection
                let (tx, rx) = bounded::<S::O>(cap);

                // Spawn workers
                let mut handles = Vec::new();
                for i in 0..workers {
                    let mut source = f(pool, i);
                    
                    let dst = tx.clone();
                    handles.push(std::thread::spawn(move || {
                        source.run(dst).unwrap();
                    }));
                }

                (handles, rx)
            }),
        }
    }
}

impl<T: Send + 'static> Pipeline<T> {

    /// Create a transform stage in the pipeline
    pub fn stage<F, S>(self, workers: usize, f: F) -> Pipeline<S::O>
    where
        S: Stage<I = T>,
        F: Fn(&MemoryPool, usize) -> S + Send + 'static,
    {
        let prev = self.f;
        Pipeline {
            f: Box::new(move |pool, cap| {
                
                // Call previous pipeline stages
                let (mut handles, src_rx) = prev(pool, cap);

                // Create channel connection
                let (tx, rx) = bounded::<S::O>(cap);

                // Spawn workers
                for i in 0..workers {
                    let mut stage = f(pool, i);

                    let src = src_rx.clone();
                    let dst = tx.clone();
                    handles.push(std::thread::spawn(move || {
                        stage.run(src, dst).unwrap();
                    }));
                }

                (handles, rx)
            }),
        }
    }

    /// Create a single-worker transform stage, consuming the factory `FnOnce`.
    /// Use this when the stage value cannot be cloned (e.g. a Python object).
    pub fn stage_once<F, S>(self, f: F) -> Pipeline<S::O>
    where
        S: Stage<I = T>,
        F: FnOnce(&MemoryPool) -> S + Send + 'static,
    {
        let prev = self.f;
        Pipeline {
            f: Box::new(move |pool, cap| {

                // Call previous pipeline stages
                let (mut handles, src_rx) = prev(pool, cap);

                // Create channel connection
                let (tx, rx) = bounded::<S::O>(cap);

                let mut stage = f(pool);
                handles.push(std::thread::spawn(move || {
                    stage.run(src_rx, tx).unwrap();
                }));

                (handles, rx)
            }),
        }
    }

    /// Create a sink for the pipeline
    pub fn sink<F, S>(self, workers: usize, f: F) -> PipelineExecutable
    where
        S: Sink<I = T>,
        F: Fn(&MemoryPool, usize) -> S + Send + 'static,
    {
        let prev = self.f;
        PipelineExecutable {
            launcher: Box::new(move |pool, cap| {
                
                // Call previous pipeline stages
                let (mut handles, src_rx) = prev(pool, cap);

                // Spawn workers
                for i in 0..workers {
                    let mut sink = f(pool, i);

                    let src = src_rx.clone();
                    handles.push(std::thread::spawn(move || {
                        sink.run(src).unwrap();
                    }));
                }

                handles
            }),
        }
    }
}

/// A fully built pipeline ready to execute.
pub struct PipelineExecutable {
    launcher: Box<dyn FnOnce(&MemoryPool, usize) -> Vec<JoinHandle<()>>>,
}

impl PipelineExecutable {
    pub fn execute(self, pool: &MemoryPool, cap: usize) -> PipelineHandle {
        let handles = (self.launcher)(pool, cap);
        PipelineHandle {
            handles
        }
    }
}

/// Handle to the pipeline
pub struct PipelineHandle {
    handles: Vec<JoinHandle<()>>
}

impl PipelineHandle {
    /// Wait for all threads to finish
    pub fn join(&mut self) {
        for h in self.handles.drain(..) {
            h.join().unwrap();
        }
    }
}

impl Drop for PipelineHandle {
    fn drop(&mut self) {
        self.join();
    }
}

// =============================================================================
// Utilities
// =============================================================================

/// Sender for a driven source
pub struct Driven<T>(Option<Sender<T>>);
impl<T: Send> Driven<T> {

    pub fn send(&self, item: T) -> Result<(), StageError> {
        self.0.as_ref()
            .expect("pipeline was closed")
            .send(item)
            .map_err(|_e| StageError)
    }

    // Close the pipeline
    pub fn close(&mut self) {
        let _ = self.0.take();
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod test {
    use crate::memory::CHUNK_SIZE;
    use super::*;

    struct SourceStage(u32);
    impl Source for SourceStage {
        type O = u32;

        fn name() -> &'static str { "SourceStage" }

        fn next(&mut self) -> Result<Option<Self::O>, StageError> {
            if self.0 == 0 { return Ok(None); }
            else { self.0 -= 1; }
            Ok(Some(self.0))
        }
    }

    struct SinkStage;
    impl Sink for SinkStage {
        type I = u32;

        fn name() -> &'static str { "SinkStage" }

        fn consume(&mut self, _item: Self::I) -> Result<(), StageError> {
            Ok(())
        }
    }

    struct DoubleStage;
    impl Stage for DoubleStage {
        type I = u32;
        type O = u32;

        fn name() -> &'static str { "DoubleStage" }

        fn process(&mut self, input: Self::I, emitter: &impl Emit<Self::O>) -> Result<(), StageError> {
            emitter.emit(input * 2)
        }
    }

    #[test]
    fn it_works() {
        let pool = MemoryPool::new(CHUNK_SIZE, |_, _| { });

        let pipeline = Pipeline::source(2, |_pool, _worker| SourceStage(10))
            .stage(2, |_pool, _worker| DoubleStage)
            .stage(3, |_pool, _worker| DoubleStage)
            .sink(1, |_pool, _worker| SinkStage);

        let _handle = pipeline
            .execute(&pool, 2);
    }

    #[test]
    fn it_works_driven() {
        let pool = MemoryPool::new(CHUNK_SIZE, |_, _| { });

        let (input, pipeline) = Pipeline::driven::<u32>(2);
        let pipeline = pipeline
            .stage(2, |_pool, _worker| DoubleStage)
            .stage(3, |_pool, _worker| DoubleStage)
            .sink(1, |_pool, _worker| SinkStage);

        let _handle = pipeline
            .execute(&pool, 2);

        for i in 0..10 {
            input.send(i)
                .unwrap();
        }

        drop(input); // close the channel so workers can finish
        // _handle drops here, joining all threads
    }
}