use crossbeam::channel::SendError;
use std::ops::{Deref, DerefMut};
use crossbeam_channel::TryRecvError;
use tracing::{error, trace, warn};

// use crate::bindings;

use super::GpuPtr;

/// A slot of the ring buffer, containing owned CPU and (optional) GPU data.
///
/// Each `RingSlot` owns both its CPU memory (as a `Vec<T>`) and an optional GPU buffer
/// (allocated with [`GpuPtr`]). It also tracks whether the CPU or GPU memory
/// has been modified since the last synchronization.
///
/// Typically, slots are managed via [`RingSlotLease`], which ensures they are
/// automatically returned to the ring when dropped.
pub struct RingSlot<T> {
    /// Buffer stored on the GPU (if allocated).
    storage_gpu: Option<GpuPtr<T>>,
    /// Buffer stored in CPU memory.
    storage_cpu: Vec<T>,
    /// True if CPU memory has been modified since last sync to GPU
    cpu_dirty: bool,
    /// True if GPU memory has been modified since last sync to CPU
    gpu_dirty: bool,
}

impl<T> RingSlot<T> {
    /// Returns the number of elements that can fit in this slot.
    ///
    /// Note: this is the capacity of the underlying `Vec`, not its current length.
    pub fn capacity(&self) -> usize {
        self.storage_cpu.capacity()
    }

    /// Returns a raw pointer to the GPU buffer.
    ///
    /// # Panics
    /// Panics if no GPU memory has been allocated for this slot.
    pub fn gpu_ptr(&self) -> *const T {
        self.storage_gpu.as_ref().unwrap().as_ptr() as *const T
    }

    /// Returns a mutable raw pointer to the GPU buffer, marking it as dirty.
    ///
    /// # Panics
    /// Panics if no GPU memory has been allocated for this slot.
    pub fn gpu_ptr_mut(&mut self) -> *mut T {
        self.gpu_dirty = true;
        self.storage_gpu.as_mut().unwrap().as_ptr()
    }

    /// Copy CPU buffer contents into GPU memory.
    ///
    /// This will only perform the copy if [`cpu_dirty`] is set.
    /// After synchronization, the CPU buffer is considered clean.
    #[tracing::instrument(name = "ring_slot", skip_all)]
    pub fn sync_cpu_to_gpu(&mut self, bytes: Option<usize>) {
        if self.cpu_dirty {
            if let Some(gpu) = &self.storage_gpu {
                trace!("buffer with dirty CPU memory, synchronizing");
                // unsafe{
                //     bindings::cuda::memcpy_to_gpu::<T>(
                //         self.storage_cpu.as_ptr(),
                //         gpu.as_ptr(),
                //         bytes.unwrap_or(self.capacity()),
                //     );
                // }
                self.cpu_dirty = false;
            }
        }
    }

    /// Copy GPU buffer contents into CPU memory.
    ///
    /// This will only perform the copy if [`gpu_dirty`] is set.
    /// After synchronization, the GPU buffer is considered clean.
    #[tracing::instrument(name = "ring_slot", skip_all)]
    pub fn sync_gpu_to_cpu(&mut self, bytes: Option<usize>) {
        if self.gpu_dirty {
            if let Some(gpu) = &self.storage_gpu {
                trace!("buffer with dirty GPU memory, synchronizing");
                // unsafe {
                //     bindings::cuda::memcpy_to_cpu::<T>(
                //         self.storage_cpu.as_mut_ptr(),
                //         gpu.as_ptr() as *const T,
                //         bytes.unwrap_or(self.capacity()),
                //     );
                // }
                self.gpu_dirty = false;
            }
        }
    }
}

impl<T> Deref for RingSlot<T> {
    type Target = [T];
    fn deref(&self) -> &Self::Target {
        if self.gpu_dirty && self.storage_gpu.is_some() {
            warn!("access to RingBuffer with dirty GPU memory");
        }
        self.storage_cpu.as_slice()
    }
}

impl<T> DerefMut for RingSlot<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        if self.gpu_dirty && self.storage_gpu.is_some() {
            warn!("access to RingBuffer with dirty GPU memory");
        }
        self.cpu_dirty = true;
        self.storage_cpu.as_mut_slice()
    }
}

/// RAII guard for a ring buffer slot.
///
/// A `RingSlotLease` ensures that a buffer slot is **always returned to the ring**
/// when dropped, even if the user forgets to explicitly release it.
///
/// - If a slot is dropped without being marked as "used", a warning is printed.
/// - If the ring buffer has been destroyed, the slot is dropped silently with a warning.
///
/// Typically, users interact with a slot through a custom [`RingAdapter`].
pub struct RingSlotLease {
    /// Sender used to return the slot to the producer pool.
    drop: crossbeam::channel::Sender<RingSlot<u8>>,
    /// The owned slot, wrapped in an `Option` so that it can be `take()`n during commit.
    slot: Option<RingSlot<u8>>,
    /// Whether this buffer was ever "used" before being returned.
    used: bool,
}

impl Deref for RingSlotLease {
    type Target = RingSlot<u8>;
    fn deref(&self) -> &Self::Target {
        self.slot.as_ref().unwrap()
    }
}

impl DerefMut for RingSlotLease {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.slot.as_mut().unwrap()
    }
}

impl Drop for RingSlotLease {
    fn drop(&mut self) {
        if let Some(slot) = self.slot.take() {
            if !self.used {
                warn!("possibly unused buffer sent back to ring");
            }
            match self.drop.send(slot) {
                Ok(()) => {}
                Err(SendError(_)) => {
                    error!("unable to drop buffer to ring, ring buffer has been destroyed");
                }
            }
        }
    }
}

/// The producer side of the ring buffer.
///
/// A `Producer` owns a pool of available slots and is responsible for:
/// - **Acquiring** slots for writing data.
/// - **Committing** slots once they are ready to be consumed.
///
/// Producers can be cloned to allow multiple producer threads.
pub struct Producer<A: RingAdapter> {
    producer_rx: crossbeam::channel::Receiver<RingSlot<u8>>,
    producer_tx: crossbeam::channel::Sender<RingSlot<u8>>,
    consumer_tx: crossbeam::channel::Sender<(A::Descr, RingSlot<u8>)>,
}

impl<A: RingAdapter> Producer<A> {
    /// Acquire a buffer slot
    #[tracing::instrument(name = "ring", skip_all)]
    pub fn acquire(&self) -> A {
        trace!("acquiring buffer from ring");

        let slot = self.producer_rx.recv().unwrap();
        A::attach(
            RingSlotLease {
                drop: self.producer_tx.clone(),
                slot: Some(slot),
                used: false,
            },
            A::Descr::default(),
        )
    }

    // Send a buffer to the consumers
    #[tracing::instrument(name = "ring", skip_all)]
    pub fn commit(&self, adapter: A) {
        trace!("committing buffer to ring");

        let (descr, mut lease) = adapter.detach();
        let slot = lease.slot.take().unwrap();
        lease.used = true;
        self.consumer_tx
            .send((descr, slot))
            .expect("unable to send buffer to ring consumer");
    }

    // Send a buffer to the consumers
    #[tracing::instrument(name = "ring", skip_all)]
    pub fn commit_with_descriptor(&self, adapter: A, descr: A::Descr) {
        trace!("committing buffer to ring");

        let (_descr, mut lease) = adapter.detach();
        let slot = lease.slot.take().unwrap();
        lease.used = true;
        self.consumer_tx
            .send((descr, slot))
            .expect("unable to send buffer to ring consumer");
    }

    /// Close the channel, no more items will be produced
    pub fn close(self) {
        drop(self.producer_tx);
        drop(self.producer_rx);
        drop(self.consumer_tx);
    }
}

impl<A: RingAdapter> Clone for Producer<A> {
    fn clone(&self) -> Self {
        Self {
            producer_tx: self.producer_tx.clone(),
            producer_rx: self.producer_rx.clone(),
            consumer_tx: self.consumer_tx.clone(),
        }
    }
}

/// The consumer side of the ring buffer.
///
/// A `Consumer` receives slots that have been committed by a [`Producer`].
/// It is responsible for:
/// - **Receiving** slots ready for processing.
/// - **Finishing** with slots and returning them to the producer pool.
///
/// Consumers can be cloned to allow multiple consumer threads.
pub struct Consumer<A: RingAdapter> {
    producer_tx: crossbeam::channel::Sender<RingSlot<u8>>,
    consumer_rx: crossbeam::channel::Receiver<(A::Descr, RingSlot<u8>)>,
}

impl<A: RingAdapter> Consumer<A> {
    /// Get a buffer and metadata ready to be processed, if the channel is closed returns None
    #[tracing::instrument(name = "ring", skip_all)]
    pub fn recv(&self) -> Option<A> {
        trace!("receiving buffer from ring");
        let result = self.consumer_rx.recv();
        match result {
            // Channel closed
            Err(_) => None,
            Ok((descr, slot)) => Some(A::attach(
                RingSlotLease {
                    drop: self.producer_tx.clone(),
                    slot: Some(slot),
                    used: false,
                },
                descr,
            )),
        }
    }

    /// Get a buffer and metadata ready to be processed
    #[tracing::instrument(name = "ring", skip_all)]
    pub fn try_recv(&self) -> Result<Option<A>, TryRecvError> {
        trace!("receiving buffer from ring");
        let result = self.consumer_rx.try_recv();
        match result {
            // Channel closed
            Err(e) => Err(e),
            Ok((descr, slot)) => Ok(Some(A::attach(
                RingSlotLease {
                    drop: self.producer_tx.clone(),
                    slot: Some(slot),
                    used: false,
                },
                descr,
            ))),
        }
    }

    /// Send a buffer to the producers
    #[tracing::instrument(name = "ring", skip_all)]
    pub fn finish(&self, adapter: A) {
        trace!("finishing buffer to ring");

        let (_descr, mut lease) = adapter.detach();
        let slot = lease.slot.take().unwrap();
        lease.used = true;
        if self.producer_tx.send(slot).is_err() {
            warn!("the producer was closed, no need for this buffer");
        }
    }
}

impl<A: RingAdapter> Clone for Consumer<A> {
    fn clone(&self) -> Self {
        Self {
            producer_tx: self.producer_tx.clone(),
            consumer_rx: self.consumer_rx.clone(),
        }
    }
}

/// Adapter trait for viewing and managing slots in the ring buffer.
///
/// An adapter provides a typed view (`Self`) over a `RingSlotLease` and an
/// associated descriptor (`Descr`). This allows users to attach metadata or
/// structured access patterns without changing the underlying buffer implementation.
///
/// For example, one could implement a `FrameAdapter` that interprets a
/// `RingSlot<u8>` as an RGB image with resolution metadata.
pub trait RingAdapter {
    /// Constructor parameters
    type Descr: Clone + Default;
    /// Create a new adapter on top of a ring buffer slot
    fn attach(slot: RingSlotLease, descr: Self::Descr) -> Self;
    /// Get the slot and descriptor back
    fn detach(self) -> (Self::Descr, RingSlotLease);
    /// Access mutably to the underlying data
    fn as_mut(&mut self) -> &mut RingSlotLease;
    /// Access to the underlying data
    fn as_ref(&self) -> &RingSlotLease;
}

/// Create a new ring buffer with the given number of slots and slot size.
///
/// Each slot contains `slot_bytes` elements of type `u8`.
/// If `gpu` is true, GPU memory is also allocated for each slot.
///
/// Returns a `(Producer<A>, Consumer<A>)` pair, which form the two ends of the ring.
///
/// # Example
/// ```
/// let (producer, consumer) = ring_buffer::<MyAdapter>(4, 1024, true);
/// let mut slot = producer.acquire();
/// slot.as_mut()[0] = 42;
/// producer.commit(slot);
///
/// if let Some(mut slot) = consumer.recv() {
///     assert_eq!(slot.as_ref()[0], 42);
///     consumer.finish(slot);
/// }
/// ```
pub fn ring_buffer<A: RingAdapter>(
    slot_count: usize,
    slot_bytes: usize,
    gpu: bool,
) -> (Producer<A>, Consumer<A>) {
    let (producer_tx, producer_rx) = crossbeam::channel::bounded::<RingSlot<u8>>(slot_count);
    let (consumer_tx, consumer_rx) =
        crossbeam::channel::bounded::<(A::Descr, RingSlot<u8>)>(slot_count);

    // At the beginning all slots are free
    for _ in 0..slot_count {
        // Also allocate on gpu if requested
        let storage_gpu = if gpu {
            Some(GpuPtr::<u8>::alloc(slot_bytes))
        } else {
            None
        };

        producer_tx
            .send(RingSlot {
                storage_cpu: vec![0u8; slot_bytes],
                storage_gpu,
                cpu_dirty: false,
                gpu_dirty: false,
            })
            .expect("unable to initialize ring buffer");
    }

    (
        Producer {
            producer_tx: producer_tx.clone(),
            producer_rx: producer_rx.clone(),
            consumer_tx,
        },
        Consumer {
            producer_tx,
            consumer_rx,
        },
    )
}
