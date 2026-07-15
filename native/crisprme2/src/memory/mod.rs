use crate::bindings;
use arena::Memory;
use bump_scope::NoDrop;
use std::{
    ops::{Deref, DerefMut},
    ptr::NonNull,
};
use tracing::trace;

pub mod arena;
pub mod batch;
pub mod ring;

pub type CpuBuffer<'s, T> = &'s mut [T];

/// Wrapper for gpu data, this is safe to send between threads
pub struct GpuPtr<T> {
    pub ptr: NonNull<T>,
}

unsafe impl<T: Send> Send for GpuPtr<T> {}
unsafe impl<T: Sync> Sync for GpuPtr<T> {}

impl<T> GpuPtr<T> {
    /// Allocate buffer on the GPU
    pub fn alloc(len: usize) -> Self {
        let ptr = bindings::cuda::malloc::<T>(len);
        trace!("allocated gpu buffer ({len} elements)");
        Self {
            ptr: NonNull::new(ptr).expect("failed CUDA malloc"),
        }
    }

    #[inline]
    pub fn as_ptr(&self) -> *mut T {
        self.ptr.as_ptr()
    }
}

/// Free the memory on the GPU
impl<T> Drop for GpuPtr<T> {
    fn drop(&mut self) {
        unsafe { bindings::cuda::free::<T>(self.ptr.as_ptr()) };
        trace!("dropped gpu buffer");
    }
}

/// Hybrid buffer storing the same data on both CPU and GPU.
///
/// Provides synchronization between CPU and GPU copies and exposes
/// safe access to CPU memory. GPU memory is accessed through raw pointers.
pub struct HybridBuffer<'s, T> {
    /// CPU-side storage (borrowed from bump allocator)
    cpu: CpuBuffer<'s, T>,
    /// GPU-side storage
    gpu: Option<GpuPtr<T>>,
    /// True if CPU memory has been modified since last sync to GPU
    cpu_dirty: bool,
    /// True if GPU memory has been modified since last sync to CPU
    gpu_dirty: bool,
    /// Number of allocated elements presents
    capacity: usize,
}

impl<'s, T: 'static> HybridBuffer<'s, T> {
    /// Creates a new `HybridBuffer` from a CPU slice.
    ///
    /// Initially, the GPU buffer is allocated.
    pub fn from_slice(data: &'s mut [T], gpu: bool) -> Self {
        let gpu_buffer = if gpu {
            Some(GpuPtr::alloc(data.len()))
        } else {
            None
        };

        let mut result = Self {
            capacity: data.len(),
            gpu: gpu_buffer,
            cpu: data,
            cpu_dirty: gpu,
            gpu_dirty: false,
        };

        result.sync_to_gpu();
        result
    }

    /// Returns a raw mutable pointer to the GPU buffer.
    ///
    /// # Safety
    /// The caller must ensure GPU memory is valid before dereferencing.
    pub fn gpu_ptr_mut(&mut self) -> Option<*mut T> {
        match &self.gpu {
            Some(gpu) => {
                self.gpu_dirty = true;
                Some(gpu.as_ptr())
            }
            None => None,
        }
    }

    /// Returns a raw const pointer to the GPU buffer
    ///
    /// # Safety
    /// The caller must ensure GPU memory is valid before dereferencing.
    pub fn gpu_ptr(&self) -> Option<*const T> {
        self.gpu.as_ref().map(|e| e.as_ptr() as *const T)
    }

    /// Synchronizes CPU data to GPU if the CPU is dirty.
    ///
    /// Currently a placeholder; actual GPU memory copy must be implemented.
    pub fn sync_to_gpu(&mut self) {
        if self.cpu_dirty {
            if let Some(gpu) = &self.gpu {
                unsafe {
                    bindings::cuda::memcpy_to_gpu::<T>(
                        self.cpu.as_ptr(), // *const T
                        gpu.as_ptr(),      // *mut T
                        self.capacity,
                    );
                }
            }
        }
    }

    /// Synchronizes GPU data to CPU if the GPU is dirty.
    ///
    /// Currently a placeholder; actual GPU memory copy must be implemented.
    pub fn sync_to_cpu(&mut self) {
        if self.gpu_dirty {
            if let Some(gpu) = &self.gpu {
                unsafe {
                    bindings::cuda::memcpy_to_cpu(
                        self.cpu.as_mut_ptr(),    // *mut T
                        gpu.as_ptr() as *const T, // *const T
                        self.capacity,
                    );
                }
            }
        }
    }
}

impl<'mem, T> HybridBuffer<'mem, T>
where
    T: 'static + Default + NoDrop,
{
    /// Allocate the `HybridBuffer` inside an arena
    pub fn new_in(mem: &'mem Memory, len: usize, gpu: bool) -> Self {
        let data = mem.alloc_slice_fill_with(len, || T::default());
        let mut result = Self {
            gpu: if gpu { Some(GpuPtr::alloc(len)) } else { None },
            cpu: data.into_mut(),
            cpu_dirty: gpu,
            gpu_dirty: false,
            capacity: len,
        };

        result.sync_to_gpu();
        result
    }
}

// ==================================================================================
// STD implementations

/// Returns a safe immutable view of the CPU data.
/// Does not perform any synchronization with the GPU.
impl<'s, T> Deref for HybridBuffer<'s, T> {
    type Target = [T];
    fn deref(&self) -> &Self::Target {
        if self.gpu_dirty {
            println!("WARN: Access to HybridBuffer with dirty GPU memory");
        }
        self.cpu
    }
}

/// Returns a safe mutable view of the CPU data.
/// Does not perform any synchronization with the GPU.
impl<'s, T> DerefMut for HybridBuffer<'s, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        if self.gpu_dirty {
            println!("WARN: Access to HybridBuffer with dirty GPU memory");
        }
        self.cpu_dirty = true;
        self.cpu
    }
}
