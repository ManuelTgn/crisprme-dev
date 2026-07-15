/// Thread pool utilities for efficient and reproducible parallel execution.
///
/// This module provides a lightweight caching layer on top of Rayon’s
/// `ThreadPool`, allowing thread pools to be reused across multiple calls
/// with the same number of threads.
///
/// Caching avoids the significant overhead of repeatedly constructing
/// and tearing down Rayon thread pools during high-frequency operations
/// such as genome-wide sequence scanning.
use rayon::ThreadPool;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

/// Global cache of Rayon thread pools indexed by thread count.
///
/// Each distinct `threads` value corresponds to exactly one `ThreadPool`,
/// which is lazily initialized on first use and reused thereafter.
///
/// The cache is protected by a `Mutex` and initialized via `OnceLock`
/// to ensure thread-safe, one-time construction.
static POOLS: OnceLock<Mutex<HashMap<usize, ThreadPool>>> = OnceLock::new();

/// Executes a closure inside a cached Rayon thread pool with a fixed
/// number of threads.
///
/// If a thread pool with the requested number of threads already exists,
/// it is reused. Otherwise, a new pool is created, cached, and used.
/// The provided closure is executed using `ThreadPool::install`, ensuring
/// that all Rayon parallel iterators spawned within the closure use the
/// specified pool.
///
/// # Arguments
/// * `threads` - Number of worker threads to use. Must be greater than zero.
/// * `f` - Closure to execute inside the thread pool.
///
/// # Returns
/// * `Ok(T)` containing the result of the closure on success
/// * `Err(String)` if the thread count is invalid or the pool cache is poisoned
///
/// # Errors
/// Returns an error if:
/// * `threads == 0`
/// * The internal thread pool cache mutex is poisoned
///
/// # Panics
/// Panics if Rayon fails to build a thread pool. This is considered unrecoverable
/// and indicates a serious system-level failure.
///
/// # Concurrency Notes
/// * Thread pools are shared across calls but never modified after creation.
/// * Multiple calls with different `threads` values may execute concurrently
///   using distinct pools.
/// * The cache lookup is serialized via a mutex, but pool execution is fully
///   parallel and lock-free.
///
/// # Rationale
/// Creating a Rayon `ThreadPool` is relatively expensive. In workloads such as
/// CRISPRme2 genome scanning, where many short-lived parallel regions are invoked,
/// caching thread pools provides a measurable performance benefit while retaining
/// deterministic thread usage.
pub fn with_pool<T, F>(threads: usize, f: F) -> Result<T, String>
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    if threads == 0 {
        return Err("threads must be > 0".to_string());
    }

    let pools = POOLS.get_or_init(|| Mutex::new(HashMap::new()));
    let mut guard = pools
        .lock()
        .map_err(|_| "ThreadPool cache lock poisoned".to_string())?;

    let pool = guard.entry(threads).or_insert_with(|| {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build Rayon ThreadPool")
    });

    Ok(pool.install(f))
}
