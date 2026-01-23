// threadpool.rs
use rayon::ThreadPool;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

static POOLS: OnceLock<Mutex<HashMap<usize, ThreadPool>>> = OnceLock::new();

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
