fn main() {
    // Only run the heavy CUDA/C++ build when the feature is enabled.
    let cuda_enabled = std::env::var("CARGO_FEATURE_CUDA").is_ok();
    if !cuda_enabled {
        // Still re-run if these change, so switching features is reliable
        println!("cargo:rerun-if-changed=build.rs");
        return;
    }

    // Optional: if nvcc is missing, fail gracefully (or warn and skip).
    let nvcc_ok = which::which("nvcc").is_ok();
    if !nvcc_ok {
        // Prefer *not* failing for analysis; just skip building CUDA.
        println!("cargo:warning=nvcc not found; skipping CUDA build (feature=cuda)");
        return;
    }

    // Your real CUDA build steps here (cxx_build, cmake, etc.)
    // ...
}
