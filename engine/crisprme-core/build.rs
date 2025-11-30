use std::{fs, path::Path};

// Specify the desired architecture version
const ARCH: &str = "compute_86";
const CODE: &str = "sm_86";

/// Where the ffi bridges are located
const BRIDGES: &[&str] = &[
    "src/bindings/miner.rs",
    "src/bindings/cuda.rs",
    "src/bindings/score.rs",
];

/// Cuda source code
const CUDA_SRC: &[&str] = &[
    "src/cuda/driver.cu",
    "src/cuda/scores.cu",
    //"src/cuda/miners/warp_reg_stack.cu",
    "src/cuda/miners/warp_reg_stack_nopam.cu",
    //"src/cuda/miners/shared_stack.cu"
];

fn build_on_change(folder: &str) {
    let dir = Path::new(folder);

    // Watch the directory itself (so new/removed files trigger rebuild)
    println!("cargo:rerun-if-changed={}", dir.display());

    // Recursive descent
    fn visit(dir: &Path) {
        for entry in fs::read_dir(dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            println!("cargo:rerun-if-changed={}", path.display());
            if path.is_dir() {
                visit(&path);
            }
        } 
    }

    visit(dir);
}

fn main() {
    // Set bridge file for Rust <-> CUDA comunication
    let mut cc = cxx_build::bridges(BRIDGES);

    // Setup CUDA
    let cc = cc
        .cuda(true)
        .std("c++17")
        .flag("-gencode")
        .flag(format!("arch={},code={}", ARCH, CODE))
        .flag("-m64") // for 64bit system
        .flag("-O2") // optimizations
        //.flag("-G").flag("-g")
        .flag("-cudart=shared");

    // Link CCCL and compile kernels
    cc.files(CUDA_SRC)
        .includes(&[
            // NVIDIA CCCL
            "extern/cccl/libcudacxx/include",
            "extern/cccl/thrust",
            "extern/cccl/cub",
            "include",
        ])
        .compile("crisprme-core-cuda.a");

    // Link CUDA runtime (libcudart.so)
    println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64");
    println!("cargo:rustc-link-lib=cudart");

    build_on_change("src/cuda");
    build_on_change("include");
}
