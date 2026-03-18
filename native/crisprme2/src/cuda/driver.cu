#include <common.cuh>

// Bindings of shared Rust/C++ structs
#include <crisprme-core/src/bindings/cuda.rs.h>

namespace cuda {
    
    /// Allocate memory on the device
    u8* malloc(u64 bytes) {
        void* result = nullptr;
        cudaMalloc(&result, bytes);
        return (u8*)result;
    }

    /// Free memory on the device
    void free(u8* memory) {
        cudaFree((void*)memory);
    }

    void memcpy_to_gpu(u8* gpu, const u8* cpu, u64 bytes) {
        cudaMemcpy((void*)gpu, (const void*)cpu, bytes, cudaMemcpyHostToDevice);
    }

    void memcpy_to_cpu(const u8* gpu, u8* cpu, u64 bytes) {
        cudaMemcpy((void*)cpu, (const void*)gpu, bytes, cudaMemcpyDeviceToHost);
    }

    void pin(const u8* ptr, u64 bytes) {
        CUDA_CHECK(cudaHostRegister((void*)ptr, bytes, 0));
    }

    void unpin(const u8* ptr) {
        CUDA_CHECK(cudaHostUnregister((void*)ptr));
    }
}
