#pragma once

#include <common.cuh>

namespace cuda {

    /// Allocate memory on the device
    u8* malloc(u64 bytes);
    
    /// Free memory on the device
    void free(u8* memory);

    void memcpy_to_gpu(u8* gpu, const u8* cpu, u64 bytes);
    void memcpy_to_cpu(const u8* gpu, u8* cpu, u64 bytes);

    /// Pin host memory
    void pin(const u8* ptr, u64 bytes);
    /// Unpin host memory
    void unpin(const u8* ptr);

    namespace miner {
        struct MinerOutput;

        /// Invoked at the beginning of the pipeline
        void initialize(u32 device);

        /// Invoked before a new batch is mined
        void pre_mine(const u8* guide, u32 glen, u32 slen, u32 ggap, u32 sgap, u32 mism, u8 strand);

        /// Mines a sequence batch and generates a single alignment batch
        MinerOutput mine(const u8* batch, u32 batch_size, u8* alignments, u32 capacity);

        /// Invoked after a batch has been mined
        void post_mine();

        /// Invoked at the end of the pipeline
        void shutdown(u32 device);
    }
}
