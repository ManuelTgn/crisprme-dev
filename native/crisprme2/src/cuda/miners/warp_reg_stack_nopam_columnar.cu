
#include <assert.h>
#include <stdio.h>

#define CRISPRME_IMPLEMENTATION
#include <common.cuh>
#include <miner.cuh>

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

// Bindings of shared Rust/C++ structs
#include <crisprme-core/src/bindings/miner.rs.h>

// Maximum length of a guide
#define KERNEL_GUIDE_CAPACITY 32

// Maximum length of a sequence
#define KERNEL_SEQ_CAPACITY 32

/// ====================================================================
/// Constant memory
/// ====================================================================

__constant__ u8 GUIDE[KERNEL_GUIDE_CAPACITY];

__constant__ u32 GLEN;
__constant__ u32 SLEN;

// Thresholds
__constant__ u32 GGAP;
__constant__ u32 SGAP;
__constant__ u32 MISM;

/// ====================================================================
/// Warp intrinsics utils
/// ====================================================================

/*
/// Sum values in the warp
/// NOTE: Only thread lane 0 has the full sum
__device__ int warp_sum(int val)
{
    unsigned mask = 0xFFFFFFFF; // full warp
    for (int offset = 16; offset > 0; offset /= 2)
    {
        val += __shfl_down_sync(mask, val, offset);
    }
    return val;
}
*/

/// Returns true if any warp has the flag to true
__device__ bool warp_any(bool flag)
{
    unsigned mask = 0xFFFFFFFF; // full warp
    for (int offset = 16; offset > 0; offset /= 2)
    {
        flag |= __shfl_down_sync(mask, flag, offset);
    }
    // Broadcast the final result from lane 0 to all threads
    flag = __shfl_sync(mask, flag, 0);
    return flag;
}

/// ====================================================================
/// Kernel
/// ====================================================================

// A complete alignment contains 3 fields stored in 3 columns:
// - Cigarx64  (u64)
// - SeqRowIdx (u32)
// - Offset    ( u8)

/// Current write index to alignment buffer
__device__ u32 write_cursor = 0;

__global__ void kmine(
    const u8 *__restrict__ sequences, // Input sequences
    u64 *__restrict__ cigarxs,        // Output cigarx column
    u32 *__restrict__ indexes,        // Output index column
    u8 *__restrict__ offsets,         // Output offset column
    u32 sequence_count,               // Sequence count
    u32 capacity                      // Output capacity
)
{
    // This warp
    namespace cg = cooperative_groups;
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(
        cg::this_thread_block());

    // Each thread of a warp has it's own stack and counters
    ThreadMiner<u64> miner;

    // TODO: This will be loaded from the state in global memory
    //       when the resumable kernel feature is done...
    const auto start_bseq = blockDim.x * blockIdx.x + threadIdx.x;
    const auto start_sidx = 0;

    /// Process all sequences and starting offset
    const auto stride = gridDim.x * blockDim.x;
    for (u32 bseq = start_bseq; bseq < sequence_count; bseq += stride)
    {
        for (u8 offset = start_sidx; offset <= (SLEN - GLEN + GGAP); offset += 1)
        {
            // Initialize stack for this sequence if not already done
            if (miner.mem.len == 0)
            {
                u8 s = sequences[bseq * KERNEL_SEQ_CAPACITY + offset + 0];
                u8 g = GUIDE[0];

                miner.push(Step::initial(), !iupac_match(g, s));
            }

            // Mine all toghether
            // NOTE: They will all do the same operation every time
            while (miner.has_work())
            {
                Step step = miner.current();

                // Backtrack
                if (step.is_backtrack())
                {
                    miner.pop();
                    if (miner.has_work())
                    {
                        // NOTE: This is necessary to exit an infinite loop
                        miner.travel();
                    }
                    continue;
                }

                // Controllers, these depends on input data
                bool inside_thresholds = miner.inside_thresholds(MISM, GGAP, SGAP);
                bool can_continue = miner.can_continue(SLEN, GLEN, offset, GGAP);
                bool is_complete = miner.is_complete(GLEN);

                // Skip initial I
                // NOTE: I will kill whoever says that goto should not be used >:(
                if (miner.mem.state.gidx == 0 && miner.current().value == Step::Inner::S)
                    goto travel;

                // Some thread has a solution
                if (warp.any(inside_thresholds && is_complete))
                {
                    // Cigarx64 cigarx = miner.cigarx();

                    // Add solution to alignment batch
                    if (inside_thresholds && is_complete)
                    {
                        assert(miner.mem.state.ggap <= GGAP);
                        assert(miner.mem.state.sgap <= SGAP);
                        assert(miner.mem.state.mism <= MISM);

                        // Write all output columns
                        u32 write_idx = atomicAdd(&write_cursor, 1);
                        if (write_idx < capacity)
                        {
                            cigarxs[write_idx] = miner.cigarx().storage();
                            indexes[write_idx] = bseq;
                            offsets[write_idx] = offset;
                        }
                        else
                        {
                            assert(false && "ERR: too many alignments!");
                        }
                    }
                }

                // If some thread can push they all push
                if (warp.any(inside_thresholds && can_continue))
                {
                    // If the target sequence is out of bound this means that we can only add deletions
                    // NOTE: The mismatch flag is not needed as the deletion is the last type of step before
                    // a backtrack, it will never be used.
                    if (miner.mem.state.sidx + offset >= SLEN)
                    {
                        miner.push(Step::deletion(), false);
                        continue;
                    }

                    // In the other cases we can proceed as normal and push a match/mismatch step
                    u8 s = sequences[bseq * KERNEL_SEQ_CAPACITY + miner.mem.state.sidx + offset];
                    u8 g = GUIDE[miner.mem.state.gidx];

                    miner.push(Step::initial(), !iupac_match(g, s));
                    continue;
                }

            travel:
                // No thread can continue, travel to next operation
                miner.travel();
            }
        }
    }
}

/// ====================================================================
/// Host
/// ====================================================================

// Module entry point
namespace cuda::miner
{
    const u32 ZERO = 0;

    void initialize(u32 device)
    {
        // Make all CUDA calls locat to the calling thread
        CUDA_CHECK(cudaSetDevice(device));
    }

    // Load all meta parameters in constant memory
    void prepare(MinerConfig config)
    {
        // Copy guide to constant memory
        CUDA_CHECK(cudaMemcpyToSymbol(GUIDE, config.guide, config.glen, 0, cudaMemcpyHostToDevice));

        // Reset write cursors
        CUDA_CHECK(cudaMemcpyToSymbol(write_cursor, &ZERO, sizeof(u32), 0, cudaMemcpyHostToDevice));

        // Copy thresholds and lengths to constant memory
        CUDA_CHECK(cudaMemcpyToSymbol(SLEN, &config.slen, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(GLEN, &config.glen, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(SGAP, &config.sgap, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(GGAP, &config.ggap, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(MISM, &config.mism, sizeof(u32), 0, cudaMemcpyHostToDevice));
    }

    // Launch columnar miner kernel
    MinerOutput launch(MinerInput input)
    {
        u32 BLOCK_SIZE = 256;
        u32 BLOCK_COUNT = (input.seq_count + BLOCK_SIZE - 1) / BLOCK_SIZE;

        // Launch kernel and pray
        kmine<<<BLOCK_COUNT, BLOCK_SIZE>>>(
            input.sequences, input.cigarx, input.index, input.offset,
            input.seq_count, input.capacity);

        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            fprintf(stderr, "CUDA kernel launch error in '%s' at line %d: %s (%d)\n",
                    __FILE__, __LINE__, cudaGetErrorString(err), err);
        }

        // Get total mined alignments
        u32 mined = 0;
        CUDA_CHECK(cudaMemcpyFromSymbol(&mined, write_cursor, sizeof(u32), 0, cudaMemcpyDeviceToHost));
        return MinerOutput{
            mined,
            true};
    }

    void post_mine() { }

    void shutdown(u32 device) { }
}