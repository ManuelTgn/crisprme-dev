
#include <assert.h>
#include <stdio.h>

#define CRISPRME_IMPLEMENTATION
#include <common.cuh>
#include <miner.cuh>

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cuda/atomic>

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
__constant__ u32 PSTOP;   // == SLEN - PLEN: PAM start / protospacer end (exclusive)

// Thresholds
__constant__ u32 GGAP;
__constant__ u32 SGAP;
__constant__ u32 MISM;

/// ====================================================================
/// Checkpoint
/// ====================================================================

template<typename Storage>
struct ThreadCheckpoint {

    // Inner state of the miner
    typename ThreadMiner<Storage>::Inner inner;
    // Current sequence
    u32 bseq;
    // Current offset
    u8 offset;
};

/// ====================================================================
/// Warp intrinsics utils
/// ====================================================================

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
/// Output buffer is full, ignores cache
__device__ volatile u32 output_full = 0;

__global__ void kmine(
    const u8 *__restrict__ sequences,   // Input sequences
    u64 *__restrict__ cigarxs,          // Output cigarx column
    u32 *__restrict__ indexes,          // Output index column
    u8 *__restrict__ offsets,           // Output offset column
    u32 sequence_count,                 // Sequence count
    u32 capacity,                       // Output capacity
    ThreadCheckpoint<u64> *checkpoints  // Per-thread checkpoint state
) {
    // This warp
    namespace cg = cooperative_groups;
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(
        cg::this_thread_block());

    // Each thread of a warp has it's own stack and counters
    ThreadMiner<u64> miner;

    const u32 tid = blockDim.x * blockIdx.x + threadIdx.x;
    const u32 stride = gridDim.x * blockDim.x;

    // Restore from checkpoint if this thread was interrupted
    ThreadCheckpoint<u64> cp = checkpoints[tid];

    u32 start_bseq;
    u8  start_offset;

    if (cp.inner.len > 0)
    {
        start_bseq   = cp.bseq;
        start_offset = cp.offset;
        miner.mem    = cp.inner;
    }
    else
    {
        start_bseq   = tid;
        start_offset = 0;
    }

    // Only offsets that can satisfy the right-anchor are worth exploring:
    //   offset = pad - ggap_used + sgap_used,  ggap_used <= GGAP, sgap_used <= SGAP
    const u32 pad = PSTOP - GLEN;                        // left padding for DNA bulges
    const u32 offset_lo = pad > GGAP ? pad - GGAP : 0;
    const u32 offset_hi = pad + SGAP;

    /// Process all sequences and starting offset
    for (u32 bseq = start_bseq; bseq < sequence_count; bseq += stride)
    {
        // A resumed thread must continue from the checkpointed offset only for the
        // interrupted sequence. All later sequences need to restart from offset_lo.
        const u8 seq_start_offset = (bseq == start_bseq) ? start_offset : (u8)offset_lo;
        for (u8 offset = seq_start_offset; offset <= offset_hi; offset += 1)
        {
            // Check if another thread signalled output full
            if (output_full)
            {
                checkpoints[tid] = { miner.mem, bseq, offset };
                return;
            }

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
                bool can_continue = miner.can_continue(GLEN, PSTOP, offset);
                bool is_complete  = miner.is_complete(GLEN, PSTOP, offset);

                // Skip initial I
                // NOTE: I will kill whoever says that goto should not be used >:(
                if (miner.mem.state.gidx == 0 && miner.current().value == Step::Inner::S)
                    goto travel;

                // Some thread has a solution
                if (warp.any(inside_thresholds && is_complete))
                {
                    // Add solution to alignment batch
                    if (inside_thresholds && is_complete)
                    {
                        assert(miner.mem.state.ggap <= GGAP);
                        assert(miner.mem.state.sgap <= SGAP);
                        assert(miner.mem.state.mism <= MISM);
                        assert(miner.mem.state.sidx + offset == PSTOP);   // <-- NEW

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
                            // Output buffer is full — clamp cursor and signal
                            atomicMin(&write_cursor, capacity);
                            atomicExch((u32*)&output_full, 1);
                            // Save state without advancing so this result
                            // is re-discovered on the next launch
                            checkpoints[tid] = { miner.mem, bseq, offset };
                            return;
                        }
                    }
                }

                // If some thread can push they all push
                if (warp.any(inside_thresholds && can_continue))
                {
                    // Guide exhausted but the alignment has not reached the PAM:
                    // only DNA bulges can close the gap. Push directly —
                    // Step::initial() is a B step and would advance gidx past GLEN.
                    if (miner.mem.state.gidx >= GLEN)
                    {
                        miner.push(Step::dna_bulge(), false);
                        continue;
                    }

                    // DNA cursor is at the PAM but the guide is not exhausted:
                    // only RNA bulges remain. Note the bound is PSTOP, not SLEN —
                    // the guide must never align into the PAM.
                    if (miner.mem.state.sidx + offset >= PSTOP)
                    {
                        miner.push(Step::rna_bulge(), false);
                        continue;
                    }

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

    // Thread completed all work — clear checkpoint
    checkpoints[tid].inner.len = 0;
}

/// ====================================================================
/// Host
/// ====================================================================

// Module entry point
namespace cuda::miner
{
    const u32 ZERO = 0;

    // Persistent checkpoint buffer across launches within a prepare/post_mine session
    static ThreadCheckpoint<u64>* d_checkpoints = nullptr;
    static u32 checkpoint_capacity = 0;

    void initialize(u32 device)
    {
        // Make all CUDA calls local to the calling thread
        CUDA_CHECK(cudaSetDevice(device));
    }

    // Load all meta parameters in constant memory
    void prepare(MinerConfig config)
    {
        // Copy guide to constant memory
        CUDA_CHECK(cudaMemcpyToSymbol(GUIDE, config.guide, config.glen, 0, cudaMemcpyHostToDevice));

        // Reset write cursor and output_full flag
        CUDA_CHECK(cudaMemcpyToSymbol(write_cursor, &ZERO, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(output_full,  &ZERO, sizeof(u32), 0, cudaMemcpyHostToDevice));

        // Copy thresholds and lengths to constant memory
        CUDA_CHECK(cudaMemcpyToSymbol(SLEN, &config.slen, sizeof(u32), 0, cudaMemcpyHostToDevice));

        // Protospacer must end exactly where the PAM begins.
        u32 pstop = config.slen - config.plen;
        CUDA_CHECK(cudaMemcpyToSymbol(PSTOP, &pstop, sizeof(u32), 0, cudaMemcpyHostToDevice));

        CUDA_CHECK(cudaMemcpyToSymbol(GLEN, &config.glen, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(SGAP, &config.sgap, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(GGAP, &config.ggap, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(MISM, &config.mism, sizeof(u32), 0, cudaMemcpyHostToDevice));
    }

    // Launch columnar miner kernel
    MinerOutput launch(MinerInput input)
    {
        u32 BLOCK_SIZE = 256;
        u32 BLOCK_COUNT = (input.seq_count + BLOCK_SIZE - 1) / BLOCK_SIZE; // 3000
        u32 TH_COUNT = BLOCK_SIZE * BLOCK_COUNT;

        // Lazy-allocate checkpoint buffer (persists across launches, freed in post_mine)
        if (d_checkpoints == nullptr || checkpoint_capacity < TH_COUNT)
        {
            printf("allocated CUDA checkpoint buffer (%d bytes)\n", 
                TH_COUNT * sizeof(ThreadCheckpoint<u64>));

            if (d_checkpoints) 
                CUDA_CHECK(cudaFree(d_checkpoints));

            CUDA_CHECK(cudaMalloc(&d_checkpoints,    TH_COUNT * sizeof(ThreadCheckpoint<u64>)));
            CUDA_CHECK(cudaMemset(d_checkpoints,  0, TH_COUNT * sizeof(ThreadCheckpoint<u64>)));
            checkpoint_capacity = TH_COUNT;
        }

        // Reset write cursor and output_full flag before each launch
        CUDA_CHECK(cudaMemcpyToSymbol(write_cursor, &ZERO, sizeof(u32), 0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpyToSymbol(output_full,  &ZERO, sizeof(u32), 0, cudaMemcpyHostToDevice));

        kmine<<<BLOCK_COUNT, BLOCK_SIZE>>>(
            input.sequences, input.cigarx, input.index, input.offset,
            input.seq_count, input.capacity, d_checkpoints);

        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            fprintf(stderr, "CUDA kernel launch error in '%s' at line %d: %s (%d)\n",
                    __FILE__, __LINE__, cudaGetErrorString(err), err);
        }

        // Read back results
        u32 mined = 0, full = 0;
        CUDA_CHECK(cudaMemcpyFromSymbol(&mined, write_cursor, sizeof(u32), 0, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpyFromSymbol(&full,  output_full,  sizeof(u32), 0, cudaMemcpyDeviceToHost));

        return MinerOutput{
            mined,
            full == 0};
    }

    void post_mine() { }

    void shutdown(u32 device) 
    {
        // Free checkpoint buffer at end of mining session
        if (d_checkpoints)
        {
            printf("free CUDA checkpoint buffer\n");
            CUDA_CHECK(cudaFree(d_checkpoints));
            d_checkpoints = nullptr;
            checkpoint_capacity = 0;
        }
    }
}
