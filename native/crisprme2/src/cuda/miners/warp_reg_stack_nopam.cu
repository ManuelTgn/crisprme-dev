#include <common.cuh>
#include <assert.h>
#include <stdio.h>

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

#include <crisprme-core/src/cuda/iupac.cu>
#include <crisprme-core/src/cuda/stack.cu>

// Bindings of shared Rust/C++ structs
#include <crisprme-core/src/bindings/miner.rs.h>

#define DEBUG 0
#define SHMEM 0

#define BLOCK_COUNT 3000
// #define BLOCK_SIZE 128
#define BLOCK_SIZE 256
// #define BLOCK_COUNT 1
// #define BLOCK_SIZE 1
#define WARP_SIZE 32

#define BLOCK_WARP_COUNT ((BLOCK_SIZE + WARP_SIZE - 1) / WARP_SIZE)

namespace cg = cooperative_groups;
using warp_t = cg::thread_block_tile<32>;

/// ====================================================================
/// Global configuration

__constant__ u8 GUIDE[64];
__constant__ u8 STRAND;

__constant__ u32 GLEN;
__constant__ u32 SLEN;
__constant__ u32 PSTOP;   // == SLEN - PLEN: PAM start / protospacer end (exclusive)

__constant__ u32 GGAP;
__constant__ u32 SGAP;
__constant__ u32 MISM;

/// Total mined alignments
__device__ unsigned long long dev_total_mined = 0;
/// Current write index to alignment buffer
__device__ u32 dev_alignment_write = 0;

/// ====================================================================

/// Result alignment
/// NOTE: This must have the same representation as the Rust struct
struct Alignment
{
  Cigarx<u64> cigarx;
  u32 position;
  u8 offset;
  u8 strand;
};

/// ====================================================================

/// Sum values in the warp
/// NOTE: Only thread lane 0 has the full sum
__device__ int warp_sum(int val)
{
  unsigned mask = 0xFFFFFFFF; // full warp
  for (int offset = 16; offset > 0; offset /= 2)
  {
    val += __shfl_down_sync(mask, val, offset);
  }
  return val; // only thread lane 0 has the full sum
}

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

// Warp centric version
// NOTE: The sequence batch MUST be ordered for this to have any sense
//
//  SLEN = 10
//  GLEN = 7
//
//  starting_sidx = SLEN - GLEN - GGAP
//  end_sidx = SLEN - GLEN
//
//  0 1 2 3 4 5 6 7 8 9 P A M
//
//  S S S S S S S S S S . . .
//  G G G G G G G . . . . . .
//  . . . G G G G G G G . . .
//
//        G G G G G G G . . . / GOOD
//      G G G G G G G - . . . / GOOD
//      G G G G - G G G . . . / GOOD
//      - G G G G G G G G . . / BAD
//
__global__ void warp_reg_stack(
    const u8 *__restrict__ batch, u32 batch_size,
    Alignment *__restrict__ result, u32 capacity)
{

  // This warp
  cg::thread_block_tile<32> warp = cg::tiled_partition<32>(cg::this_thread_block());

  // Each thread of a warp has it's own stack and counters
  ThreadMiner<u64> miner;

  // TODO: This will be loaded from the state in global memory
  const u32 start_bseq = blockDim.x * blockIdx.x + threadIdx.x;
  const u32 start_sidx = 0;
  // const u32 offset = 0;

#if SHMEM
  // Contains all required sequences for this block
  // NOTE: For now use a maximum SLEN of 32, we need less
  __shared__ u8 batch_shared[BLOCK_SIZE][32];
#endif

  /// Process all sequences and starting offset
  for (u32 bseq = start_bseq; bseq < batch_size; bseq += gridDim.x * blockDim.x)
  {

#if SHMEM
    // TODO: Load data to shared memory after the barrier
    for (u32 i = 0; i < SLEN; ++i)
    {
      batch_shared[i][threadIdx.x] = batch[bseq * SLEN + i];
    }
    __syncthreads();
#endif

    // Only offsets that can satisfy the right-anchor are worth exploring:
    //   offset = pad - ggap_used + sgap_used,  ggap_used <= GGAP, sgap_used <= SGAP
    const u32 pad = PSTOP - GLEN;                        // left padding for DNA bulges
    const u32 offset_lo = pad > GGAP ? pad - GGAP : 0;
    const u32 offset_hi = pad + SGAP;

    for (u32 offset = offset_lo; offset <= offset_hi; offset += 1)
    {

      // Initialize stack for this sequence if not already done
      if (miner.mem.len == 0)
      {
#if SHMEM
        u8 s = batch_shared[threadIdx.x][offset + 0];
#else
        u8 s = batch[bseq * SLEN + offset + 0];
#endif
        u8 g = GUIDE[0];

        miner.push(Step::initial(), !iupac_match(g, s));
      }

      // Mine all toghether
      // NOTE: They will all do the same operation every time
      while (miner.has_work())
      {
        Step step = miner.current();

#if DEBUG
        printf("============================\n");
        printf("current: ");
        miner.print();
        printf("offset: %d\n", offset);
        printf("next_sidx: %d/%d\n", miner.mem.state.sidx, SLEN - 1);
        printf("next_gidx: %d/%d\n", miner.mem.state.gidx, GLEN - 1);
        printf("sgap: %d\n", miner.mem.state.sgap);
        printf("ggap: %d\n", miner.mem.state.ggap);
        printf("mism: %d\n", miner.mem.state.mism);
#endif

        // Backtrack
        if (step.is_backtrack())
        {
#if DEBUG
          printf("pop\n");
#endif

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
        bool is_complete = miner.is_complete(GLEN, PSTOP, offset);

#if DEBUG
        printf("inside_thresholds: %d\n", inside_thresholds);
        printf("can_continue: %d\n", can_continue);
        printf("is_complete: %d\n", is_complete);
#endif

        // Skip initial I
        // NOTE: I will kill whoever says that goto should not be used >:(
        if (miner.mem.state.gidx == 0 && miner.current().value == Step::Inner::S)
          goto travel;

        // Some thread has a solution
        if (warp.any(inside_thresholds && is_complete))
        {
          Cigarx<u64> cigarx = miner.cigarx();

          /*
          // Calculate how many solutions there are
          u32 solutions = (!invalid && end);
          solutions = warp_sum(solutions);
          if (warp.thread_rank() == 0) {
            u32 write_idx = atomicAdd(&dev_alignment_write, 1);
            // TODO
          }
          */

          // Add solution to alignment batch
          if (inside_thresholds && is_complete)
          {

            assert(miner.mem.state.ggap <= GGAP);
            assert(miner.mem.state.sgap <= SGAP);
            assert(miner.mem.state.mism <= MISM);
            assert(miner.mem.state.sidx + offset == PSTOP);   // <-- NEW

#if DEBUG
            char sol[64];
            cigarx.extract(sol);
            printf("solution: %s (sgap: %d, ggap: %d, mism: %d)\n", sol, miner.mem.state.sgap, miner.mem.state.ggap, miner.mem.state.mism);
#endif

            u32 write_idx = atomicAdd(&dev_alignment_write, 1);
            assert(write_idx < capacity && "ERR: too many alignments!");
            result[write_idx] = Alignment{
                cigarx,
                bseq,
                static_cast<u8>(offset),
                STRAND};
          }
        }

        // If some thread can push they all push
        if (warp.any(inside_thresholds && can_continue))
        {
#if DEBUG
          printf("push\n");
#endif

          // Guide exhausted but the alignment has not reached the PAM:
          // only DNA bulges can close the gap.
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

#if SHMEM
          u8 s = batch_shared[threadIdx.x][miner.mem.state.sidx + offset];
#else
          u8 s = batch[bseq * SLEN + miner.mem.state.sidx + offset];
#endif
          u8 g = GUIDE[miner.mem.state.gidx];

#if DEBUG
          printf("match (%c:%d vs %c:%d)? %d\n",
                 iupac_decode(s), miner.mem.state.sidx + offset, iupac_decode(g), miner.mem.state.gidx, iupac_match(g, s));
#endif
          miner.push(Step::initial(), !iupac_match(g, s));
          continue;
        }

      travel:
        // No thread can continue, travel to next operation
#if DEBUG
        printf("travel\n");
#endif
        miner.travel();
      }
    }
  }
}

namespace cuda::miner
{

  /// Invoked at the beginning of the program
  void initialize(u32 device)
  {
    CUDA_CHECK(cudaSetDevice(device));
  }

  /// Invoked before a new batch is mined
  void pre_mine(const u8 *guide, u32 glen, u32 slen, u32 plen, u32 ggap, u32 sgap, u32 mism, u8 strand)
  {

    // Copy guide to constant memory
    CUDA_CHECK(cudaMemcpyToSymbol(GUIDE, guide, glen, 0, cudaMemcpyHostToDevice));

    // Copy strand to constant memory
    CUDA_CHECK(cudaMemcpyToSymbol(STRAND, &strand, sizeof(u8), 0, cudaMemcpyHostToDevice));

    // Copy thresholds and lengths to constant memory
    CUDA_CHECK(cudaMemcpyToSymbol(SLEN, &slen, sizeof(u32), 0, cudaMemcpyHostToDevice));

    // Protospacer must end exactly where the PAM begins.
    u32 pstop = slen - plen;
    CUDA_CHECK(cudaMemcpyToSymbol(PSTOP, &pstop, sizeof(u32), 0, cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpyToSymbol(GLEN, &glen, sizeof(u32), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(SGAP, &sgap, sizeof(u32), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(GGAP, &ggap, sizeof(u32), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(MISM, &mism, sizeof(u32), 0, cudaMemcpyHostToDevice));

    // Reset alignment counter
    u32 zero = 0;
    CUDA_CHECK(cudaMemcpyToSymbol(dev_alignment_write, &zero, sizeof(u32), 0, cudaMemcpyHostToDevice));

    cudaDeviceSynchronize();
  }

  /// Mines a sequence batch and generates a single alignment batch
  MinerOutput mine(const u8 *batch, u32 batch_size, u8 *alignments, u32 capacity)
  {

    Alignment *output = reinterpret_cast<Alignment *>(alignments);
    warp_reg_stack<<<BLOCK_COUNT, BLOCK_SIZE>>>(batch, batch_size, output, capacity);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "CUDA kernel launch error in '%s' at line %d: %s (%d)\n",
              __FILE__, __LINE__, cudaGetErrorString(err), err);
    }

    u32 alignments_count = 0;
    CUDA_CHECK(cudaMemcpyFromSymbol(&alignments_count, dev_alignment_write, sizeof(u32), 0, cudaMemcpyDeviceToHost));

    // TODO: implement resume capability
    return MinerOutput{
        alignments_count,
        true};
  }

  /// Invoked after a batch has been mined
  void post_mine() {}

  /// Invoked at the end of the program
  void shutdown(u32 device) {}
}
