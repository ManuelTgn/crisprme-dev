#include <common.cuh>
#include <assert.h>
#include <stdio.h>

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

#include <crisprme-core/src/cuda/iupac.cu>
#include <crisprme-core/src/cuda/step.cu>

// Bindings of shared Rust/C++ structs
#include <crisprme-core/src/bindings/miner.rs.h>

#define DEBUG 0

#define BLOCK_COUNT 3000
#define BLOCK_SIZE 128
//#define BLOCK_COUNT 1
//#define BLOCK_SIZE 1
#define WARP_SIZE 32

#define BLOCK_WARP_COUNT ((BLOCK_SIZE + WARP_SIZE - 1) / WARP_SIZE)

namespace cg = cooperative_groups;
using warp_t = cg::thread_block_tile<32>;

/// ====================================================================
/// Global configuration

__constant__ u8 GUIDE[64];

__constant__ u32 GLEN;
__constant__ u32 SLEN;

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
struct Alignment {
  Cigarx<u64> cigarx;
  u32 position;
  u8 offset;
};

/// State of a single thread
/// NOTE: This uses 8 registers
struct ThreadState {

  /// Contains the current exploration path
  /// NOTE: This uses 1+4 registers
  StepStack128 stack;

  /// Next location for both guide and reference
  u32 next_gidx;
  u32 next_sidx;

  /// Counters
  u32 curr_mism;
  u32 curr_ggap;
  u32 curr_sgap;

  /// Remove top step and update counters
  __device__ __forceinline__ void pop() {
    dec(stack.current());
    stack.pop();
  }

  /// Add a new step and update counters
  __device__ __forceinline__ void push(Step step) {
    stack.push(step);
    inc(step);
  }

  /// Go the next operation of top step and update counters
  __device__ __forceinline__ void travel() {
    Step step = stack.current();
    dec(step);

    step = step.next();
    stack.replace(step);

    inc(step);
  }

  /// Check if state is invalid
  __device__ __forceinline__ bool invalid(u32 slen, u32 max_mism, u32 max_ggap, u32 max_sgap) {
    return (next_sidx >= slen || curr_mism > max_mism || curr_ggap > max_ggap || curr_sgap > max_sgap);
  }

  /// Check if state is inside the thresholds
  __device__ __forceinline__ bool inside_thresholds(u32 slen, u32 max_mism, u32 max_ggap, u32 max_sgap) {
    return (curr_mism <= max_mism && curr_ggap <= max_ggap && curr_sgap <= max_sgap);
  }

  /// Check if we can proceed to the next index
  __device__ __forceinline__ bool can_continue(u32 slen, u32 glen, u32 offset) {
    return (next_sidx + offset < slen || next_gidx < glen);
  }

  /// Check if state is final
  __device__ __forceinline__ bool can_be_solution(u32 glen, u32 slen) {
    return (next_gidx == glen && next_sidx <= slen);
  }

  /// Check if state is final
  __device__ __forceinline__ bool final(u32 glen) {
    return (next_gidx == glen);
  }

  /// Decrement counters based on the step
  __device__ __forceinline__ void dec(Step step) {
    next_gidx -= step.gidx_dt();
    next_sidx -= step.sidx_dt();
    curr_mism -= step.mism_dt();
    curr_ggap -= step.ggap_dt();
    curr_sgap -= step.sgap_dt();
  }

  /// Increment counters based on the step
  __device__ __forceinline__ void inc(Step step) {
    next_gidx += step.gidx_dt();
    next_sidx += step.sidx_dt();
    curr_mism += step.mism_dt();
    curr_ggap += step.ggap_dt();
    curr_sgap += step.sgap_dt();
  }
};

/// ====================================================================

/// Sum values in the warp
/// NOTE: Only thread lane 0 has the full sum
__device__ int warp_sum(int val) {
    unsigned mask = 0xFFFFFFFF; // full warp
    for (int offset = 16; offset > 0; offset /= 2) {
        val += __shfl_down_sync(mask, val, offset);
    }
    return val; // only thread lane 0 has the full sum
}

/// Returns true if any warp has the flag to true
__device__ bool warp_any(bool flag) {
    unsigned mask = 0xFFFFFFFF; // full warp
    for (int offset = 16; offset > 0; offset /= 2) {
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
//  S S S S ? // This is valid only if we have a sequence gap
//  G G G - G
//
__global__ void warp_reg_stack(
  const u8* __restrict__ batch, u32 batch_size,
  Alignment* __restrict__ result, u32 capacity
) {

  // This warp
  cg::thread_block_tile<32> warp = cg::tiled_partition<32>(cg::this_thread_block());

  // Each thread of a warp has it's own stack and counters
  ThreadState state = { };

  // TODO: This will be loaded from the state in global memory
  const u32 start_bseq = blockDim.x * blockIdx.x + threadIdx.x;
  const u32 start_sidx = 0;
  //const u32 offset = 0;

  /// Process all sequences and starting offset
  for (u32 bseq = start_bseq; bseq < batch_size; bseq += gridDim.x * blockDim.x) {
    // TODO: Here the ending offset must take into consideration the gaps allowed! (for now +3)
    for (u32 offset = start_sidx; offset <= (SLEN - GLEN + GGAP + SGAP); offset += 1) {

      // Initialize stack for this sequence if not already done
      if (state.stack.len() == 0) {
        u8 s = batch[bseq * SLEN + offset + 0];
        u8 g = GUIDE[0];

        Step initial = Step::initial(iupac_match(g, s));
        state.push(initial);
      }
    
      // Mine all toghether
      // NOTE: They will all do the same operation every time
      while(state.stack.len() != 0) {
        Step step = state.stack.current();

#if DEBUG
        printf("============================\n");
        printf("current: "); state.stack.print();
	printf("offset: %d\n", offset);
        printf("next_sidx: %d\n", state.next_sidx);
        printf("next_gidx: %d\n", state.next_gidx);
        printf("sgap: %d\n", state.curr_sgap);
        printf("ggap: %d\n", state.curr_ggap);
        printf("mism: %d\n", state.curr_mism);
#endif

        // Backtrack
        if (step.is_backtrack()) {
#if DEBUG
          printf("pop\n");
#endif

          state.pop();
          if (state.stack.len() != 0) {
            state.travel();
          }
          continue;
        }

        // Controllers, these depends on input data
        bool inside_thresholds = state.inside_thresholds(SLEN, MISM, GGAP, SGAP);
        bool can_continue = state.can_continue(SLEN, GLEN, offset);
        bool can_be_solution = state.can_be_solution(GLEN, SLEN);

#if DEBUG
	printf("inside_thresholds: %d\n", inside_thresholds);
	printf("can_continue: %d\n", can_continue);
	printf("can_be_solution: %d\n", can_be_solution);
#endif

        // Skip initial I
        // NOTE: I will kill whoever says that goto should not be used >:(
	if (state.next_gidx == 0 && state.stack.current().to_cigarx_op() == Cigarx<u64>::Inner::I)
          goto travel;

	/*
	// Skip match/mismatch out-of-bound
        // NOTE: I will kill whoever says that goto should not be used >:(
	if (state.next_sidx + offset >= SLEN && cur_operation != Cigarx<u64>::Inner::D) {
#if DEBUG
	  printf("skipped invalid step\n");
#endif
	  goto travel;
	}
	*/

        // Some thread has a solution
        if (warp.any(inside_thresholds && can_be_solution)) {
          Cigarx<u64> cigarx = state.stack.cigarx();

#if DEBUG
          printf("solution\n");
          cigarx.print();
#endif

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
          if (inside_thresholds) {
            u32 write_idx = atomicAdd(&dev_alignment_write, 1);
            assert(write_idx < capacity && "ERR: too many alignments!");
            result[write_idx] = Alignment {
              cigarx, bseq, static_cast<u8>(offset)
            };
          }

        }

        // If some thread can still push they all push
        if (warp.any(inside_thresholds && can_continue)) {
#if DEBUG
          printf("push\n");
#endif
          
	  // If the target sequence is out-of-bounds we can only add deletions
	  if (state.next_sidx + offset >= SLEN) {
#if DEBUG
	  	printf("out-of-bound, push deletion\n");
#endif
		step = Step::deletion();
		state.push(step);
		continue;
	  }

	  // Otherwise we can continue as normal
          assert(state.next_gidx < GLEN);
          assert(state.next_sidx + offset < SLEN);

	  u8 s = batch[bseq * SLEN + state.next_sidx + offset];
          u8 g = GUIDE[state.next_gidx];

#if DEBUG
	  printf("match (%c:%d vs %c:%d)? %d\n", 
	      iupac_decode(s), state.next_sidx + offset, iupac_decode(g), state.next_gidx, iupac_match(g, s));
#endif
          step = Step::initial(iupac_match(g, s));
          state.push(step);
          continue;
        }

  travel:
        // No thread can continue, travel to next operation
#if DEBUG
	printf("travel\n");
#endif
        state.travel();
      }
    }
  }
}

namespace cuda::miner {

  /// Invoked at the beginning of the program
  void initialize(u32 device) {
    CUDA_CHECK(cudaSetDevice(device)); 
  }

  /// Invoked before a new batch is mined
  void pre_mine(const u8* guide, u32 glen, u32 slen, u32 ggap, u32 sgap, u32 mism) {
  
    // Copy guide to constant memory
    CUDA_CHECK(cudaMemcpyToSymbol(GUIDE, guide, glen, 0, cudaMemcpyHostToDevice));

    // Copy thresholds and lengths to constant memory
    CUDA_CHECK(cudaMemcpyToSymbol(SLEN, &slen, sizeof(u32), 0, cudaMemcpyHostToDevice));
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
  MinerOutput mine(const u8* batch, u32 batch_size, u8* alignments, u32 capacity) {
    
    Alignment* output = reinterpret_cast<Alignment*>(alignments);
    warp_reg_stack<<<BLOCK_COUNT, BLOCK_SIZE>>>(batch, batch_size, output, capacity);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();                                
    if (err != cudaSuccess) {                                            
      fprintf(stderr, "CUDA kernel launch error in '%s' at line %d: %s (%d)\n",
                    __FILE__, __LINE__, cudaGetErrorString(err), err);
    }

    u32 alignments_count = 0;
    CUDA_CHECK(cudaMemcpyFromSymbol(&alignments_count, dev_alignment_write, sizeof(u32), 0, cudaMemcpyDeviceToHost));
    
    // TODO: implement resume capability
    return MinerOutput {
      alignments_count,
      true
    };
  }

  /// Invoked after a batch has been mined
  void post_mine() { }

  /// Invoked at the end of the program
  void shutdown(u32 device) { }
}

