#include <common.cuh>
#include <assert.h>
#include <stdio.h>

#include <crisprme-core/src/cuda/iupac.cu>
#include <crisprme-core/src/cuda/step.cu>

// Bindings of shared Rust/C++ structs
#include <crisprme-core/src/bindings/miner.rs.h>

#define BLOCK_COUNT 3000
#define BLOCK_SIZE 128
#define TH_STACK 7

/// ====================================================================
/// Global configuration

__constant__ u8 GUIDE[64];

__constant__ u32 GLEN;
__constant__ u32 SLEN;

__constant__ u32 GGAP;
__constant__ u32 SGAP;
__constant__ u32 MISM;
__constant__ u8 STATES_SEQ;

/// Total mined alignments
__device__ unsigned long long dev_total_mined = 0;
/// Current write index to alignment buffer
__device__ u32 dev_alignment_write = 0;
/// Current global sequence index
__device__ u32 dev_idx = BLOCK_COUNT*BLOCK_SIZE;

/// ====================================================================

/// Result alignment
/// NOTE: This must have the same representation as the Rust struct
struct Alignment {
  Cigarx<u64> cigarx;
  u32 position;
  u8 offset;
};

struct State {
    u8 qidx;
    u8 tidx;
    u8 qgap : 4;
    u8 tgap : 4;
    u8 mgap;
    u32 soluz1;
    u32 soluz2;
};

struct State_ {
    u32 qidx : 6;
    u32 tidx : 6;
    u32 qgap : 5;
    u32 tgap : 5;
    u32 mgap : 6;
    u32 soluz1;
    u32 soluz2;
};

__global__ void alignments_array(const u8* T, u32 tot_seq, Alignment* soluz){
    __shared__ State sh_stack[BLOCK_SIZE*TH_STACK];
    u8 i_sh_stack = 0;
    u32 idx = (threadIdx.x + (blockIdx.x * blockDim.x));// * states_th;
    State current_state = (idx < tot_seq) ? State{0, (u8)(idx%STATES_SEQ), 0, 0, 0, 0} : State{0};

    //while (current_state.qidx < UINT8_MAX){
    while (idx < tot_seq){

      //current_state.soluz2 = (current_state.soluz2<<2) + (current_state.soluz1>>30);
      current_state.soluz2 += current_state.soluz1>>30;
      current_state.soluz1 = current_state.soluz1<<2;

      if (current_state.tidx < SLEN){
        // Match - Mismatch
        u8 next_mism = !(GUIDE[current_state.qidx] & T[(idx/STATES_SEQ)*SLEN + current_state.tidx]);
        if (current_state.mgap + next_mism <= MISM){
          if ((current_state.qidx+1 < GLEN))
            sh_stack[(threadIdx.x*TH_STACK)+i_sh_stack++] = State{
              u8(current_state.qidx + 1),
              u8(current_state.tidx + 1),
              current_state.qgap,
              current_state.tgap,
              u8(current_state.mgap + next_mism),
              current_state.soluz1 + next_mism,
              current_state.soluz2<<2
            };
          else
            soluz[atomicInc(&dev_alignment_write, UINT32_MAX)] = {
              {
                (uint64_t(current_state.soluz2)<<32) + current_state.soluz1 + next_mism,
                u8((current_state.qidx + current_state.qgap + 1)<<1)
              },
              idx/STATES_SEQ,
              u8(idx%STATES_SEQ)
            };
        }

        // Gap in query (deletion in target)
        if (current_state.qgap < GGAP)
          sh_stack[(threadIdx.x*TH_STACK)+i_sh_stack++] = State{
            current_state.qidx,
            u8(current_state.tidx + 1),
            u8(current_state.qgap + 1),
            current_state.tgap,
            current_state.mgap,
            current_state.soluz1 + 3,
            current_state.soluz2<<2
          };
      }

      // Gap in target (insertion in query)
      if (current_state.tgap < SGAP){
        current_state.soluz1 += 2;
        current_state.tgap++;
        current_state.soluz2 = current_state.soluz2<<2;
        if (++current_state.qidx == GLEN){
          soluz[atomicInc(&dev_alignment_write, UINT32_MAX)] = {
            {
              (uint64_t(current_state.soluz2)<<30) + current_state.soluz1,
              u8((current_state.qidx + current_state.qgap)<<1)
            },
            idx/STATES_SEQ,
            u8(idx%STATES_SEQ)
          };
          if (i_sh_stack) current_state = sh_stack[(threadIdx.x*TH_STACK)+(--i_sh_stack)];
          else if ((idx = atomicInc(&dev_idx, UINT32_MAX)) < tot_seq) current_state = State{0, u8(idx%STATES_SEQ), 0, 0, 0, 0};
        }
      }
      else if (i_sh_stack) current_state = sh_stack[(threadIdx.x*TH_STACK)+(--i_sh_stack)];
      else if ((idx = atomicInc(&dev_idx, UINT32_MAX)) < tot_seq) current_state = State{0, u8(idx%STATES_SEQ), 0, 0, 0, 0};
    }
}

namespace cuda::miner {
u8 s_seq;

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
    CUDA_CHECK(cudaMemcpyToSymbol(SGAP, &sgap, sizeof(u8), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(GGAP, &ggap, sizeof(u8), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(MISM, &mism, sizeof(u32), 0, cudaMemcpyHostToDevice));

    s_seq = slen-(glen-(sgap+1));
    CUDA_CHECK(cudaMemcpyToSymbol(STATES_SEQ, &s_seq, sizeof(u8), 0, cudaMemcpyHostToDevice));

    // Reset alignment counter
    u32 zero = 0;
    CUDA_CHECK(cudaMemcpyToSymbol(dev_alignment_write, &zero, sizeof(u32), 0, cudaMemcpyHostToDevice));

    u32 val_idx = BLOCK_COUNT * BLOCK_SIZE;
    CUDA_CHECK(cudaMemcpyToSymbol(dev_idx, &val_idx, sizeof(u32), 0, cudaMemcpyHostToDevice));

    cudaDeviceSynchronize();
  }

  /// Mines a sequence batch and generates a single alignment batch
  MinerOutput mine(const u8* batch, u32 batch_size, u8* alignments, u32 capacity) {

    Alignment* output = reinterpret_cast<Alignment*>(alignments);

    //u32 states_th = (batch_size*s_seq)/(BLOCK_COUNT*BLOCK_SIZE) + (((batch_size*s_seq)%(BLOCK_COUNT*BLOCK_SIZE)) ? 1 : 0);
    u32 tot_seq = batch_size*s_seq;

    //printf("---\tstates_th: %i\tTlen: %i\tstate_seq: %i\n", states_th, batch_size, s_seq);

    alignments_array<<<BLOCK_COUNT, BLOCK_SIZE>>>(batch, tot_seq, output);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, "CUDA kernel launch error in '%s' at line %d: %s (%d)\n",
                    __FILE__, __LINE__, cudaGetErrorString(err), err);
    }

    u32 alignments_count = 0;
    CUDA_CHECK(cudaMemcpyFromSymbol(&alignments_count, dev_alignment_write, sizeof(u32), 0, cudaMemcpyDeviceToHost));
    //printf("+++\tmined: %i\t%f%\n", alignments_count, (alignments_count*100.0)/capacity);

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