#include "crisprme-core/include/scores.cuh"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#define DEBUG_STRINGS 0
#define DEBUG_QUERY 0

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

constexpr u32 warp_size = 32;
constexpr u32 cost_gaps = 1;


/// Compare two IUPAC characters for a mismatch
/// that is when they don't have any single bit in common
__device__ __forceinline__ u8 mismatch(u8 gchar, u8 tchar) {
  return ((gchar & tchar) == 0) ? 1 : 0;
}

/// This version uses shared memory only for the query and the strings
/// all tmp values are stored in registers
template <u32 warp_count, u32 qlen, u32 tlen>
__global__ void scores1(u8 *gquery, u8 *gstrings, u8 *result, u32 N) {

  cg::thread_block block = cg::this_thread_block();
  cg::thread_block_tile warp = cg::tiled_partition<warp_size>(block);
  auto warp_id = warp.meta_group_rank();
  auto lane = warp.thread_rank();

  // Cache stored in registers, each thread requires qlen x 2 bytes,
  // that is qlen x 2 / 4 registers of 32 bits each.
  u8 cache[2][qlen];
  u8 cache_curr_page = 0;

  // Keep strings and query in shared memory for fast access
  __shared__ u8 strings[warp_count][warp_size][tlen];
  __shared__ u8 query[qlen];

// Load strings, each warp load its strings
#pragma unroll(warp_size)
  for (auto sidx = 0; sidx < warp_size; ++sidx)
    if (lane < tlen) {
      const auto global_idx = (blockDim.x * blockIdx.x * tlen) +
                              (warp_id * warp_size * tlen) + (sidx * tlen);
      strings[warp_id][sidx][lane] = gstrings[global_idx + lane];
    }

  // Load query for all the warps
  if (warp_id == 0 && lane < qlen)
    query[lane] = gquery[lane];

// Initialize cache with base values
#pragma unroll(qlen)
  for (auto qidx = 0; qidx < qlen; ++qidx) {
    cache[0][qidx] = qidx + 1;
  }

  // SAFETY: Guard for strings and query in shared memory
  block.sync();

#if DEBUG
  block.sync();
  if (warp_id == 0 && lane == 0) {

    // Print strings
    printf("strings:\n");
    for (u32 i = 0; i < warp_size; ++i) {
      printf(" %2d ", i);
      for (u32 s = 0; s < tlen; ++s) {
        printf("%c", strings[warp_id][i][s]);
      }
      printf("\n");
    }

    // Print query
    printf("query:\n    ");
    for (u32 s = 0; s < qlen; ++s) {
      printf("%c", query[s]);
    }
    printf("\n");
  }
#endif

  // End free alignment
  u8 best_global_score = 255;
  bool invalid = false;

#pragma unroll(tlen)
  for (u32 tidx = 0; tidx < tlen; ++tidx) {

    // Advance cache circular buffer
    // NOTE: underflow is OK
    u8 cache_prev_page = cache_curr_page;
    cache_curr_page = (cache_curr_page + 1) % 2;

#pragma unroll(qlen)
    for (u32 qidx = 0; qidx < qlen; ++qidx) {

      // 0: left, 1: up, 2: diag
      u8 parents[3] = {0};

      // NOTE: Is the index calculation optimized?
      parents[0] = cache[cache_prev_page][qidx];
      parents[1] = 0; // No gap penalties at first row
      parents[2] = 0;

      // If we are at the first row
      if (qidx != 0) {
        parents[1] = cache[cache_curr_page][qidx - 1];
        parents[2] = cache[cache_prev_page][qidx - 1];
      }

#if DEBUG
      if (warp_id == 0 && lane == 20) {
        printf("loaded (left: %d, up: %d, diag: %d)\n", parents[0], parents[1],
               parents[2]);
      }
#endif

      // Update costs

      parents[0] += cost_gaps;
      parents[1] += cost_gaps;

      // TODO: change memory layout of strings, have 'lane' last
      char tchar = strings[warp_id][lane][tidx];
      parents[2] += mismatch(query[qidx], tchar);

      // If target contains Ns then it is invalid
      if (tchar == 'N') invalid = true;

#if DEBUG
      if (warp_id == 0 && lane == 20) {
        printf("update (left: %d, up: %d, diag: %d) for q:%c vs t:%c\n",
               parents[0], parents[1], parents[2], query[qidx],
               strings[warp_id][lane][tidx]);
      }
#endif

      // Find best cell
      u8 best_score = 255;

#pragma unroll(3)
      for (u32 i = 0; i < 3; ++i)
        if (parents[i] < best_score) {
          best_score = parents[i];
        }

#if DEBUG
      if (warp_id == 0 && lane == 20) {
        printf("store %d\n", best_score);
      }
#endif

      // Store results
      cache[cache_curr_page][qidx] = best_score;

      // Store best global result
      if (qidx == qlen - 1 && best_score < best_global_score)
        best_global_score = best_score;

      warp.sync();
    }

    warp.sync();
  }

  // Store the best global score
  if (invalid) best_global_score = 200; 
  result[blockDim.x * blockIdx.x + block.thread_rank()] = best_global_score;

#if DEBUG
  printf("score %2d: %d\n", lane, best_global_score);
#endif
}

void scores(const u8 *query, const u8 *strings, u8 *result, int qlen, int slen, int n) {

  constexpr u32 warp_size = 32;
  constexpr u32 warp_count = 4;

  constexpr u32 kernel_qlen = 20;
  constexpr u32 kernel_slen = 22;

  assert(kernel_qlen == qlen);
  assert(kernel_slen == slen);

  unsigned char *dev_query, *dev_strings, *dev_result;
  CUDA_CHECK(cudaMalloc(&dev_result, n));
  CUDA_CHECK(cudaMalloc(&dev_strings, (u64)slen * n));
  CUDA_CHECK(cudaMalloc(&dev_query, qlen));

  // Copy CPU memory to GPU
  CUDA_CHECK(cudaMemcpy(dev_strings, strings, sizeof(u8) * slen * n,
             cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(dev_query, query, sizeof(u8) * qlen, cudaMemcpyHostToDevice));

  u32 blocks = (n + warp_count * warp_size - 1) / (warp_count * warp_size);
  //printf("launching kernel with %d blocks\n", blocks);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  scores1<warp_count, kernel_qlen, kernel_slen>
      <<<blocks, warp_count * warp_size>>>(dev_query, dev_strings, dev_result, n);

  cudaError_t err = cudaGetLastError();                                
  if (err != cudaSuccess) {                                            
      fprintf(stderr, "CUDA kernel launch error in '%s' at line %d: %s (%d)\n",
                    __FILE__, __LINE__, cudaGetErrorString(err), err);
  } 

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  float time;
  cudaEventElapsedTime(&time, start, stop);
  //printf("elapsed kernel time: %.2f ms\n", time);

  // Get result for GPU
  CUDA_CHECK(cudaMemcpy(result, dev_result, n, cudaMemcpyDeviceToHost));

  CUDA_CHECK(cudaFree(dev_query));
  CUDA_CHECK(cudaFree(dev_strings));
  CUDA_CHECK(cudaFree(dev_result));
}

/// Allocate memory on the device
u8* cuda_malloc(u64 bytes) {
  void* result = nullptr;
  cudaMalloc(&result, bytes);
  return (u8*)result;
}

/// Free memory on the device
void cuda_free(u8* memory) {
  cudaFree((void*)memory);
}

void cuda_memcpy_to_gpu(u8* gpu, const u8* cpu, u64 bytes) {
  cudaMemcpy((void*)gpu, (const void*)cpu, bytes, cudaMemcpyHostToDevice);
}

void cuda_memcpy_to_cpu(const u8* gpu, u8* cpu, u64 bytes) {
  cudaMemcpy((void*)cpu, (const void*)gpu, bytes, cudaMemcpyDeviceToHost);

}
