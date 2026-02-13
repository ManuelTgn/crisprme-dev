#pragma once

#include <stdint.h>

using u64 = uint64_t;
using u32 = uint32_t;
using s32 = int32_t;
using u8  = uint8_t;

#define CUDA_CHECK(call)                                                      \
    do {                                                                      \
        cudaError_t err = call;                                               \
        if (err != cudaSuccess) {                                             \
            fprintf(stderr, "CUDA error in file '%s' at line %d: %s (%d)\n",  \
                    __FILE__, __LINE__, cudaGetErrorString(err), err);        \
            exit(EXIT_FAILURE);                                               \
        }                                                                     \
    } while (0)
