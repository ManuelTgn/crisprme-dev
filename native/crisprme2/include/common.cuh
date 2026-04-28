#pragma once

#include <stdint.h>
#include <stdio.h>

using u64 = uint64_t;
using u32 = uint32_t;
using s32 = int32_t;
using u8  = uint8_t;

#define CUDA_CHECK(call)                                                     \
    do                                                                       \
    {                                                                        \
        cudaError_t err = call;                                              \
        if (err != cudaSuccess)                                              \
        {                                                                    \
            fprintf(stderr, "CUDA error in file '%s' at line %d: %s (%d)\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err), err);       \
            exit(EXIT_FAILURE);                                              \
        }                                                                    \
    } while (0)

/// ====================================================================
/// IUPAC
/// ====================================================================

// A IUPAC is represented using a single byte
using iupac = u8;

#ifdef CRISPRME_IMPLEMENTATION

__device__ __forceinline__ char iupac_decode(iupac v)
{
    switch (v)
    {
    case 0b0001:
        return 'A';
    case 0b0010:
        return 'C';
    case 0b0100:
        return 'G';
    case 0b1000:
        return 'T';
    case 0b0101:
        return 'R';
    case 0b1010:
        return 'Y';
    case 0b0110:
        return 'S';
    case 0b1001:
        return 'W';
    case 0b1100:
        return 'K';
    case 0b0011:
        return 'M';
    case 0b1110:
        return 'B';
    case 0b1101:
        return 'D';
    case 0b1011:
        return 'H';
    case 0b0111:
        return 'V';
    case 0b1111:
        return 'N';
    }
    printf("invalid iupac encoding!");
    return '?';
}

// Check if two IUPAC characters match
__device__ __forceinline__ bool iupac_match(iupac a, iupac b) {
  return (a & b) != 0;
}

#endif