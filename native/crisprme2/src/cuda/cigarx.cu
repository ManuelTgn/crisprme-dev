#include <core/src/cuda/iupac.cu>

template<typename Storage>
struct Cigarx {

    enum class Inner : u8 {
        M = 0b00,
        X = 0b01,
        D = 0b10,
        I = 0b11
    };

    Storage m_storage;
    u8 m_bits;

    __device__ __forceinline__ void push(Inner op) {
        m_storage = (m_storage << 2) | static_cast<u8>(op);
        m_bits += 2;
    }

    __device__ __forceinline__ static char as_str(Inner op) {
        switch (op) {
            case Inner::M: return '=';
            case Inner::X: return 'X';
            case Inner::D: return 'D';
            case Inner::I: return 'I';
            default: 
                return '?';
        }
    }

    /// Print the content of the cigarx
    __device__ __forceinline__ void print() {
        for(int i = (m_bits / 2) - 1; i >= 0; --i) {
            Inner op = static_cast<Inner>((m_storage >> (i * 2)) & 0b11);
            printf("%c", as_str(op));
        }
        printf("\n");
    }

    __device__ __forceinline__ void extract(char* result) {
        u32 write = 0;
        for(int i = (m_bits / 2) - 1; i >= 0; --i) {
            Inner op = static_cast<Inner>((m_storage >> (i * 2)) & 0b11);
            result[write] = as_str(op);
            write += 1;
        }
        result[write] = '\0';
    }
};

__device__ const char* STEP_CIGARX = "EIDXEID=";
__device__ const Cigarx<u64>::Inner CIGARX_ENCODE_64[8] = {
    Cigarx<u64>::Inner::X, // invalid!
    Cigarx<u64>::Inner::I,
    Cigarx<u64>::Inner::D,
    Cigarx<u64>::Inner::X,
    Cigarx<u64>::Inner::X, // invalid!
    Cigarx<u64>::Inner::I,
    Cigarx<u64>::Inner::D,
    Cigarx<u64>::Inner::M,
};

#define CIGARXOP_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define CIGARXOP_TO_BINARY(byte)  \
  ((byte) & 0x80 ? '1' : '0'), \
  ((byte) & 0x40 ? '1' : '0'), \
  ((byte) & 0x20 ? '1' : '0'), \
  ((byte) & 0x10 ? '1' : '0'), \
  ((byte) & 0x08 ? '1' : '0'), \
  ((byte) & 0x04 ? '1' : '0'), \
  ((byte) & 0x02 ? '1' : '0'), \
  ((byte) & 0x01 ? '1' : '0') 

#if 0

struct CigarxState {
  u8 valid, sidx, gidx, pad;
};

enum class CigarxOp : u8 {
    Match     = 0b00, // =
    Mismatch  = 0b01, // X
    Deletion  = 0b10, // D
    Insertion = 0b11  // I
};

__device__ __forceinline__ char cigarxop_decode(CigarxOp op) {
  switch (op) {
    case CigarxOp::Match:     return '=';
    case CigarxOp::Mismatch:  return 'X';
    case CigarxOp::Deletion:  return 'D';
    case CigarxOp::Insertion: return 'I';
  }

  printf("invalid cigarxop value: %u\n", static_cast<u8>(op));
  return '?';
}

template<typename Storage>
struct Cigarx {

    Storage storage; // bit-packed storage
    u8 bits; // current number of bits used

    __device__ __forceinline__ Cigarx() 
      : storage(0), bits(0) {}

    __device__ __forceinline__ bool empty() {
      return bits == 0;
    }

    __device__ __forceinline__ void push(CigarxOp op) {
        //printf("push %c\n", cigarxop_decode(op));
        //assert(count * 2 + 2 <= sizeof(T) * 8); // prevent overflow
        storage = (storage << 2) | (static_cast<u8>(op) & 0b11);
        bits += 2;

        //printf("result, ");
        //print();
    }

    __device__ __forceinline__ bool last_is_diagonal() {
      return (peek() & 0b10) == 0;
    }

    // Pop last operation
    __device__ __forceinline__ CigarxOp pop() {
        //assert(count > 0);

        CigarxOp val = static_cast<CigarxOp>(storage & 0b11);
        storage = storage >> 2;
        bits -= 2;

        //printf("pop %c\n", cigarxop_decode(val));
        return val;
    }

    __device__ __forceinline__ void pop_discard() {
      storage = storage >> 2;
      bits -= 2;
    }

    // Peek last operation
    __device__ __forceinline__ CigarxOp peek() {
      return storage & 0b11;
    }

    __device__ void print_solution(u32 bseq, u32 offset) {
      char buffer[20];

      u32 write = 0;
      for (int i = bits - 2; i >= 0; i -= 2) {
        u8 op = (storage >> i) & 0b11;
        buffer[write] = cigarxop_decode(static_cast<CigarxOp>(op));
        write += 1;
      }

      buffer[write] = '\0';
      printf("Solution(seq: %u, offset: %u, cigarx: %s)\n", bseq, offset, buffer);
    }

    /// Check if solution
    __device__ __forceinline__ bool valid(s32 sgap, s32 ggap, s32 mism) {
      #pragma unroll
      for(int i = 0; i < bits; i += 2) {
        CigarxOp op = (storage >> i) & 0b11;

        if (op == CigarxOp::Mismatch)  mism -= 1;
        if (op == CigarxOp::Deletion)  sgap -= 1;
        if (op == CigarxOp::Insertion) ggap -= 1;
      }

      return sgap >= 0 && ggap >= 0 && mism >= 0; 
    }

    /// Get the current sidx, gidx and validity state in a single function
    __device__ __forceinline__ CigarxState state(s32 sgap, s32 ggap, s32 mism) {
      CigarxState result = { 0 };
      for (int i = 0; i < bits; i += 2) {
        

      }
      return result;
    }

    /*
    /// Travel
    __device__ __forceinline__ bool travel() {
      CigarxOp last_op = pop();
      last_op += 1;

      // Skip mismatch, we are interested only on general movements
      if (last_op == CigarxOp::Mismatch) 
        last_op += 1;

      bool exhausted = (last_op == 4);
      if (!exhausted) push(last_op);
      return exhausted;
    }
    */
};

/// Up to 16 operations
using Cigarx32 = Cigarx<u32>;
/// Up to 32 operations
using Cigarx64 = Cigarx<u64>;

#endif
