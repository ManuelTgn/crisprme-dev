#include <common.cuh>

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

/*
struct Step2 {

    static constexpr u8 MASK      = static_cast<u8>(0b111);
    static constexpr u8 MASK_OPER = static_cast<u8>(0b011);
    static constexpr u8 MASK_MISM = static_cast<u8>(0b100);

    // Operation bitfield
    enum class Oper : u8 {
	M = 0b11, // Match or mismatch, both guide and sequence moved
	G = 0b10, // Insertion, only guide moved
	S = 0b01, // Deletion, only sequence moved
	E = 0b00, // Exhausted, no movement
    } m_oper;

    /// Is this a mismatch?
    bool m_mism;

    /// Create an initial step as a deletion
    __device__ __forceinline__ static Step initial() {
	Element elem;
	elem.oper = Oper::S & MASK_OPER;
	elem.mism = 0;
        return Step { elem }; // Deletion
    }

    /// Create new step from bits
    __device__ __forceinline__ Step(Element e) {
        m_elem = e;
    }

    __device__ __forceinline__ u8 bits() {
	u8 bits = m_elem.oper | (m_elem.mism << 2);
        return bits & MASK;
    }

    __device__ __forceinline__ u8 oper() {
        return m_elem.oper;
    }
    
    /// Create next step from this
    __device__ __forceinline__ Step next(u8* g, u8* s) {
	u8 next_oper = oper() + 1 & MASK_OPER;
        u32 mmism = bits() & 0b100;
        u32 o = oper() - 1;
        return Step { mmism | o };
    }
   
    __device__ __forceinline__ bool is_backtrack() {
        return oper() == 0;
    }

    __device__ __forceinline__ u32 sidx_dt() {
        return (bits() & 0b01) != 0;
    }

    __device__ __forceinline__ u32 gidx_dt() {
        return (bits() & 0b10) != 0;
    }

    __device__ __forceinline__ u32 mism_dt() {
        return m_value == Inner::XB;
    }

    __device__ __forceinline__ u32 ggap_dt() {
        return oper() == 0b01;
    }

    __device__ __forceinline__ u32 sgap_dt() {
        return oper() == 0b10;
    }

    __device__ __forceinline__ Cigarx<u64>::Inner to_cigarx_op() {
        return CIGARX_ENCODE_64[bits()];
    }

    __device__ __forceinline__ char as_str() {
      return STEP_CIGARX[bits()];
    }
}
*/

struct Step {

    /// Get only the required bits of the step
    static constexpr u8 MASK = static_cast<u8>(0b111);

    /// Exploration step, 1 bit for match and 2 bits for oper
    enum class Inner : u8 {
        MB = 0b111, // Both guide and sequence (match)
        XB = 0b011, // Both guide and sequence (mismatch)
        MG = 0b110, // Guide (match)
        XG = 0b010, // Guide (mismatch)
        MS = 0b101, // Sequence (match)
        XS = 0b001, // Sequence (mismatch)
        ME = 0b100, // Exhausted (match)
        XE = 0b000, // Exhausted (mismath)

    } m_value;

    /// Create an initial step as a match/mismatch
    __device__ __forceinline__ static Step initial(bool match) {
        u32 bits = (match << 2) | 0b011;
        return Step { bits };
    }

    /// Create a deletion step
    __device__ __forceinline__ static Step deletion() {
	return Step { 0b001 };
    }

    /// Create new step from bits
    __device__ __forceinline__ Step(u32 bits) {
        m_value = static_cast<Inner>(bits & MASK);
    }

    __device__ __forceinline__ u8 bits() {
        return static_cast<u8>(m_value);
    }

    __device__ __forceinline__ u8 oper() {
        return bits() & 0b11;
    }
    
    /// Create next step from this
    __device__ __forceinline__ Step next() {
        u32 mmism = bits() & 0b100;
        u32 o = oper() - 1;
        return Step { mmism | o };
    }
   
    __device__ __forceinline__ bool is_backtrack() {
        return oper() == 0;
    }

    __device__ __forceinline__ u32 sidx_dt() {
        return (bits() & 0b01) != 0;
    }

    __device__ __forceinline__ u32 gidx_dt() {
        return (bits() & 0b10) != 0;
    }

    __device__ __forceinline__ u32 mism_dt() {
        return m_value == Inner::XB;
    }

    __device__ __forceinline__ u32 ggap_dt() {
        return oper() == 0b01;
    }

    __device__ __forceinline__ u32 sgap_dt() {
        return oper() == 0b10;
    }

    __device__ __forceinline__ Cigarx<u64>::Inner to_cigarx_op() {
        return CIGARX_ENCODE_64[bits()];
    }

    __device__ __forceinline__ char as_str() {
      return STEP_CIGARX[bits()];
    }
};

typedef struct {
    unsigned long long lo;
    unsigned long long hi;
} uint128_t;

__device__ inline uint128_t shl128(uint128_t x, unsigned int shift) {
    uint128_t r;

    if (shift > 64) {
	r.lo = x.lo << (shift - 64);
	r.hi = 0ULL;
	return r;
    }

    r.hi = (x.hi << shift) | (x.lo >> (64 - shift));
    r.lo = x.lo << shift;
    return r;
}

__device__ inline uint128_t shr128(uint128_t x, unsigned int shift) {
    uint128_t r;

    if (shift >= 64) {
	r.hi = x.hi << (shift - 64);
	r.lo = 0;
    }

    r.lo = (x.lo >> shift) | (x.hi << (64 - shift));
    r.hi = x.hi >> shift;
    return r;
}

struct StepStack128 {

    uint128_t m_storage;
    u8 m_len;

    /// Create an empty stack
    __device__ __forceinline__ StepStack128() {
	m_storage.hi = 0ULL;
	m_storage.lo = 0ULL;
	m_len = 0;
    }

    __device__ __forceinline__ u8 len() {
        return m_len;
    }

    /// Reserve space for a new Step
    __device__ __forceinline__ void reserve() {
        m_storage = shl128(m_storage, 3);
        m_len += 1;
    }

    /// Replace top of stack
    __device__ __forceinline__ void replace(Step s) {
        m_storage.lo = (m_storage.lo & ~0b111ULL) | s.bits(); 
    }

    /// Push a new step
    __device__ __forceinline__ void push(Step s) {
        reserve(); replace(s);
    }

    /// Pop the current element
    __device__ __forceinline__ void pop() {
        m_storage = shr128(m_storage, 3);
        m_len -= 1;
    }

    /// Get top of stack
    __device__ __forceinline__ Step current() {
        u8 bits = static_cast<u8>(m_storage.lo & Step::MASK);
        return Step { bits };
    }

    /// Print the current stack
    __device__ __forceinline__ void print() {
      uint128_t storage = m_storage;
      for(int i = m_len - 1; i >= 0; --i) {
        Step step { static_cast<u32>(shr128(m_storage, i * 3).lo) };
        printf("%c", step.as_str());
      }
      printf("\n");
    }

    /// Get the compacted cigarx
    __device__ __forceinline__ Cigarx<u64> cigarx() {
	Cigarx<u64> result = { 0 };
        for(int i = m_len - 1; i >= 0; --i) {
            Step step { static_cast<u32>(shr128(m_storage, i * 3).lo) };
            result.push(step.to_cigarx_op());
        }
        return result;
    }
};

template<typename T>
struct StepStack {

    /// Maximum number of steps that can be stored
    static constexpr int MAX_LEN = (sizeof(T) * 8) / 3;

    T m_storage;
    u8 m_len;

    /// Create an empty stack
    __device__ __forceinline__ StepStack()
        : m_storage{0}, m_len{0} { }

    __device__ __forceinline__ u8 len() {
        return m_len;
    }

    /// Reserve space for a new Step
    __device__ __forceinline__ void reserve() {
        m_storage = m_storage << 3;
        m_len += 1;

	//assert(m_len < MAX_LEN);
    }

    /// Replace top of stack
    __device__ __forceinline__ void replace(Step s) {
        m_storage = (m_storage & ~0b111) | s.bits(); 
    }

    /// Push a new step
    __device__ __forceinline__ void push(Step s) {
        reserve(); replace(s);
    }

    /// Pop the current element
    __device__ __forceinline__ void pop() {
        m_storage = m_storage >> 3;
        m_len -= 1;
    }

    /// Get top of stack
    __device__ __forceinline__ Step current() {
        u8 bits = static_cast<u8>(m_storage & Step::MASK);
        return Step { bits };
    }

    /// Print the current stack
    __device__ __forceinline__ void print() {
      for(int i = m_len - 1; i >= 0; --i) {
        Step step { m_storage >> static_cast<T>(i * 3) };
        printf("%c", step.as_str());
      }
      printf("\n");
    }

    /// Get the compacted cigarx
    __device__ __forceinline__ Cigarx<u64> cigarx() {
        Cigarx<u64> result = { 0 };
        for(int i = m_len - 1; i >= 0; --i) {
            Step step { m_storage >> static_cast<T>(i * 3) };
            result.push(step.to_cigarx_op());
        }
        return result;
    }

    /// Get current gidx
    __device__ __forceinline__ u32 gidx() {
        u32 result = 0;
        for (int i = 0; i < MAX_LEN; ++i) {
            if (i < m_len) {
                Step step { m_storage >> static_cast<T>(i * 3) };
                result += step.gidx_dt();
            }
        }
        return result;
    }

    /// Get current sidx
    __device__ __forceinline__ u32 sidx() {
        u32 result = 0;
        for (int i = 0; i < MAX_LEN; ++i) {
            if (i < m_len) {
                Step step { m_storage >> static_cast<T>(i * 3) };
                result += step.sidx_dt();
            }
        }
        return result;
    }

    /// Get current mism
    __device__ __forceinline__ u32 mism() {
        u32 result = 0;
        for (int i = 0; i < MAX_LEN; ++i) {
            if (i < m_len) {
                Step step { m_storage >> static_cast<T>(i * 3) };
                result += step.mism_dt();
            }
        }
        return result;
    }

    /// Get current ggap
    __device__ __forceinline__ u32 ggap() {
        u32 result = 0;
        for (int i = 0; i < MAX_LEN; ++i) {
            if (i < m_len) {
                Step step { m_storage >> static_cast<T>(i * 3) };
                result += step.ggap_dt();
            }
        }
        return result;
    }

    /// Get current sgap
    __device__ __forceinline__ u32 sgap() {
        u32 result = 0;
        for (int i = 0; i < MAX_LEN; ++i) {
            if (i < m_len) {
                Step step { m_storage >> static_cast<T>(i * 3) };
                result += step.sgap_dt();
            }
        }
        return result;
    }
};
