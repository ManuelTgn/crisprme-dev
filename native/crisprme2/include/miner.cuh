#pragma once

#include <common.cuh>

#ifdef CRISPRME_IMPLEMENTATION

// Follows same implementation of rust `Cigarx64`
struct Cigarx64
{
    // Number of bits for the type `S`
    static constexpr u32 BITS = sizeof(u64) * 8;
    // Number of operations upported
    static constexpr u32 CAPACITY = (BITS - 1) / 2;

    // All supported operations
    enum class Oper : u8
    {
        M = 0b00, // `=` (sequence match) 
        X = 0b01, // `X` (sequence mismatch)
        D = 0b10, // `D` (deletion from reference)
        I = 0b11  // `I` (insertion to reference)
    };

    // We always set the sentinel bit at first position
    u64 m_storage = 1;

    // Push a new operation
    __device__ __forceinline__ void push(Oper value)
    {
#if DEBUG
        assert(len() < CAPACITY);
#endif
        m_storage = (m_storage << 2) | static_cast<u64>(value);
    }

    // Pop the last operation
    __device__ __forceinline__ void pop()
    {
#if DEBUG
        assert(len() > 0);
#endif
        m_storage = m_storage >> 2;
    }

    // Number of operations
    __device__ __forceinline__ u32 len()
    {
        u32 bits = BITS - __clzll(m_storage) - 1;
        return bits / 2;
    }

    // Returns the inner value
    __device__ __forceinline__ u64 storage()
    {
        return m_storage;
    }
};

__device__ const Cigarx64::Oper CIGARX_ENCODE_64[8] = {
    Cigarx64::Oper::X, // 0b0'00 invalid!
    Cigarx64::Oper::D, // 0b0'01 OK
    Cigarx64::Oper::I, // 0b0'10 OK
    Cigarx64::Oper::M, // 0b0'11 OK
    Cigarx64::Oper::X, // 0b1'00 invalid!
    Cigarx64::Oper::D, // 0b1'01 OK
    Cigarx64::Oper::I, // 0b1'10 OK
    Cigarx64::Oper::X, // 0b1'11 OK
};

/// Single exploration step in the stack machine
struct Step
{
    static constexpr u8 MASK = static_cast<u8>(0b11);

    /// All possible step types
    enum Inner : u8
    {
        B = 0b11, // Increment index of guide AND sequence (=/X)
        S = 0b10, // Increment index of sequence (I)
        G = 0b01, // Increment index of guide (D)
        E = 0b00  // Backtrack

    } value;

    /// Create step with initial type
    static __device__ __forceinline__ Step initial()
    {
        return Step{Inner::B};
    }

    /// Create step as a deletion (last possible type)
    static __device__ __forceinline__ Step deletion()
    {
        return Step{Inner::G};
    }

    /// Advance the step type
    __device__ __forceinline__ Step next()
    {
        return Step{(Inner)(value - 1)};
    }

    /// Returns true if we are a backtrack
    __device__ __forceinline__ bool is_backtrack()
    {
        return value == Inner::E;
    }

    /// Returns true if this step has a bit for match/mismatch
    __device__ __forceinline__ bool is_match_or_mismatch()
    {
        return value == Inner::B;
    }

    /// Returns the guide index delta caused by this step
    __device__ __forceinline__ u32 gidx_dt()
    {
        return (value & Inner::G) != 0;
    }

    /// Returns the sequence index delta caused by this step
    __device__ __forceinline__ u32 sidx_dt()
    {
        return (value & Inner::S) != 0;
    }

    __device__ __forceinline__ u32 ggap_dt()
    {
        return value == Inner::S;
    }

    __device__ __forceinline__ u32 sgap_dt()
    {
        return value == Inner::G;
    }

    /// Convert the step to a CIGARX value
    __device__ __forceinline__ Cigarx64::Oper to_cigarx_op_u64(bool mismatch)
    {
        u32 index = (mismatch << 2) | value;
        return CIGARX_ENCODE_64[index];
    }

    __device__ __forceinline__ char as_char()
    {
        static const char *letters = "EGSB";
        return letters[value];
    }
};

/// Hold the current running state of a thread
/// NOTE: Implemented as a single u32 to fit inside a single register
struct ThreadState
{
    u32 gidx : 8;
    u32 sidx : 8;
    u32 mism : 8;
    u32 ggap : 4;
    u32 sgap : 4;
};

/// Thread miner that keeps track of steps and match/mismatches
template <typename Storage>
struct ThreadMiner
{
    /// Persistent data across launches
    struct Inner
    {
        /// Contains a bitpacked version of the current explorations steps
        /// NOTE: This is the same inside each thread of the warp
        Storage stack;

        /// Bitset for mismatches, on bit for each step present in the stack (1 = mismatch, 0 = match/other)
        /// NOTE: This differs inside each thread of the warp
        u32 mismatches;

        /// Current running state
        /// NOTE: This differs inside each thread of the warp
        ThreadState state;

        /// Number of steps
        /// NOTE: This is the same inside each thread of the warp
        u8 len;

    } mem;

    __device__ __forceinline__ ThreadMiner()
    {
        mem.stack = 0;
        mem.mismatches = 0;
        mem.state = {0, 0, 0, 0, 0};
        mem.len = 0;
    }

    __device__ __forceinline__ bool has_work()
    {
        return mem.len != 0;
    }

    __device__ __forceinline__ void decrement_state(Step s)
    {
        mem.state.gidx -= s.gidx_dt();
        mem.state.sidx -= s.sidx_dt();
        mem.state.ggap -= s.ggap_dt();
        mem.state.sgap -= s.sgap_dt();
    }

    __device__ __forceinline__ void increment_state(Step s)
    {
        mem.state.gidx += s.gidx_dt();
        mem.state.sidx += s.sidx_dt();
        mem.state.ggap += s.ggap_dt();
        mem.state.sgap += s.sgap_dt();
    }

    /// Replace the current step without changing the mismatch bit
    __device__ __forceinline__ void replace(Step s)
    {
        mem.stack = (mem.stack & ~0b11) | (s.value & 0b11);
    }

    /// Push a new step on top of the stack
    __device__ __forceinline__ void push(Step s, bool set_mismatch)
    {
        assert(mem.len <= sizeof(mem.stack) * 8 / 2 && "StepStack too small!");

        // Add state to memory
        mem.stack = (mem.stack << 2) | s.value;
        mem.mismatches = (mem.mismatches << 1) | set_mismatch;
        mem.len += 1;

        // Update running state
        mem.state.mism = __popc(mem.mismatches);
        increment_state(s);
    }

    /// Pop the top of the stack and update running state
    __device__ __forceinline__ Step pop()
    {
        Step s = current();

        // Remove state from memory
        mem.mismatches = mem.mismatches >> 1;
        mem.stack = mem.stack >> 2;
        mem.len -= 1;

        // Update running state
        mem.state.mism = __popc(mem.mismatches);
        decrement_state(s);

        return s;
    }

    /// Travel the current step
    __device__ __forceinline__ void travel()
    {
        Step step = current();
        decrement_state(step);

        // Set mismatch bit at zero after a travel!
        // Only the initial B step uses the bit, all other step types do not use it.
        mem.mismatches = mem.mismatches & ~0b1;
        mem.state.mism = __popc(mem.mismatches);

        step = step.next();
        increment_state(step);
        replace(step);
    }

    /// Get the current step
    __device__ __forceinline__ Step current()
    {
        return Step{(Step::Inner)(mem.stack & 0b11)};
    }

    /// Check if we are inside the thresholds
    __device__ __forceinline__ bool inside_thresholds(u32 max_mism, u32 max_ggap, u32 max_sgap)
    {
        return (
            mem.state.mism <= max_mism &&
            mem.state.ggap <= max_ggap &&
            mem.state.sgap <= max_sgap);
    }

    /// Check if we can continue with the exploration, that is if the guide is not completelly matched
    /// and the target sequences index is not outside the relaxed bounds
    __device__ __forceinline__ bool can_continue(u32 slen, u32 glen, u32 offset, u32 max_ggap)
    {
        return (
            mem.state.gidx < glen &&
            mem.state.sidx + offset < slen + max_ggap);
    }

    /// Check if the current state is final, the guide has been completelly matched
    __device__ __forceinline__ bool is_complete(u32 glen)
    {
        return mem.state.gidx == glen;
    }

    /// Get the current CIGARX
    __device__ __forceinline__ Cigarx64 cigarx()
    {
        Cigarx64 result{1};
        for (int i = mem.len - 1; i >= 0; --i)
        {
            Step step{(Step::Inner)((mem.stack >> i * 2) & 0b11)};
            bool mismatch = (mem.mismatches >> i & 0b1);
            result.push(step.to_cigarx_op_u64(mismatch));
        }
        return result;
    }

    /// Print the current stack
    __device__ __forceinline__ void print()
    {
        for (int i = mem.len - 1; i >= 0; --i)
        {
            Step step{(Step::Inner)((mem.stack >> i * 2) & 0b11)};
            printf("%c", step.as_char());
        }
        printf("\n");
    }

    /// Print inner state
    __device__ __forceinline__ void print_state()
    {
        printf("next_sidx: %d\n", mem.state.sidx);
        printf("next_gidx: %d\n", mem.state.gidx);
        printf("sgap: %d\n", mem.state.sgap);
        printf("ggap: %d\n", mem.state.ggap);
        printf("mism: %d\n", mem.state.mism);
    }
};

#endif