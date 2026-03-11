

static constexpr __device__ u32 MASK(u32 bits) { return (1u << bits) - 1; }

/// Alignment state packed inside a single u32
/// [sidx:5][qidx:5][sgap:4][qgap:4][mism:4][oper:2][reserved:8]
struct AlignmentState {
    u32 storage; 

    // --- Getters ---  
    __device__ __forceinline__ u32 get_sidx() const { return (storage >> 27) & MASK(5); }
    __device__ __forceinline__ u32 get_gidx() const { return (storage >> 22) & MASK(5); }
    __device__ __forceinline__ u32 get_sgap() const { return (storage >> 18) & MASK(4); }
    __device__ __forceinline__ u32 get_ggap() const { return (storage >> 14) & MASK(4); }
    __device__ __forceinline__ u32 get_mism() const { return (storage >> 10) & MASK(4); }
    __device__ __forceinline__ u32 get_oper() const { return (storage >>  8) & MASK(2); }

    // --- Setters ---
    __device__ __forceinline__ void set_sidx(u32 val) { storage = (storage & ~(MASK(5) << 27)) | ((val & MASK(5)) << 27); }
    __device__ __forceinline__ void set_gidx(u32 val) { storage = (storage & ~(MASK(5) << 22)) | ((val & MASK(5)) << 22); }
    __device__ __forceinline__ void set_sgap(u32 val) { storage = (storage & ~(MASK(4) << 18)) | ((val & MASK(4)) << 18); }
    __device__ __forceinline__ void set_ggap(u32 val) { storage = (storage & ~(MASK(4) << 14)) | ((val & MASK(4)) << 14); }
    __device__ __forceinline__ void set_mism(u32 val) { storage = (storage & ~(MASK(4) << 10)) | ((val & MASK(4)) << 10); }
    __device__ __forceinline__ void set_oper(u32 val) { storage = (storage & ~(MASK(2) <<  8)) | ((val & MASK(2)) <<  8); }

    // --- Utility ---
    __device__ __forceinline__ AlignmentState match() {
        AlignmentState result = AlignmentState { storage };
        result.set_gidx(get_gidx() + 1);
        result.set_sidx(get_sidx() + 1);
        result.set_oper(0);
        return result;
    }

    __device__ __forceinline__ AlignmentState mismatch() {
        AlignmentState result = AlignmentState { storage };
        result.set_gidx(get_gidx() + 1);
        result.set_sidx(get_sidx() + 1);
        result.set_mism(get_mism() + 1);
        result.set_oper(0);
        return result;
    }
  
    __device__ __forceinline__ AlignmentState insertion() {
        AlignmentState result = AlignmentState { storage };
        result.set_gidx(get_gidx() + 1);
        result.set_ggap(get_ggap() + 1);
        result.set_oper(0);
        return result;
    }

    __device__ __forceinline__ AlignmentState deletion() {
        AlignmentState result = AlignmentState { storage };
        result.set_sidx(get_sidx() + 1);
        result.set_sgap(get_sgap() + 1);
        result.set_oper(0);
        return result;
    }

    /// Check is the state is invalid compared to the thresholds
    __device__ __forceinline__ bool invalid(u32 ggap, u32 sgap, u32 mism) {
        return get_mism() > mism || get_sgap() > sgap || get_ggap() > sgap;
    }

    __device__ __forceinline__ AlignmentState travel() {
        AlignmentState result = AlignmentState { storage };
        result.set_oper(get_oper() + 1);
        return result;
    }

    __device__ void print() {
        printf("AlignmentState(gidx: %u, sidx: %u, ggap: %u, sgap: %u, mism: %u, oper: %u)\n",
               get_gidx(), get_sidx(), get_ggap(), get_sgap(), get_mism(), get_oper());
    }
};

