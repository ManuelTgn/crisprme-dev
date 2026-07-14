//! Protospacer Adjacent Motif (PAM) parsing, matching, and variant indexing.
//!
//! A PAM is a short, possibly-degenerate IUPAC motif (e.g. SpCas9 `NGG`).
//! This module provides [`PAM`], which precomputes everything the
//! scanner and batcher need to work with a PAM at high throughput:
//!
//! * the forward IUPAC bitmasks ([`PAM::bytes`]),
//! * the reverse-complement bitmasks ([`PAM::revcomp`]), so the
//!   scanner never reverse-complements inside the hot loop,
//! * a fast-path flag for fully-unconstrained PAMs, and
//! * a **finite variant enumeration**: because each degenerate position
//!   admits a fixed set of concrete bases, the set of concrete PAMs a
//!   degenerate motif matches is finite and can be addressed by a single
//!   integer index.
//!
//! # PAM variant indexing
//!
//! For `NGG`, the concrete variants are `AGG, CGG, GGG, TGG` — four of
//! them, indexed `0..4`. In general the count is the product over
//! positions of the number of bases each position allows
//! (`popcount(mask)`). The index is a **mixed-radix** number with
//! **position 0 as the most significant digit**:
//!
//! ```text
//! index = ((rank_0 · r_1 + rank_1) · r_2 + rank_2) · …
//! ```
//!
//! where `r_i = popcount(bytes[i])` is the radix at position `i` and
//! `rank_i` is the position of the concrete base within the ascending
//! (A < C < G < T) list of bases allowed at `i`. This makes `AGG → 0`,
//! `CGG → 1`, `GGG → 2`, `TGG → 3`, and the mapping is a bijection onto
//! `0..variant_count`, so it round-trips ([`PAM::pam_index`] ∘
//! [`PAM::pam_variant`] is the identity).
//!
//! The index is represented as a `u16` downstream; [`PAM::new`]
//! rejects PAMs whose variant count would overflow that (see
//! [`PamError::TooManyVariants`]), which no real PAM approaches.

use crate::error::crisprme_errors::PamError;
use crate::sequence::iupac::Iupac;

/// IUPAC bitmask for `N` (any base).
///
/// In this representation, `N` is `0b1111` (A|C|G|T). It is used both as a
/// degenerate PAM code and as an ambiguity marker in the input sequence.
/// In the scanner, k-mers containing `N` are skipped (see `scan_targets`).
const N_MASK: u8 = 0b1111;

/// Parsed representation of a Protospacer Adjacent Motif (PAM).
///
/// Stores the PAM as IUPAC bitmasks together with its reverse complement,
/// a fast-path flag, and the metadata needed to encode/decode the PAM's
/// finite set of concrete variants (see the [module docs](self)).
///
/// The bitmask representation enables extremely fast matching using
/// bitwise operations and supports both exact and degenerate PAMs.
pub struct PAM {
    /// PAM sequence encoded as IUPAC bitmasks (forward orientation).
    ///
    /// Each element is a 4-bit mask representing the set of allowed bases
    /// at that position.
    pub bytes: Vec<u8>,

    /// Reverse-complement of the PAM sequence, also encoded as IUPAC
    /// bitmasks. Precomputed to allow strand-aware scanning without
    /// reverse-complementing inside the hot loop.
    pub revcomp: Vec<u8>,

    /// `true` iff every PAM position is `N` (`0b1111`), i.e. the PAM
    /// imposes no constraint. Enables skipping PAM checks entirely.
    pub unconstrained: bool,

    /// PAM length in bases (`== bytes.len()`). Cached for convenience.
    plen: usize,

    /// Number of concrete variants this (degenerate) PAM matches, i.e.
    /// the product of `popcount(mask)` over all positions. Guaranteed to
    /// fit in a `u16` index (`<= u16::MAX + 1`).
    variant_count: u32,

    /// The upper-cased IUPAC motif exactly as supplied (e.g. `"NGG"`, `"TTTV"`).
    /// Single source of truth for anything that renders the PAM as text.
    motif: Box<str>,
}

impl PAM {
    /// Parse a PAM string into its bitmask representation and precompute
    /// its reverse complement and variant enumeration.
    ///
    /// Emits a `debug`/`info` trace on success and an `error` trace on
    /// failure (routed to `verbose.log`/`basic.log`/`errors.log`
    /// respectively via the logging bridge).
    ///
    /// # Errors
    /// * [`PamError::InvalidCharacter`] — a byte is not a valid IUPAC code.
    /// * [`PamError::TooManyVariants`] — the PAM is so degenerate its
    ///   variant count would not fit a `u16` index (unreachable for real
    ///   PAMs).
    pub fn new(pam: &str) -> Result<Self, PamError> {
        // Normalise once, then parse the normalised form: `Iupac::try_from_ascii`
        // is case-sensitive, so this also makes lowercase input legal.
        let motif = pam.to_ascii_uppercase();

        // Convert each ASCII nucleotide to its IUPAC bitmask.
        let mut bytes = Vec::with_capacity(motif.len());
        for (i, &b) in motif.as_bytes().iter().enumerate() {
            match Iupac::try_from_ascii(b) {
                Some(code) => bytes.push(code.as_u8()),
                None => {
                    let err = PamError::InvalidCharacter {
                        position: i,
                        byte: b,
                    };
                    tracing::error!("{err}");
                    return Err(err);
                }
            }
        }

        let plen = bytes.len();
        let unconstrained = bytes.iter().all(|&m| m == N_MASK);

        // Reverse-complement via Iupac::complement (single source of truth).
        let revcomp: Vec<u8> = bytes
            .iter()
            .rev()
            .map(|&m| Iupac::new(m).complement().as_u8())
            .collect();

        // Variant count = product of per-position radices (popcounts).
        // Accumulate in u64 with saturation and bail early if it exceeds
        // the u16-index ceiling, so the downstream `u16` is provably safe.
        const MAX_VARIANTS: u64 = u16::MAX as u64 + 1; // 65_536
        let mut variant_count: u64 = 1;
        for &m in &bytes {
            variant_count = variant_count.saturating_mul(m.count_ones() as u64);
            if variant_count > MAX_VARIANTS {
                let err = PamError::TooManyVariants {
                    count: variant_count,
                    plen,
                    max: MAX_VARIANTS as u32,
                };
                tracing::error!("{err}");
                return Err(err);
            }
        }
        let variant_count = variant_count as u32;

        tracing::debug!(
            "parsed PAM {pam:?}: plen={plen}, variant_count={variant_count}, \
             unconstrained={unconstrained}"
        );
        tracing::info!("PAM {pam:?} ready ({variant_count} concrete variant(s))");

        Ok(Self {
            bytes,
            revcomp,
            unconstrained,
            plen,
            variant_count,
            motif: motif.into_boxed_str(),
        })
    }

    /// The degenerate IUPAC motif as ASCII, e.g. `"NGG"` or `"TTTV"`.
    #[inline]
    pub fn motif(&self) -> &str {
        &self.motif
    }

    /// PAM length in bases.
    #[inline(always)]
    pub fn plen(&self) -> usize {
        self.plen
    }

    /// Number of concrete variants this PAM matches (the addressable index
    /// range is `0..variant_count`).
    #[inline(always)]
    pub fn variant_count(&self) -> u32 {
        self.variant_count
    }

    /// Encode a concrete PAM occurrence into its variant index.
    ///
    /// `concrete` must be `plen` single-base IUPAC bitmasks (one bit set
    /// each) that are consistent with this PAM. In the pipeline this
    /// invariant is guaranteed: the scanner only accepts windows whose
    /// bases are pure (no `N`) and which already matched the PAM, so this
    /// call is infallible on the hot path and returns the index directly.
    ///
    /// The invariants are checked with `debug_assert!` (compiled out in
    /// release) rather than returning a `Result`, to keep the hot path
    /// branch-free of error handling. See the [module docs](self) for the
    /// mixed-radix scheme.
    #[inline]
    pub fn pam_index(&self, concrete: &[u8]) -> u16 {
        debug_assert_eq!(
            concrete.len(),
            self.plen,
            "concrete PAM length must equal plen"
        );

        let mut index: u32 = 0;
        for i in 0..self.plen {
            let mask = self.bytes[i];
            let base = concrete[i];

            debug_assert!(
                base.count_ones() == 1,
                "concrete PAM base must be a single pure base"
            );
            debug_assert!(
                mask & base != 0,
                "concrete base is not allowed by the PAM at this position"
            );

            let radix = mask.count_ones(); // number of allowed bases here
            let rank = (mask & (base - 1)).count_ones(); // allowed bases below `base`
            index = index * radix + rank; // position 0 == most significant
        }

        debug_assert!(index <= u16::MAX as u32);
        index as u16
    }

    /// Decode a variant index back into its concrete PAM bitmasks.
    ///
    /// Inverse of [`pam_index`](Self::pam_index). Intended for cold paths
    /// (output/annotation), so it is fallible and validates the range.
    ///
    /// # Errors
    /// [`PamError::IndexOutOfRange`] if `index >= variant_count`.
    pub fn pam_variant(&self, index: u16) -> Result<Vec<u8>, PamError> {
        if index as u32 >= self.variant_count {
            return Err(PamError::IndexOutOfRange {
                index,
                count: self.variant_count,
            });
        }

        let mut rem = index as u32;
        let mut out = vec![0u8; self.plen];
        // Unwind the mixed-radix number: position 0 is most significant,
        // so recover digits from the least significant (last) position up.
        for i in (0..self.plen).rev() {
            let mask = self.bytes[i];
            let radix = mask.count_ones();
            let rank = rem % radix;
            rem /= radix;
            out[i] = nth_set_base(mask, rank);
        }
        Ok(out)
    }

    /// Decode a variant index directly to an ASCII PAM string (e.g. `"CGG"`).
    ///
    /// Convenience wrapper over [`pam_variant`](Self::pam_variant) for
    /// human-readable output.
    ///
    /// # Errors
    /// [`PamError::IndexOutOfRange`] if `index >= variant_count`.
    pub fn pam_variant_ascii(&self, index: u16) -> Result<String, PamError> {
        Ok(self
            .pam_variant(index)?
            .into_iter()
            .map(|mask| Iupac::new(mask).to_ascii() as char)
            .collect())
    }
}

/// Return the `n`-th set bit of `mask` (0-indexed, ascending bit order) as
/// a single-bit mask. Used to turn a per-position rank back into a base.
///
/// Requires `n < popcount(mask)`; returns `0` otherwise.
#[inline]
fn nth_set_base(mask: u8, n: u32) -> u8 {
    let mut remaining = mask;
    let mut count = 0u32;
    while remaining != 0 {
        let lowest = remaining & remaining.wrapping_neg(); // isolate lowest set bit
        if count == n {
            return lowest;
        }
        remaining &= remaining - 1; // clear lowest set bit
        count += 1;
    }
    0
}

/// Computes the complement of an IUPAC nucleotide bitmask.
///
/// Swaps the bit positions of complementary bases (A↔T, C↔G), handling
/// degenerate codes correctly (e.g. `N` complements to `N`).
fn complement_bitmask(mask: u8) -> u8 {
    let a = mask & 0b0001;
    let c = mask & 0b0010;
    let g = mask & 0b0100;
    let t = mask & 0b1000;

    let complement_t = a << 3;
    let complement_g = c << 1;
    let complement_c = g >> 1;
    let complement_a = t >> 3;

    complement_t | complement_g | complement_c | complement_a
}

/// Builds a sparse representation of a PAM pattern, keeping only the
/// *informative* (non-`N`) positions.
///
/// Returns `(idx, mask)` where `idx[i]` is the position of the `i`-th
/// informative base and `mask[i]` is its IUPAC bitmask. If every position
/// is `N` both vectors are empty; if none is `N`, `idx.len() == pam.len()`.
#[inline]
pub fn build_sparse(pam: &[u8]) -> (Vec<usize>, Vec<u8>) {
    let mut idx: Vec<usize> = Vec::with_capacity(pam.len());
    let mut mask: Vec<u8> = Vec::with_capacity(pam.len());
    for (i, &m) in pam.iter().enumerate() {
        if m != N_MASK {
            idx.push(i);
            mask.push(m);
        }
    }
    (idx, mask)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Concrete single-base masks, for building test PAM occurrences.
    const A: u8 = 0b0001;
    const C: u8 = 0b0010;
    const G: u8 = 0b0100;
    const T: u8 = 0b1000;

    #[test]
    fn ngg_variant_enumeration() {
        let pam = PAM::new("NGG").unwrap();
        assert_eq!(pam.variant_count(), 4);
        assert_eq!(pam.pam_index(&[A, G, G]), 0);
        assert_eq!(pam.pam_index(&[C, G, G]), 1);
        assert_eq!(pam.pam_index(&[G, G, G]), 2);
        assert_eq!(pam.pam_index(&[T, G, G]), 3);
    }

    #[test]
    fn two_degenerate_positions_roundtrip() {
        // R = {A,G}, Y = {C,T}, G -> 2 * 2 * 1 = 4 variants.
        let pam = PAM::new("RYG").unwrap();
        assert_eq!(pam.variant_count(), 4);
        // Round-trip every index: decode -> encode must be the identity.
        for idx in 0..pam.variant_count() as u16 {
            let variant = pam.pam_variant(idx).unwrap();
            assert_eq!(pam.pam_index(&variant), idx, "roundtrip failed at {idx}");
        }
        // Position 0 is most significant: idx 0 == A.C == smallest bases.
        assert_eq!(pam.pam_variant_ascii(0).unwrap(), "ACG");
        assert_eq!(pam.pam_variant_ascii(3).unwrap(), "GTG");
    }

    #[test]
    fn decode_rejects_out_of_range() {
        let pam = PAM::new("NGG").unwrap();
        assert!(matches!(
            pam.pam_variant(4),
            Err(PamError::IndexOutOfRange { index: 4, count: 4 })
        ));
    }

    #[test]
    fn invalid_character_is_reported_with_position() {
        assert!(matches!(
            PAM::new("NXG"),
            Err(PamError::InvalidCharacter {
                position: 1,
                byte: b'X'
            })
        ));
    }
}
