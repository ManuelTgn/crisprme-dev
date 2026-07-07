//! Target extraction: splitting a scan window into a canonical protospacer
//! and a PAM variant index.
//!
//! The scanner ([`crate::sequence::scanner`]) reports *where* candidate
//! targets are (a local position + a strand bit) but leaves the window bytes
//! untouched — every hit still carries the full `size`-base window in the
//! genome's **forward** orientation, PAM included. This module owns the next
//! step: given that forward window and the strand it was found on, produce
//!
//! 1. the **protospacer** (the window with the PAM removed), in canonical
//!    5'->3' orientation for that strand, and
//! 2. the **PAM variant index** ([`ParsedPAM::pam_index`]), i.e. which
//!    concrete member of the PAM's finite variant set this occurrence used.
//!
//! # Window geometry
//!
//! A window is `size` bases. The PAM occupies a contiguous sub-slice at one
//! end, determined (as in the scanner) by the `upstream` flag:
//!
//! ```text
//!  upstream == true            upstream == false
//!  ┌─────┬────────────┐      ┌────────────┬─────┐
//!  │ PAM │ protospacer│      │ protospacer│ PAM │
//!  └─────┴────────────┘      └────────────┴─────┘
//!   0   plen        size     0        size-plen size
//! ```
//!
//! `pam_start_fwd = if upstream { 0 } else { size - plen }` — identical to the
//! scanner's convention, so the two modules agree on where the PAM lives.
//!
//! # Strand canonicalization
//!
//! * **Forward strand** (`strand_bit == 1`): the window is already 5'->3';
//!   slice out the PAM and protospacer directly.
//! * **Reverse strand** (`strand_bit == 0`): the target reads 5'->3' on the
//!   *other* strand, so we reverse-complement the whole window first. Under
//!   reverse-complement the PAM region maps back onto the canonical forward
//!   position `pam_start_fwd` (`pam_start_rev` <-> `pam_start_fwd`), so after
//!   the flip the **same** slicing and the **same** forward enumeration table
//!   apply — no separate reverse index scheme is needed.
//!
//! The upshot is a strand-invariant key: a forward target and its
//! reverse-strand twin collapse to the *same* protospacer bytes and the
//! *same* PAM index, which is exactly what lets the batcher deduplicate
//! across strands. The occurrence's original genomic strand is recorded
//! separately by the batcher and is not lost.
//!
//! # Hot-path discipline
//!
//! This module sits on the per-hit path, so it carries **no tracing**:
//! per-occurrence logging would both dominate runtime and flood the logs.
//! The extraction geometry is logged once, at [`TargetExtractor::new`]
//! (called once per batcher). Per-hit invariants are checked with
//! `debug_assert!` and compiled out of release builds.

use crate::crispr::pam::ParsedPAM;
use crate::error::crisprme_errors::TargetError;
use crate::model::input::SEQ_MAX_LEN;
use crate::sequence::iupac::Iupac;

/// Precomputed window geometry for splitting hits into (protospacer, PAM).
///
/// Constructed once per [`crate::batching::batching::TargetBatcher`] and then
/// shared (it is [`Copy`]) across every `feed_chunk` call. It owns no data and
/// borrows nothing, so it is trivially `Send + Sync` and cheap to pass by
/// value into the hot loop.
///
/// All ranges are half-open `[lo, hi)` offsets into a forward window of length
/// [`Self::window_len`].
#[derive(Clone, Copy, Debug)]
pub struct TargetExtractor {
    /// Full window length in bases (`size`).
    size: usize,
    /// PAM length in bases (`plen`).
    plen: usize,
    /// PAM slice start (canonical forward orientation).
    pam_lo: usize,
    /// PAM slice end (`pam_lo + plen`).
    pam_hi: usize,
    /// Protospacer slice start.
    proto_lo: usize,
    /// Protospacer slice end.
    proto_hi: usize,
}

impl TargetExtractor {
    /// Build the extraction geometry for a `(plen, size, upstream)` configuration.
    ///
    /// # Arguments
    /// * `plen`  — PAM length in bases (must be `1..=size`).
    /// * `size`  — full scan-window width in bases (must be `<= SEQ_MAX_LEN`).
    /// * `upstream` — `true` if the PAM sits at the 5' end of the forward window
    ///             (offset `0`), `false` if it sits at the 3' end
    ///             (offset `size - plen`). This matches the scanner's flag.
    ///
    /// # Errors
    /// * [`TargetError::PamOutOfRange`] if `plen == 0` or `plen > size`.
    /// * [`TargetError::WindowTooLong`] if `size > SEQ_MAX_LEN` (the
    ///   reverse-complement scratch buffer is a fixed `SEQ_MAX_LEN` stack
    ///   array, and stored protospacers must fit a `SeqFrame` row).
    ///
    /// Emits one `DEBUG` trace describing the resolved geometry; because this
    /// runs once per batcher it is safe to log here without flooding.
    pub fn new(plen: usize, size: usize, upstream: bool) -> Result<Self, TargetError> {
        if plen == 0 || plen > size {
            return Err(TargetError::PamOutOfRange { plen, size });
        }
        if size > SEQ_MAX_LEN {
            return Err(TargetError::WindowTooLong { size, max: SEQ_MAX_LEN });
        }

        let pam_lo = if upstream { 0 } else { size - plen };
        let pam_hi = pam_lo + plen;
        let (proto_lo, proto_hi) = if upstream { (plen, size) } else { (0, size - plen) };

        tracing::debug!(
            "target geometry: size={size}, plen={plen}, placement={}, \
             pam=[{pam_lo}..{pam_hi}), protospacer=[{proto_lo}..{proto_hi}) \
             (proto_len={})",
            if upstream { "5'" } else { "3'" },
            proto_hi - proto_lo,
        );

        Ok(Self { size, plen, pam_lo, pam_hi, proto_lo, proto_hi })
    }

    /// Length of the extracted protospacer (`size - plen`).
    #[inline(always)]
    pub fn proto_len(&self) -> usize {
        self.proto_hi - self.proto_lo
    }

    /// Length of the full scan window (`size`).
    #[inline(always)]
    pub fn window_len(&self) -> usize {
        self.size
    }

    /// PAM length (`plen`).
    #[inline(always)]
    pub fn pam_len(&self) -> usize {
        self.plen
    }

    /// Split one hit into its canonical protospacer and PAM variant index.
    ///
    /// Writes the protospacer bytes (IUPAC bitmasks, canonical 5'->3' for the
    /// hit's strand) into `out`, clearing it first, and returns the PAM
    /// variant index for the occurrence.
    ///
    /// `out` is a caller-owned scratch buffer so the batcher can reuse a
    /// single allocation across all hits in a chunk and only pay for a boxed
    /// map key when a genuinely new unique window is inserted.
    ///
    /// # Arguments
    /// * `pam`         — the parsed PAM (source of the variant enumeration).
    /// * `window`      — the forward-orientation window, length == `size`.
    /// * `strand_bit`  — `1` = forward (+), `0` = reverse (−).
    /// * `out`         — protospacer output buffer (cleared and refilled).
    ///
    /// # Invariants (from the scanner, `debug_assert`ed)
    /// * `window.len() == size`.
    /// * every window base is pure (no `N`), so each PAM base is a single bit.
    /// * for a reverse hit, the reverse-complemented PAM region is consistent
    ///   with `pam.bytes`, which [`ParsedPAM::pam_index`] relies on.
    #[inline]
    pub fn extract(
        &self,
        pam: &ParsedPAM,
        window: &[u8],
        strand_bit: u8,
        out: &mut Vec<u8>,
    ) -> u16 {
        debug_assert_eq!(window.len(), self.size, "window length must equal size");
        out.clear();

        if strand_bit == 1 {
            // Forward strand: the window is already canonical 5'->3'.
            out.extend_from_slice(&window[self.proto_lo..self.proto_hi]);
            pam.pam_index(&window[self.pam_lo..self.pam_hi])
        } else {
            // Reverse strand: canonicalize by reverse-complementing the window
            // into a fixed stack buffer (window fits by construction). After
            // the flip the PAM sits at the canonical forward offset, so the
            // same slicing and forward enumeration apply.
            let mut rc = [0u8; SEQ_MAX_LEN];
            let rc = &mut rc[..self.size];
            revcomp_into(window, rc);

            out.extend_from_slice(&rc[self.proto_lo..self.proto_hi]);
            pam.pam_index(&rc[self.pam_lo..self.pam_hi])
        }
    }
}

/// Reverse-complement `src` into `dst` (both IUPAC bitmask slices of equal
/// length), using [`Iupac::complement`] as the single source of truth for
/// base complementation (degenerate codes included).
#[inline]
fn revcomp_into(src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    let n = src.len();
    for i in 0..n {
        dst[i] = Iupac::new(src[n - 1 - i]).complement().as_u8();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Concrete single-base IUPAC masks.
    const A: u8 = 0b0001;
    const C: u8 = 0b0010;
    const G: u8 = 0b0100;
    const T: u8 = 0b1000;

    fn revcomp_vec(src: &[u8]) -> Vec<u8> {
        let mut out = vec![0u8; src.len()];
        revcomp_into(src, &mut out);
        out
    }

    /// A forward target and its reverse-strand genomic image must yield the
    /// SAME protospacer and the SAME PAM index. This is the core
    /// strand-canonicalization guarantee.
    #[test]
    fn strand_roundtrip_pam_5_prime() {
        let pam = ParsedPAM::new("NGG").unwrap();
        // upstream = true  -> PAM at the 5' end of the window.
        let ex = TargetExtractor::new(3, 5, true).unwrap();

        // Canonical target 5'->3': PAM = CGG, protospacer = AT.
        let target = [C, G, G, A, T];

        let mut fwd = Vec::new();
        let idx_fwd = ex.extract(&pam, &target, 1, &mut fwd);
        assert_eq!(idx_fwd, 1, "CGG is variant 1 of NGG");
        assert_eq!(fwd, vec![A, T]);

        // The same target seen on the reverse strand is stored (forward
        // orientation) as its reverse complement.
        let genomic = revcomp_vec(&target);
        let mut rev = Vec::new();
        let idx_rev = ex.extract(&pam, &genomic, 0, &mut rev);

        assert_eq!(idx_rev, idx_fwd, "PAM index must be strand-invariant");
        assert_eq!(rev, fwd, "protospacer must be strand-invariant");
    }

    #[test]
    fn strand_roundtrip_pam_3_prime() {
        let pam = ParsedPAM::new("NGG").unwrap();
        // upstream = false -> PAM at the 3' end of the window.
        let ex = TargetExtractor::new(3, 5, false).unwrap();

        // Canonical target 5'->3': protospacer = AT, PAM = CGG.
        let target = [A, T, C, G, G];

        let mut fwd = Vec::new();
        let idx_fwd = ex.extract(&pam, &target, 1, &mut fwd);
        assert_eq!(idx_fwd, 1);
        assert_eq!(fwd, vec![A, T]);

        let genomic = revcomp_vec(&target);
        let mut rev = Vec::new();
        let idx_rev = ex.extract(&pam, &genomic, 0, &mut rev);
        assert_eq!(idx_rev, idx_fwd);
        assert_eq!(rev, fwd);
    }

    /// Unconstrained PAMs (all-N, up to the parse-time cap) still record the
    /// index — it simply encodes the observed concrete PAM bases.
    #[test]
    fn unconstrained_pam_records_observed_bases() {
        let pam = ParsedPAM::new("NNN").unwrap();
        let ex = TargetExtractor::new(3, 5, true).unwrap();

        // PAM = CGT (variant of NNN), protospacer = AA.
        let target = [C, G, T, A, A];
        let mut proto = Vec::new();
        let idx = ex.extract(&pam, &target, 1, &mut proto);

        // Decoding the index must recover the exact PAM bases.
        assert_eq!(pam.pam_variant(idx).unwrap(), vec![C, G, T]);
        assert_eq!(proto, vec![A, A]);
    }

    #[test]
    fn construction_rejects_bad_geometry() {
        assert!(matches!(
            TargetExtractor::new(0, 5, true),
            Err(TargetError::PamOutOfRange { plen: 0, size: 5 })
        ));
        assert!(matches!(
            TargetExtractor::new(6, 5, true),
            Err(TargetError::PamOutOfRange { plen: 6, size: 5 })
        ));
        assert!(matches!(
            TargetExtractor::new(3, SEQ_MAX_LEN + 4, true),
            Err(TargetError::WindowTooLong { .. })
        ));
    }
}