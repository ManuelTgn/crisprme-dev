use std::borrow::Cow;

use super::iupac::Iupac;

/// A lightweight view (borrowed or owned) over a sequence of IUPAC-encoded bases.
///
/// `Sequence<'b>` is used as an ergonomic adapter type throughout the alignment
/// and scanning code: it can wrap either
///
/// - a borrowed slice (`&'b [Iupac]`) **without allocations**, or
/// - an owned buffer (`Vec<Iupac>`) created on demand (e.g., from a `&str`).
///
/// This pattern is especially useful in a high-throughput aligner:
/// - hot paths can pass references to pre-encoded buffers cheaply
/// - user-facing/debug/test code can still construct sequences from strings
///
/// Internally, the representation is:
/// ```text
/// Cow<'b, [Iupac]>
/// ```
/// which can be either `Borrowed(&'b [Iupac])` or `Owned(Vec<Iupac>)`.
///
/// # Encoding assumption
/// `Iupac` is expected to represent a 4-bit ambiguity mask. Many methods (e.g.
/// mutation score) interpret the mask as "number of possible bases".
#[derive(Clone, PartialEq, Eq)]
pub struct Sequence<'b>(Cow<'b, [Iupac]>);

impl<'b> Sequence<'b> {
    /// Construct a borrowed sequence view from a slice.
    ///
    /// This is the preferred constructor for hot paths: it performs no allocations
    /// and simply wraps the provided slice.
    pub fn new(bytes: &'b [Iupac]) -> Self {
        Self(Cow::Borrowed(bytes))
    }

    /// Construct an owned sequence from ASCII IUPAC text, **lossy**.
    ///
    /// Invalid characters are mapped to `N` (wildcard) to avoid panics.
    ///
    /// Use this if you want robustness and don't care about rejecting invalid input.
    ///
    /// Note: A Rust `&str` is UTF-8. We interpret it as ASCII bytes here, which is
    /// what IUPAC codes are.
    #[inline]
    pub fn from_ascii_lossy(seq: &str) -> Self {
        let mut v = Vec::with_capacity(seq.len());
        for &b in seq.as_bytes() {
            v.push(Iupac::from_ascii_lossy(b));
        }
        Self(Cow::Owned(v))
    }

    /// Construct an owned sequence from ASCII IUPAC text, **strict**.
    ///
    /// Returns `None` if any character is not a valid IUPAC code.
    ///
    /// Use this in validation-heavy code paths.
    #[inline]
    pub fn try_from_ascii(seq: &str) -> Option<Self> {
        let mut v = Vec::with_capacity(seq.len());
        for &b in seq.as_bytes() {
            v.push(Iupac::try_from_ascii(b)?);
        }
        Some(Self(Cow::Owned(v)))
    }

    /// Constructor equivalent to [`from_ascii_lossy`]. If you want strict behavior,
    /// use [`try_from_ascii`].
    #[inline]
    pub fn from_utf8(seq: &str) -> Self {
        Self::from_ascii_lossy(seq)
    }

    /// Return the inner slice view.
    ///
    /// This never allocates: if the sequence is owned, it returns a slice into
    /// the owned buffer.
    #[inline]
    pub fn as_slice(&self) -> &[Iupac] {
        self.0.as_ref()
    }

    /// Length of the sequence in bases
    #[inline]
    pub fn len(&self) -> usize {
        self.0.as_ref().len()
    }

    /// Returns 'true' if the sequence is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate over the sequence
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, Iupac> {
        self.as_slice().iter()
    }

    /// Compute the "mutation score" of the sequence.
    ///
    /// This score is defined as the sum over positions of the number of possible
    /// bases represented by each IUPAC mask:
    /// - `A/C/G/T` contribute 1
    /// - `R` (A|G) contributes 2
    /// - `B` (C|G|T) contributes 3
    /// - `N` (A|C|G|T) contributes 4
    ///
    /// In other words, if your `Iupac` mask is a 4-bit set, this is:
    /// `sum(popcount(mask))`.
    ///
    /// This is useful as a quick measure of how "degenerate" a sequence is
    /// (higher = more ambiguity).
    #[inline]
    pub fn mutation_score(&self) -> u32 {
        self.0
            .iter()
            .map(|&e| Iupac::mutation_score(e))
            .sum::<u32>()
    }

    /// Render as a human-readable IUPAC string.
    ///
    /// Intended for logging/debugging and tests.
    /// Avoid calling in tight loops as it allocates a `String`.
    #[inline]
    pub fn as_string(&self) -> String {
        let mut s = String::with_capacity(self.len());
        for &e in self.as_slice() {
            s.push(e.to_utf8())
        }
        s
    }

    /// Convert into an owned `Vec<Iupac>` (clones if borrowed).
    #[inline]
    pub fn to_owned_vec(&self) -> Vec<Iupac> {
        self.0.to_vec()
    }
}

// =============================================================================
// STD implementations

/// Treat `Sequence` as a slice of `Iupac` for convenient read-only access.
impl<'s> std::ops::Deref for Sequence<'s> {
    type Target = [Iupac];

    #[inline]
    fn deref(&self) -> &Self::Target {
        self.0.as_ref()
    }
}

/// Indexing inside a sequence.
/// Panics on out-of-bounds indexing (standard slice semantics).
impl<'s> std::ops::Index<usize> for Sequence<'s> {
    type Output = Iupac;

    #[inline]
    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

/// Pretty-print as DNA string
impl<'s> std::fmt::Display for Sequence<'s> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for c in self.0.as_ref() {
            write!(f, "{}", c.to_utf8())?;
        }
        Ok(())
    }
}

/// Debug output includes mutation score and the IUPAC string.
impl<'s> std::fmt::Debug for Sequence<'s> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Sequence(mutation: {}, iupac: \"", self.mutation_score())?;
        for c in self.0.as_ref() {
            write!(f, "{}", c.to_utf8())?;
        }
        write!(f, "\")")
    }
}

impl<'s> std::hash::Hash for Sequence<'s> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // delegate hashing to the slice of Iupac
        self.0.hash(state);
    }
}

impl<'s> Ord for Sequence<'s> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.as_ref().cmp(other.0.as_ref())
    }
}

impl<'s> PartialOrd for Sequence<'s> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<'s> From<Sequence<'s>> for Vec<Iupac> {
    fn from(seq: Sequence<'s>) -> Self {
        seq.0.to_vec()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn mutation_score() {
        let elements: &[Iupac; 4] = &[
            Iupac::from_utf8('A'),
            Iupac::from_utf8('N'),
            Iupac::from_utf8('A'),
            Iupac::from_utf8('B'),
        ];

        let sequence = Sequence::new(elements);
        // A(1) + N(4) + A(1) + B(3) = 9
        assert_eq!(sequence.mutation_score(), 9);
    }

    #[test]
    fn strict_rejects_invalid() {
        assert!(Sequence::try_from_ascii("ACZ").is_none());
    }

    #[test]
    fn lossy_maps_invalid_to_n() {
        let s = Sequence::from_ascii_lossy("ACZ").as_string();
        assert_eq!(s, "ACN");
    }
}
