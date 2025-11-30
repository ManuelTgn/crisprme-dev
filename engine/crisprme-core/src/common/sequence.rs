use std::borrow::Cow;

use super::iupac::Iupac;

/// A view into a sequece of Iupac characters
#[derive(Clone, PartialEq, Eq)]
pub struct Sequence<'b>(Cow<'b, [Iupac]>);

impl<'b> Sequence<'b> {
    /// Create new sequence from a slice
    pub fn new(bytes: &'b [Iupac]) -> Self {
        Self(Cow::Borrowed(bytes))
    }

    pub fn from_utf8(seq: &str) -> Self {
        Self(Cow::Owned(
            seq.chars().map(Iupac::from_utf8).collect(),
        ))
    }

    /// Get the inner slice
    pub fn as_slice(&self) -> &[Iupac] {
        self.0.as_ref()
    }

    /// Get mutation score of the sequence
    pub fn mutation_score(&self) -> u32 {
        self.0
            .iter()
            .map(|&e| Iupac::mutation_score(e))
            .sum::<u32>()
    }

    /// Get string representation
    pub fn as_string(&self) -> String {
        self.as_slice().into_iter()
            .map(|e| Iupac::to_utf8(*e))
            .collect()
    }
}

// =================================================================================
// STD implementations

impl<'s> std::ops::Deref for Sequence<'s> {
    type Target = [Iupac];
    fn deref(&self) -> &Self::Target {
        self.0.as_ref()
    }
}

/// Indexing inside a sequence
impl<'s> std::ops::Index<usize> for Sequence<'s> {
    type Output = Iupac;
    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

/// Pretty print as a DNA string
impl<'s> std::fmt::Display for Sequence<'s> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for c in self.0.as_ref() {
            write!(f, "{}", c.to_utf8())?;
        }
        Ok(())
    }
}

/// Debug show as DNA string instead of Vec<Iupac>
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
        assert_eq!(sequence.mutation_score(), 9);
    }
}
