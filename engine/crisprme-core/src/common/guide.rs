use super::iupac::Iupac;

/// A guide for mining
#[derive(Clone, PartialEq, Eq)]
pub struct Guide(Vec<Iupac>);

impl Guide {
    /// Create the reverse complement
    pub fn reverse_complement(&self) -> Self {
        Self(self.0.iter().cloned()
            .map(Iupac::complement).rev()
            .collect())
    } 

    /// String representation
    pub fn as_string(&self) -> String {
        self.0.iter()
            .map(|e| Iupac::to_utf8(*e))
            .collect()
    }
}

/// Create from a string
impl From<&str> for Guide {
    fn from(value: &str) -> Self {
        Self(value.chars()
            .map(Iupac::from_utf8)
            .collect())
    }
}

/// Indexin inside a guide
impl std::ops::Index<usize> for Guide {
    type Output = Iupac;
    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

impl std::fmt::Debug for Guide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for c in self.0.iter() {
            write!(f, "{}", c.to_utf8())?;
        }
        Ok(())
    }
}

impl std::fmt::Display for Guide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for c in self.0.iter() {
            write!(f, "{}", c.to_utf8())?;
        }
        Ok(())
    }
}

/// Transparent over the inner vector
impl std::ops::Deref for Guide {
    type Target = Vec<Iupac>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn reverse_complement_normal() {
        
        const GUIDE_SRC: &str = "ATTGAGATAGTGGNGG";
        const GUIDE_REV: &str = "CCNCCACTATCTCAAT";

        let guide: Guide = GUIDE_SRC.into();
        let guide = guide.reverse_complement();

        for (a, b) in guide.iter().zip(GUIDE_REV.chars()) {
            assert_eq!(a.to_utf8(), b);
        }
    }

    #[test]
    fn reverse_complement_iupac() {
 
        const GUIDE_SRC: &str = "RSH";
        const GUIDE_REV: &str = "DSY";

        let guide: Guide = GUIDE_SRC.into();
        let guide = guide.reverse_complement();

        let guide: String = guide.iter().cloned().map(Iupac::to_utf8).collect();
        assert_eq!(guide, GUIDE_REV, "output = {guide}, correct = {GUIDE_REV}");
    }
}
