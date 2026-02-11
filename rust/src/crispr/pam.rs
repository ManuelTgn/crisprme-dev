use crate::sequence::iupac::Iupac;


/// IUPAC bitmask for `N` (any base).
///
/// In this representation, `N` is `0b1111` (A|C|G|T). It is used both as a
/// degenerate PAM code and as an ambiguity marker in the input sequence.
/// In the scanner, k-mers containing `N` are skipped (see `scan_targets`).
const N_MASK: u8 = 0b1111;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedPAMError {
    pub position: usize,
    pub byte: u8,
}

impl std::fmt::Display for ParsedPAMError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Invalid PAM character '{}' (ASCII {}) at position {}",
            self.byte as char, self.byte, self.position
        )
    }    
}

impl std::error::Error for ParsedPAMError { } 

/// Parsed representation of a Protospacer Adjacent Motif (PAM).
///
/// This structure stores a PAM sequence encoded as IUPAC bitmasks, together
/// with its reverse complement and a flag indicating whether the PAM is fully
/// unconstrained (i.e., composed exclusively of `N` characters).
///
/// The bitmask representation enables extremely fast matching using bitwise
/// operations and supports both exact and degenerate PAM definitions.
pub struct ParsedPAM {
    /// PAM sequence encoded as IUPAC bitmasks (forward orientation).
    ///
    /// Each element is a 4-bit mask representing the set of allowed bases
    /// at that position.
    pub bytes: Vec<u8>,

    /// Reverse-complement of the PAM sequence, also encoded as IUPAC bitmasks.
    ///
    /// This is precomputed to allow strand-aware scanning without performing
    /// reverse-complement operations during the hot scanning loop.
    pub revcomp: Vec<u8>,

    /// Flag indicating whether the PAM is fully unconstrained.
    ///
    /// This is `true` if and only if all PAM positions are `N`
    /// (`0b1111` in IUPAC encoding), meaning the PAM imposes no constraints
    /// on matching.
    ///
    /// This flag enables fast-path optimizations during scanning, where PAM
    /// checks can be skipped entirely.
    pub unconstrained: bool,
}


impl ParsedPAM{
    pub fn new(pam: &str) -> Result<Self, ParsedPAMError> {
        // convert each ASCII nucleotide to its IUPAC bitmask.
        let mut bytes = Vec::with_capacity(pam.len()); 
        
        for (i, &b) in pam.as_bytes().iter().enumerate() {
            let code = Iupac::try_from_ascii(b)
                .ok_or(ParsedPAMError { position: i, byte: b })?
                .as_u8();
            bytes.push(code);
        }

        let unconstrained = bytes.iter().all(|&m| m == N_MASK);

        // Reverse-complement via Iupac::complement (keeps one source of thruth)
        let revcomp: Vec<u8> = bytes
            .iter()
            .rev()
            .map(|&m| Iupac::new(m).complement().as_u8())
            .collect();

        Ok(Self { bytes, revcomp, unconstrained })
    }
}


/// Computes the complement of an IUPAC nucleotide bitmask.
///
/// This function swaps the bit positions corresponding to complementary bases:
/// * A ↔ T (`0001` ↔ `1000`)
/// * C ↔ G (`0010` ↔ `0100`)
///
/// The operation correctly handles both standard and degenerate IUPAC codes
/// (e.g., `N` complements to `N`).
///
/// # Arguments
/// * `mask` - 4-bit IUPAC nucleotide bitmask
///
/// # Returns
/// The complemented IUPAC bitmask.
fn complement_bitmask(mask: u8) -> u8 {
    // extract each base's bit
    let a = mask & 0b0001;
    let c = mask & 0b0010;
    let g = mask & 0b0100;
    let t = mask & 0b1000;
    
    // perform the bit swaps for complementation:
    let complement_t = a << 3; 
    let complement_g = c << 1; 
    let complement_c = g >> 1; 
    let complement_a = t >> 3;
    
    // combine the resulting complement bits
    complement_t | complement_g | complement_c | complement_a
}

/// Builds a sparse representation of a PAM pattern by retaining only *informative* positions.
///
/// In IUPAC encoding, the mask `0b1111` (`N`) matches any nucleotide and therefore
/// does not constrain matching. This function filters out such positions and returns:
///   - the indices of PAM positions that impose constraints, and
///   - the corresponding IUPAC bitmasks.
///
/// This representation reduces per-k-mer matching work for partially-degenerate PAMs
/// (e.g., `NNGRRT`, `GGNRG`) by checking only informative positions.
///
/// # Arguments
/// * `pam` - Slice of IUPAC bitmasks representing the PAM sequence.
///
/// # Returns
/// A tuple `(idx, mask)` where:
/// * `idx[i]` is the position within the PAM of the `i`-th informative base.
/// * `mask[i]` is the corresponding IUPAC bitmask at that position.
///
/// # Notes
/// * If all PAM positions are unconstrained (`N`), both vectors will be empty.
/// * If no positions are unconstrained, `idx.len() == pam.len()`.
#[inline]
pub fn build_sparse(pam: &[u8]) -> (Vec<usize>, Vec<u8>) {
    // define vectors of indeexes and masks
    let mut idx: Vec<usize> = Vec::with_capacity(pam.len());
    let mut mask: Vec<u8> = Vec::with_capacity(pam.len());

    // iterate over pam nts
    for (i, &m) in pam.iter().enumerate() {
        if m != N_MASK {
            idx.push(i);
            mask.push(m);
        }
    }
    (idx, mask)
}

