use crate::iupac::{Iupac, matches_iupac, iupac_to_char};

/// Represents a Protospacer Adjacent Motif (PAM) sequence parsed into IUPAC bitmasks.
///
/// This struct holds both the forward PAM sequence and its reverse complement
/// as efficient vectors of 4-bit IUPAC masks.
pub struct ParsedPAM {
    // the PAM sequence converted into IUPAC bitmasks
    pub bytes: Vec<u8>,
    // the reverse complement of the PAM sequence, also in IUPAC bitmasks
    pub revcomp: Vec<u8>,
}

impl ParsedPAM{
    /// Creates a new `ParsedPAM` instance from an ASCII PAM string.
    ///
    /// The function converts the input string into forward bitmasks and calculates
    /// the reverse complement bitmasks.
    ///
    /// # Arguments
    /// * `pam` - The ASCII string representation of the PAM (e.g., "NGG").
    ///
    /// # Returns
    /// * `Ok(Self)` if parsing is successful.
    /// * `Err(String)` if any character in `pam` is not a valid IUPAC code.
    pub fn new(pam: &str) -> Result<Self, String> {
        // 1. Convert each ASCII nucleotide to its IUPAC bitmask.
        let bytes: Result<Vec<u8>, String> = pam.as_bytes()
            .iter()
            .map(|&b| Iupac::from_ascii(b).map(|iupac| iupac.0))
            .collect();

        // Use the '?' operator to extract the Vec<u8> or return the error String immediately
        let bytes: Vec<u8> = bytes?;

        // 2. Create reverse complement using bit masks: reverse order and complement each mask
        let revcomp: Vec<u8> = bytes.iter()
            .rev()
            .map(|&b| complement_bitmask(b))
            .collect();

        Ok(Self { bytes, revcomp })
    }

    /// Checks if a sequence fragment's bitmasks match the PAM pattern bitmasks.
    ///
    /// Matching uses the `matches_iupac` logic, checking for overlap in base possibilities.
    ///
    /// # Arguments
    /// * `seq` - A slice of `u8` IUPAC bitmasks representing the sequence fragment.
    ///
    /// # Returns
    /// * `true` if the sequence slice has the same length as the PAM and matches
    ///   the pattern at all positions; `false` otherwise.
    pub fn matches(&self, seq: &[u8]) -> bool {
        if seq.len() != self.bytes.len() {
            return false;
        }

        seq.iter()
            .zip(self.bytes.iter())
            .all(|(&nt_mask, &pattern_mask)| {
                matches_iupac(nt_mask, pattern_mask)
            })
    }
    
    /// Decodes the stored IUPAC bitmask vector back into a readable ASCII string.
    ///
    /// # Arguments
    /// * `rc` - If `true`, decodes the reverse complement (`self.revcomp`); 
    ///          otherwise, decodes the forward sequence (`self.bytes`).
    ///
    /// # Returns
    /// A `String` containing the standard IUPAC nucleotide characters (e.g., "NGG")
    pub fn decode(&self, rc: bool) -> String {
        if rc {
            self.revcomp.iter()
                // map each bitmask to its corresponding char
                .map(|&bitmask| iupac_to_char(bitmask))
                .collect()
        } else {
            self.bytes.iter()
                // map each bitmask to its corresponding char
                .map(|&bitmask| iupac_to_char(bitmask))
                .collect()
        }
    }
}

/// Complements a nucleotide IUPAC bitmask.
///
/// The function swaps the bit positions corresponding to:
/// A <-> T (0001 <-> 1000) and C <-> G (0010 <-> 0100).
/// This correctly handles both standard and degenerate IUPAC codes (e.g., 'N' complements to 'N').
///
/// # Arguments
/// * `mask` - The 4-bit IUPAC code (`u8`).
///
/// # Returns
/// The complemented IUPAC bitmask (`u8`).
fn complement_bitmask(mask: u8) -> u8 {
    // extract each base's bit
    let a = mask & 0b0001;
    let c = mask & 0b0010;
    let g = mask & 0b0100;
    let t = mask & 0b1000;
    
    // perform the bit swaps for complementation:
    // 1. A (0001) shifts left 3 to become T (1000)
    let complement_t = a << 3; 
    
    // 2. C (0010) shifts left 1 to become G (0100)
    let complement_g = c << 1; 
    
    // 3. G (0100) shifts right 1 to become C (0010)
    let complement_c = g >> 1; 
    
    // 4. T (1000) shifts right 3 to become A (0001)
    let complement_a = t >> 3;
    
    // combine the resulting complement bits
    complement_t | complement_g | complement_c | complement_a
}