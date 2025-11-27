use std::result::Result;  // explicitely import Result for clarity

/// Represents a nucleotide using the 4-bit IUPAC ambiguity code bitmask.
/// 
/// The standard mapping is:
/// - A: 0b0001 (1)
/// - C: 0b0010 (2)
/// - G: 0b0100 (4)
/// - T: 0b1000 (8)
/// 
/// The contained `u8` value is the bitmask.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Iupac(pub u8);  // bit mask wrapper

impl Iupac {
    /// Converts a single ASCII nucleotide character (case-insensitive) into its 
    /// corresponding 4-bit IUPAC bitmask, wrapped in `Iupac`.
    /// 
    /// # Arguments
    /// * `nt` - The ASCII byte representing the nucleotide (e.g., `b'A'`, `b'N'`).
    /// 
    /// # Returns
    /// * `Ok(Iupac)` on success, containing the bitmask.
    /// * `Err(String)` if the character is not a recognized IUPAC code.
    pub fn from_ascii(nt: u8) -> Result<Self, String> {
        // convert to uppercase before matching to ensure case-insensitivity
        let code = match nt.to_ascii_uppercase() {
            // standard bases
            b'A' => 0b0001,
            b'C' => 0b0010,
            b'G' => 0b0100,
            b'T' => 0b1000,

            // two-base ambiguities
            b'R' => 0b0101, // A or G
            b'Y' => 0b1010, // C or T
            b'S' => 0b0110, // G or C
            b'W' => 0b1001, // A or T
            b'K' => 0b1100, // G or T
            b'M' => 0b0011, // A or C
            
            // three-base ambiguities (Not T, Not A, Not G, Not C)
            b'B' => 0b1110, // C or G or T (Not A)
            b'D' => 0b1101, // A or G or T (Not C)
            b'H' => 0b1011, // A or C or T (Not G)
            b'V' => 0b0111, // A or C or G (Not T)

            // any base
            b'N' => 0b1111, // A, C, G, or T

            // handle unknown characters
            _ => { 
                return Err(format!(
                    "Unknown nucleotide: '{}' (ASCII: {}). Valid IUPAC codes are: A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N",
                    nt as char, nt
                ));
            }
        };

        Ok(Iupac(code))
    }
}

/// Converts an IUPAC bitmask back into its corresponding single-character 
/// representation.
/// 
/// This is the inverse of the mapping performed by `Iupac::from_ascii`.
/// 
/// # Arguments
/// * `bitmask` - The 4-bit IUPAC code (`u8`).
/// 
/// # Returns
/// * The single ASCII character representing the bitmask, or '?' for unknown codes
pub fn iupac_to_char(bitmask: u8) -> char {
    match bitmask {
        // standard bases
        0b0001 => 'A',
        0b0010 => 'C',
        0b0100 => 'G',
        0b1000 => 'T',

        // two-base ambiguities
        0b0101 => 'R', 
        0b1010 => 'Y', 
        0b0110 => 'S', 
        0b1001 => 'W', 
        0b1100 => 'K', 
        0b0011 => 'M', 

        // three-base ambiguities
        0b1110 => 'B', 
        0b1101 => 'D', 
        0b1011 => 'H', 
        0b0111 => 'V', 

        // any base
        0b1111 => 'N', 

        // safety fallback for non-IUPAC or zero masks
        _ => '?',
    }
}

/// Checks if a nucleotide bitmask matches a pattern bitmask.
/// 
/// In IUPAC coding, a match occurs if the set of possible bases in the nucleotide (`nt`) 
/// overlaps with the set of possible bases in the pattern (`pattern`).
/// 
/// This is achieved by checking if the bitwise AND operation yields a non-zero result.
/// 
/// # Arguments
/// * `nt` - The bitmask of the sequence nucleotide.
/// * `pattern` - The bitmask of the pattern/template nucleotide.
/// 
/// # Returns
/// * `true` if there is at least one common base (match), `false` otherwise
pub fn matches_iupac(nt: u8, pattern: u8) -> bool {
    (nt & pattern) != 0
}