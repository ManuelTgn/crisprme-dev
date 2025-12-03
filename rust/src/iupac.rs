use std::result::Result;  // explicitely import Result for clarity

const IUPAC_LOOKUP_TABLE: [u8; 256] = {
    // 0b1111 is the 'N' mask (Any base). We use 0b0000 (0) as the default for 
    // invalid/unrecognized characters to mark them distinctly.
    // The lookup function will handle the 0b0000 case for error reporting.
    let mut table = [0u8; 256];

    // Helper macro/function to set both cases for case-insensitivity
    macro_rules! set_iupac {
        ($char:expr, $mask:expr) => {
            table[$char as usize] = $mask;
            table[$char.to_ascii_lowercase() as usize] = $mask;
        };
    }
    
    // Standard Bases
    set_iupac!(b'A', 0b0001);
    set_iupac!(b'C', 0b0010);
    set_iupac!(b'G', 0b0100);
    set_iupac!(b'T', 0b1000);
    
    // Ambiguity Codes (R, Y, S, W, K, M, B, D, H, V)
    set_iupac!(b'R', 0b0101); // A or G
    set_iupac!(b'Y', 0b1010); // C or T
    set_iupac!(b'S', 0b0110); // G or C
    set_iupac!(b'W', 0b1001); // A or T
    set_iupac!(b'K', 0b1100); // G or T
    set_iupac!(b'M', 0b0011); // A or C
    
    set_iupac!(b'B', 0b1110); // C, G, or T (Not A)
    set_iupac!(b'D', 0b1101); // A, G, or T (Not C)
    set_iupac!(b'H', 0b1011); // A, C, or T (Not G)
    set_iupac!(b'V', 0b0111); // A, C, or G (Not T)
    
    // Any Base
    set_iupac!(b'N', 0b1111); // A, C, G, or T

    table
};

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
pub fn from_ascii(nt: u8) -> Result<Self, String> {
        // Direct table lookup. This is the fastest possible conversion.
        let code = IUPAC_LOOKUP_TABLE[nt as usize];
        
        // Check if the lookup resulted in the sentinel value (0b0000)
        // which we defined as an unknown character.
        if code == 0b0000 {
            return Err(format!(
                "Unknown nucleotide: '{}' (ASCII: {}). Valid IUPAC codes are: A, C, G, T, N, and ambiguity codes.",
                nt as char, nt
            ));
        }

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
        0b1111 => 'N',  // never match it

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