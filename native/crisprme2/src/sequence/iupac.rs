use bytemuck::{Pod, Zeroable};

/// Lookup table mapping ASCII nucleotide characters to 4-bit IUPAC bitmasks.
///
/// Each entry encodes the set of possible bases represented by an IUPAC code
/// using the following bit assignments:
///
/// - A -> 0b0001
/// - C -> 0b0010
/// - G -> 0b0100
/// - T -> 0b1000
///
/// Ambiguous IUPAC codes are represented as the bitwise OR of their possible
/// bases (e.g., `R = A|G = 0b0101`).
///
/// The value `0b1111` corresponds to `N` (any base).
///
/// All non-IUPAC or unrecognized characters map to `0b0000`, which is used as
/// a sentinel value to signal invalid input during parsing.
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
    
    set_iupac!(b'A', 0b0001);
    set_iupac!(b'C', 0b0010);
    set_iupac!(b'G', 0b0100);
    set_iupac!(b'T', 0b1000);
    
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


/// Represents a nucleotide encoded as a 4-bit IUPAC ambiguity mask.
///
/// This compact representation allows constant-time matching via bitwise
/// operations and supports both standard and degenerate nucleotide codes.
///
/// # Bitmask encoding
/// - A: `0b0001`
/// - C: `0b0010`
/// - G: `0b0100`
/// - T: `0b1000`
///
/// Ambiguous codes are encoded as the bitwise OR of their possible bases.
/// For example:
/// - R (A or G): `0b0101`
/// - N (any base): `0b1111`
#[repr(transparent)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Iupac(u8);  // bit mask wrapper

unsafe impl Zeroable for Iupac { }
unsafe impl Pod for Iupac { }
unsafe impl Send for Iupac { }

impl Iupac {
    #[inline(always)]
    pub const fn new(mask: u8) -> Self {
        Self(mask)
    }

    #[inline(always)]
    pub const fn as_u8(self) -> u8 {
        self.0
    }

    #[inline(always)]
    pub fn try_from_ascii(value: u8) -> Option<Self> {
        let code = IUPAC_LOOKUP_TABLE[value as usize];
        if code == 0 { None } else { Some(Self(code)) }
    }

    #[inline(always)]
    pub fn from_ascii_lossy(value: u8) -> Self {
        Self::try_from_ascii(value).unwrap_or_else(|| Self(0b1111))
    }

    #[inline(always)]
    pub fn from_ascii_strict(value: u8) -> Self {
        Self::try_from_ascii(value).expect("invalid IUPAC ASCII character")
    }

    #[inline(always)]
    pub fn complement(self) -> Self {
        let mut code = 0;
        if self.0 & 0b0001 != 0 { code |= 0b1000; }  // A -> T
        if self.0 & 0b0010 != 0 { code |= 0b0100; }  // C -> G 
        if self.0 & 0b0100 != 0 { code |= 0b0010; }  // G -> C
        if self.0 & 0b1000 != 0 { code |= 0b0001; }  // T -> A

        if self.0 & 0b0101 != 0 { code |= 0b1010; }  // R -> Y 
        if self.0 & 0b1010 != 0 { code |= 0b0101; }  // Y -> R 
        if self.0 & 0b0110 != 0 { code |= 0b0110; }  // S -> S 
        if self.0 & 0b1001 != 0 { code |= 0b1001; }  // W -> W
        if self.0 & 0b1100 != 0 { code |= 0b0011; }  // K -> M 
        if self.0 & 0b0011 != 0 { code |= 0b1100; }  // M -> K 

        if self.0 & 0b1110 != 0 { code |= 0b0111; }  // B -> V
        if self.0 & 0b1101 != 0 { code |= 0b1011; }  // D -> H
        if self.0 & 0b1011 != 0 { code |= 0b1101; }  // H -> D
        if self.0 & 0b0111 != 0 { code |= 0b1110; }  // V -> B

        if self.0 & 0b1111 != 0 { code |= 0b1111; }  // N -> N
        
        Self(code)
    }

    #[inline(always)]
    pub fn to_ascii(self) -> u8 {
       match self.0 {

            0b0001 => b'A',
            0b0010 => b'C',
            0b0100 => b'G',
            0b1000 => b'T',

            0b0101 => b'R',
            0b1010 => b'Y',
            0b0110 => b'S',
            0b1001 => b'W',
            0b1100 => b'K',
            0b0011 => b'M',
            
            0b1110 => b'B',
            0b1101 => b'D',
            0b1011 => b'H',
            0b0111 => b'V',
            
            0b1111 => b'N',
            _ => b'?',  // invalid / unknown code
        }
    }

    #[inline(always)]
    pub fn to_ascii_lowercase(self) -> u8 {
       match self.0 {

            0b0001 => b'a',
            0b0010 => b'c',
            0b0100 => b'g',
            0b1000 => b't',

            0b0101 => b'r',
            0b1010 => b'y',
            0b0110 => b's',
            0b1001 => b'w',
            0b1100 => b'k',
            0b0011 => b'm',
            
            0b1110 => b'b',
            0b1101 => b'd',
            0b1011 => b'h',
            0b0111 => b'v',
            
            0b1111 => b'n',
            _ => b'?',  // invalid / unknown code
        }
    }

    #[inline(always)]
    fn from_utf8_lossy(value: char) -> Self {
        Self::from_ascii_lossy(value as u8)
    }

    #[inline(always)]
    pub fn from_utf8(value: char) -> Self {
        Self::from_utf8_lossy(value)
    } 

    #[inline(always)]
    pub fn to_utf8(self) -> char {
        self.to_ascii() as char
    }

    #[inline(always)]
    pub fn matches(self, other: Self) -> bool {
        (self.0 & other.0) != 0
    }

    /// Returns true if exactly one base (A/C/G/T).
    #[inline(always)]
    pub fn is_pure(self) -> bool {
        self.0.count_ones() == 1
    }

    /// Returns true if this is `N` (wildcard).
    #[inline(always)]
    pub fn is_wildcard(self) -> bool {
        self.0 == 0b1111
    }

    /// Number of possible bases represented by this code (1..4).
    #[inline(always)]
    pub fn mutation_score(self) -> u32 {
        self.0.count_ones()
    }
}


/// Checks whether a nucleotide bitmask matches a pattern bitmask.
///
/// In IUPAC semantics, a match occurs if the set of possible bases encoded
/// by the nucleotide overlaps with the set encoded by the pattern.
///
/// This is implemented as a bitwise AND operation and is the core primitive
/// used for both PAM and guide matching.
///
/// # Arguments
/// * `nt` - Bitmask of the sequence nucleotide
/// * `pattern` - Bitmask of the pattern/template nucleotide
///
/// # Returns
/// * `true` if there is at least one common base
/// * `false` otherwise
pub fn matches_iupac(nt: u8, pattern: u8) -> bool {
    (nt & pattern) != 0
}


pub fn sequence_encoder(sequence: &str) -> Vec<u8> {
    sequence
        .as_bytes()
        .iter()
        .map(|&b| Iupac::from_ascii_lossy(b).as_u8())
        .collect()
}

pub fn sequence_encoder_strict(sequence: &str) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(sequence.len());
    for &b in sequence.as_bytes() {
        out.push(Iupac::try_from_ascii(b)?.as_u8());
    }
    Some(out)
}