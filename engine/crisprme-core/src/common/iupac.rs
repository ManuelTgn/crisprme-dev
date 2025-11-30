use bytemuck::NoUninit;

/// Represents a single IUPAC character
#[repr(transparent)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Iupac(u8);

/// No problem sending a byte between threads
unsafe impl Send for Iupac { }

impl Iupac {

    /// Create Iupac from already encoded byte
    pub fn new(byte: u8) -> Self {
        Self(byte)
    }

    /// Convert an ascii character to IUPAC
    pub fn from_ascii(value: u8) -> Self {
        let code = match value.to_ascii_uppercase() {

            b'A' => 0b0001, // A
            b'C' => 0b0010, // C
            b'G' => 0b0100, // G
            b'T' => 0b1000, // T

            b'R' => 0b0101, // A or G
            b'Y' => 0b1010, // C or T
            b'S' => 0b0110, // G or C
            b'W' => 0b1001, // A or T
            b'K' => 0b1100, // G or T
            b'M' => 0b0011, // A or C
            
            b'B' => 0b1110, // C or G or T
            b'D' => 0b1101, // A or G or T
            b'H' => 0b1011, // A or C or T
            b'V' => 0b0111, // A or C or G
            
            b'N' => 0b1111, // any base
            v => {
                println!("unknown ASCII: {}", v as char);
                unimplemented!()
            },
        };
        Iupac(code)
    }
    
    /// Complement of the IUPAC code (A <-> T, C <-> G).
    pub fn complement(self) -> Iupac {
        let mut code = 0;
        if self.0 & 0b0001 != 0 { code |= 0b1000; } // A -> T
        if self.0 & 0b0010 != 0 { code |= 0b0100; } // C -> G
        if self.0 & 0b0100 != 0 { code |= 0b0010; } // G -> C
        if self.0 & 0b1000 != 0 { code |= 0b0001; } // T -> A
        Iupac(code)
    }

    /// Convert Iupac to ascii character
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
            _ => {
                println!("invalid conversion from Iupac to ASCII: {:08b}", self.0);
                unimplemented!()
            }, // invalid / unknown code
        }
    }

    /// Convert an utf8 character to IUPAC
    #[inline(always)]
    pub fn from_utf8(value: char) -> Self {
        Iupac::from_ascii(value as u8)
    }

    /// Convert Iupac to ascii character
    #[inline(always)]
    pub fn to_utf8(self) -> char {
        self.to_ascii() as char
    }

    /// Returns true if `self` and `other` share at least one base.
    #[inline(always)]
    pub fn matches(self, other: Iupac) -> bool {
        (self.0 & other.0) != 0
    }

    /// Returns true if this is exactly one base (A, C, G, or T).
    #[inline(always)]
    pub fn is_pure(self) -> bool {
        self.0.count_ones() == 1
    }

    /// Returns true if `self` is a wildcard (N)
    #[inline(always)]
    pub fn is_wildcard(self) -> bool {
        self.0 == 0b1111
    }
    
    /// Get the mutation score
    #[inline(always)]
    pub fn mutation_score(self) -> u32 {
        self.0.count_ones()
    }
}

#[cfg(test)]
mod test {
    use super::Iupac;

    const BITS: [u8; 15] = [
        0b0001,
        0b0010,
        0b0100,
        0b1000,
        0b0101,
        0b1010,
        0b0110,
        0b1001,
        0b1100,
        0b0011,
        0b1110,
        0b1101,
        0b1011,
        0b0111,
        0b1111,
    ];

    const ASCII: [u8; 15] = [
        b'A',
        b'C',
        b'G',
        b'T',
        b'R',
        b'Y',
        b'S',
        b'W',
        b'K',
        b'M',
        b'B',
        b'D',
        b'H',
        b'V',
        b'N',
    ];

    const UTF8: [char; 15] = [
        'A',
        'C',
        'G',
        'T',
        'R',
        'Y',
        'S',
        'W',
        'K',
        'M',
        'B',
        'D',
        'H',
        'V',
        'N',
    ];

    #[test]
    fn encode_ascii() {
        for (i, ascii) in ASCII.iter().enumerate() {
            let iupac = Iupac::from_ascii(*ascii);
            assert_eq!(iupac.0, BITS[i]);
        }
    }

    #[test]
    fn encode_utf8() {
        for (i, c) in UTF8.iter().enumerate() {
            let iupac = Iupac::from_utf8(*c);
            assert_eq!(iupac.0, BITS[i]);
        }
    }

    #[test]
    fn decode_ascii() {
        for (i, bits) in BITS.iter().enumerate() {
            let iupac = Iupac(*bits);
            assert_eq!(iupac.to_ascii(), ASCII[i]);
        }
    }

    #[test]
    fn decode_utf8() {
        for (i, bits) in BITS.iter().enumerate() {
            let iupac = Iupac(*bits);
            assert_eq!(iupac.to_utf8(), UTF8[i]);
        }
    }
}
