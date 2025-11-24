use crate::iupac::matches_iupac;

pub struct ParsedPAM {
    pub bytes: Vec<u8>,
    pub revcomp: Vec<u8>,
}

impl ParsedPAM{
    pub fn new(pam: &str) -> Result<Self, &'static str> {
        let bytes = pam.as_bytes().to_vec();
        let revcomp = bytes.iter()
            .rev()
            .map(|b| complement(*b))
            .collect();

        Ok(Self { bytes, revcomp })
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        // IUPAC complements:
        b'R' => b'Y', b'Y' => b'R',
        b'S' => b'S', b'W' => b'W',
        b'K' => b'M', b'M' => b'K',
        b'B' => b'V', b'V' => b'B',
        b'D' => b'H', b'H' => b'D',
        b'N' => b'N',
        _ => b'N',
    }
}