pub mod sequence {
    use crate::common::Seq;
    use rand::{Rng, SeedableRng};

    /// Create a random sequence as a vector of bytes
    pub fn generate(length: usize, nucleotides: &[u8], seed: u64) -> Seq {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let take_one = || nucleotides[rng.random_range(0..nucleotides.len())]; // Uniform random sequences
        std::iter::repeat_with(take_one).take(length).collect()
    }

    /// Materialize in linear memory a set of windows over the sequence
    pub fn materialize_windows(sequence: &[u8], lenght: usize, delta: usize) -> (Seq, usize) {
        let mut result: Seq = Vec::new();

        let mut beg = 0;
        let mut n = 0;

        while beg < sequence.len() - lenght - 1 {
            result.extend(&sequence[beg..beg + lenght]);

            beg += delta;
            n += 1;
        }

        (result, n)
    }

    /// Calculate the reverse complement of a strand
    pub fn reverse_complement(sequence: &[u8]) -> Seq {
        sequence
            .iter()
            .rev()
            .map(|e| match e {
                
                b'N' => b'N',
                
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',

                b'R' => b'Y',
                b'Y' => b'R',
                b'S' => b'S',
                b'W' => b'W',
                b'K' => b'M',
                b'M' => b'K',

                b'B' => b'V',
                b'D' => b'H',
                b'H' => b'D',
                b'V' => b'B',
                
                _ => unimplemented!(),
            })
            .collect()
    }
}

use std::{io::{self}, path::Path};
use memmap2::Mmap;

/// Threholds for the filtering and mining processes
#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    /// Max allowed gaps in query
    pub qgap: u32,
    /// Max allowed gaps in target
    pub tgap: u32,
    /// Max allowed mismatches
    pub mism: u32,
}

impl Thresholds {
    /// Calculate the max edit distance based on the thresholds
    pub fn ed(&self) -> u32 {
        self.qgap + self.tgap + self.mism
    }
}

#[derive(Debug, Clone, Copy)]
pub struct IUPAC(pub u8);

/// Collection of sequences and ids
pub struct Genome {
    pub sequences: Vec<IUPAC>,
    pub seq_len: usize,
    pub ids: Vec<u32>,
    pub n: usize
}

/// Convert a character in ASCII to a IUPAC
pub fn ascii_to_iupac(value: u8) -> IUPAC {
    match value.to_ascii_uppercase() {
        
        b'A' => IUPAC(0b0001),
        b'C' => IUPAC(0b0010),
        b'G' => IUPAC(0b0100),
        b'T' => IUPAC(0b1000),

        b'R' => IUPAC(0b0101), // A or G
        b'Y' => IUPAC(0b1010), // C or T
        b'S' => IUPAC(0b0110), // G or C
        b'W' => IUPAC(0b1001), // A or T
        b'K' => IUPAC(0b1100), // G or T
        b'M' => IUPAC(0b0011), // A or C
        
        b'B' => IUPAC(0b1110), // C or G or T
        b'D' => IUPAC(0b1101), // A or G or T
        b'H' => IUPAC(0b1011), // A or C or T
        b'V' => IUPAC(0b0111), // A or C or G
        
        b'N' => IUPAC(0b1111), // any base

        _    => unimplemented!(), // invalid
    }
}

/// Convert a slice of IUPAC characters to ASCII 
pub fn iupac_to_ascii_str(iupac: &[IUPAC]) -> String {
    iupac.iter().map(|bits| {
        match bits.0 {
            
            0b0001 => 'A',
            0b0010 => 'C',
            0b0100 => 'G',
            0b1000 => 'T',
            
            0b0101 => 'R', // A or G
            0b1010 => 'Y', // C or T
            0b0110 => 'S', // G or C
            0b1001 => 'W', // A or T
            0b1100 => 'K', // G or T
            0b0011 => 'M', // A or C
            
            0b1110 => 'B', // C or G or T
            0b1101 => 'D', // A or G or T
            0b1011 => 'H', // A or C or T
            0b0111 => 'V', // A or C or G
            
            0b1111 => 'N', // any base

            _      => unimplemented!(), // fallback for invalid bit patterns
        }
    }).collect()
}

/// Returns true if the IUPAC characters are a match
pub fn iupac_match(guide: IUPAC, target: IUPAC) -> bool {
    (guide.0 & target.0) != 0
}

/// A valid guide contains only A,C,G,T and N
pub fn is_valid_guide(guide: &[IUPAC]) -> bool {
    guide.iter().all(|b| 
        b.0 == 0b0001 || 
        b.0 == 0b0010 || 
        b.0 == 0b0100 || 
        b.0 == 0b1000 || 
        b.0 == 0b1111)
}

/// Load a standard <id><sequence> file into memory
pub fn load_standard_from_file<P: AsRef<Path>>(path: P) -> io::Result<Genome> {
    
    let mut sequences: Vec<u8> = vec![];
    let mut ids: Vec<u32> = vec![];

    let file = std::fs::File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };

    let mut seq_len = 0;
    for line in mmap.split(|&b| b == b'\n') {
        if line.is_empty() { continue; }

        // Find tab separator
        if let Some(tab_pos) = line.iter().position(|&b| b == b'\t') {
            
            let id_bytes = &line[..tab_pos];
            let id_str = unsafe { str::from_utf8_unchecked(id_bytes) };
            let id: u32 = id_str.parse().expect("invalid id");

            // Save sequence length
            if seq_len == 0 {
                seq_len = line.len() - (tab_pos + 1);
            }

            // NOTE: Allocation
            sequences.extend_from_slice(&line[tab_pos + 1 ..]);
            ids.push(id);
        }

    }

    // Convert everything to uppercase and IUPAC
    Ok(Genome { 
        sequences: sequences.iter().map(|b| ascii_to_iupac(*b)).collect(), 
        seq_len, 
        n: ids.len(),
        ids, 
    })
}

pub mod visualize {

    pub fn cigar(query: &[u8], target: &[u8], cigar: &[u8], start_pos: usize) {
        let mut qline = String::new();
        let mut mline = String::new();
        let mut cline = String::new();
        let mut tline = String::new();

        let mut qidx = 0;
        let mut tidx = 0;

        // Step 1: Add target prefix (unaligned)
        for i in 0..start_pos {
            tline.push(target[i] as char);
            qline.push(' ');
            cline.push(' ');
            mline.push(' ');
            tidx += 1;
        }

        // Step 2: Alignment
        for op in cigar {
            cline.push(*op as char);
            match op {
                b'M' | b'=' | b'X' => {
                    let qc = query[qidx] as char;
                    let tc = target[tidx] as char;
                    qline.push(qc);
                    tline.push(tc);
                    mline.push(if qc == tc || qc == 'N' { '|' } else { ' ' });
                    qidx += 1;
                    tidx += 1;
                }
                b'D' => {
                    let qc = query[qidx] as char;
                    qline.push(qc);
                    tline.push('-');
                    mline.push(' ');
                    qidx += 1;
                }
                b'I' => {
                    let tc = target[tidx] as char;
                    qline.push('-');
                    tline.push(tc);
                    mline.push(' ');
                    tidx += 1;
                }
                _ => unimplemented!(),
            }
        }

        // Step 3: Add target suffix (unaligned)
        while tidx < target.len() {
            tline.push(target[tidx] as char);
            qline.push(' ');
            mline.push(' ');
            tidx += 1;
        }

        println!("cigarx: {}", cline);
        println!("target: {}", tline);
        println!("        {}", mline);
        println!("query:  {}", qline);
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn reverse_complement() {
        let seq = b"ATTGAGATAGTGTGGGGAAGNGG";
        let rev = super::sequence::reverse_complement(seq);
        assert_eq!(rev, b"CCNCTTCCCCACACTATCTCAAT");
    }
}
