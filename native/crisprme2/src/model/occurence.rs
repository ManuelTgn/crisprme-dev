use bytemuck::{Pod, Zeroable};


#[repr(transparent)]
#[derive(Debug, Clone, Copy)]
pub struct Occurence(pub u64);

impl Occurence {
    pub fn new(contig: u32, position: u32, strand: u8) -> Self {
        Self(((contig as u64) << 33) | ((position as u64) << 1) | ((strand as u64) & 1))
    }

    pub fn contig(&self) -> u32 {
        (self.0 >> 33) as u32
    }

    pub fn position(&self) -> u32 {
        ((self.0 >> 1) & 0x7FFF_FFFF) as u32
    }

    pub fn strand(&self) -> u8 {
        (self.0 & 1) as u8
    }
}

unsafe impl Zeroable for Occurence { }
unsafe impl Pod      for Occurence { }