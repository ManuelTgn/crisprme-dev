use bytemuck::{Pod, Zeroable};


#[repr(transparent)]
#[derive(Debug, Clone, Copy)]
pub struct Occurence(u64);

impl Occurence {
    pub fn new(contig: u32, position: u32, strand: u8) -> Self {
        Self(((contig as u64) << 33) | ((position as u64) << 1) | ((strand as u64) & 1))
    }
}

unsafe impl Zeroable for Occurence { }
unsafe impl Pod      for Occurence { }