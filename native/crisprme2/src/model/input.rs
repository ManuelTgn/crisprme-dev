use columnar::{macros::Columnar, pool::BatchRef};
use crate::{crispr::guide::Guide, sequence::iupac::Iupac};

/// Maximum length of a sequence
pub const SEQ_MAX_LEN: usize = 30;

/// Type that encodes an occurrence (contig_id, position, strand)
pub type Occurrence = u64;

/// Packs an occurrent from parts
pub fn occurence(contig_id: u32, position: u32, strand: u8) -> Occurrence {
    ((contig_id as u64) << 33) | ((position as u64) << 1) | ((strand as u64) & 1)
}

/// Type that defines the unique sequence Id
pub type SeqId = u32;

pub mod sequences {
    use columnar::buffer::Schema;
    use super::*;

    /// Definition of a single unique sequence
    #[derive(Debug, Columnar)]
    pub struct Seq {

        /// Unique identifier for this particular sequence
        pub id: SeqId,

        /// IUPAC elements that compose the sequence
        /// TODO: for now row-major, move to column-major with `columnar(group)`
        pub content: [Iupac; SEQ_MAX_LEN]
    }
}

pub mod occurences {
    use columnar::buffer::Schema;
    use super::*;

    /// Definition of a single sequence occurrence
    #[derive(Debug, Columnar)]
    pub struct SeqOcc {

        /// Identifier for the owning sequence
        pub seq_id: SeqId,

        /// Where this sequence occures, packed (contig_id, position, strand)
        pub occurence: Occurrence,
    }
}

/// Metadata for a batch of sequences
#[derive(Debug)]
pub struct SeqBatchMetadata {

    /// Length of the sequences
    pub seq_len: u32,

    /// Guide used for the alignment process
    pub guide: Guide,

    /// All batches of occurences relative to this sequence batch
    pub occurences: Vec<BatchRef<occurences::SeqOccSchema, ()>>,
}

pub use occurences::SeqOccSchema;
pub use sequences::SeqSchema;