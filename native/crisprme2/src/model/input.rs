use columnar::{Columnar, Schema};

use crate::{crispr::guide::Guide, sequence::iupac::Iupac};
use crate::model::occurence::Occurence;

/// Maximum length of a sequence
pub const SEQ_MAX_LEN: usize = 26;

/// Type that defines the unique sequence Id
pub type SeqId = u32;

/// Definition of a single unique sequence
#[derive(Debug, Columnar)]
pub struct Seq {

    /// Unique identifier for this particular sequence
    pub id: SeqId,
    /// IUPAC elements that compose the sequence
    pub content: [Iupac; SEQ_MAX_LEN],
}

/// Definition of a single sequence occurrence
#[derive(Debug, Columnar)]
pub struct SeqOcc {

    /// Identifier for the owning sequence
    pub seq_id: SeqId,
    /// Where this sequence occures, packed (contig_id, position, strand)
    pub occurence: Occurence,
}

/// Metadata for a batch of sequences
pub struct SeqBatch {

    /// Length of the sequences
    pub seq_len: usize,
    /// Guide used for the alignment process
    pub guide: Guide,

    /// All sequences of this batch
    pub sequences: SeqFrame,
    /// All occurences relative to this sequence batch
    pub occurences: SeqOccFrame,
}
