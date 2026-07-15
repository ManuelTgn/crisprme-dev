use columnar::{Columnar, Schema};

use crate::alignment::thresholds::Thresholds;
use crate::model::occurence::Occurence;
use crate::{crispr::guide::Guide, sequence::iupac::Iupac};

/// Maximum length of a sequence.
/// Must be a power of 2 so that sizeof([Iupac; SEQ_MAX_LEN]) divides CHUNK_SIZE (65536) evenly.
pub const SEQ_MAX_LEN: usize = 32;

/// Type that defines the unique sequence row inside a SeqFrame
pub type SeqRowIdx = u32;

/// Definition of a single unique sequence
#[derive(Debug, Columnar)]
pub struct Seq {
    /// IUPAC elements that compose the sequence
    pub content: [Iupac; SEQ_MAX_LEN],
}

/// Definition of a single sequence occurrence
#[derive(Debug, Columnar)]
pub struct SeqOcc {
    /// Identifier for the owning sequence row
    pub seq_row_idx: SeqRowIdx,
    /// Where this sequence occures, packed (contig_id, position, strand)
    pub occurence: Occurence,
}

/// Metadata for a batch of sequences
pub struct SeqBatch {
    /// Length of the sequences
    pub seq_len: usize,
    /// PAM length; protospacer ends at seq_len - pam_len
    pub pam_len: usize,
    /// Guide used for the alignment process
    pub guide: Guide,

    /// Thresholds to use to mine this batch
    pub thresholds: Thresholds,

    /// All sequences of this batch
    pub sequences: SeqFrame,
    /// All occurences relative to this sequence batch
    pub occurences: SeqOccFrame,
}
