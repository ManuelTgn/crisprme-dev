use columnar::macros::Columnar;
use columnar::pool::BatchRef;
use crate::model::input::{SeqBatchMetadata, SeqId, SeqOccSchema};
use crate::sequence::iupac::Iupac;

/// Max length of the resolved guides and sequence
pub const ALIGN_RESOLVED_MAX_LEN: usize = 30;
/// Maximum number of features of an alignment
pub const ALIGN_MAX_FEATURES: usize = 10;
/// Maximum number of scores of an alignment
pub const ALIGN_MAX_SCORES: usize = 4;

/// Unique identifier for an alignment
pub type AlignmentId = u32;

pub mod mined {
    use columnar::buffer::Schema;
    use super::*;

    /// Definition of a mined alignment
    #[derive(Debug, Columnar)]
    pub struct Mined {

        /// Unique identifier of the source sequence
        pub seq_id: SeqId,

        /// Offset from the start of the sequence
        pub offset: u8,
    }
}

pub mod resolved {
    use columnar::buffer::Schema;
    use super::*;

    /// Definition of a resolved alignment
    #[derive(Debug, Columnar)]
    pub struct Resolved {

        /// Unique identifier of the source sequence
        pub seq_id: SeqId,

        /// Resolved guide
        pub rguide: [Iupac; ALIGN_RESOLVED_MAX_LEN],
        /// Resolved sequence
        pub rseq: [Iupac; ALIGN_RESOLVED_MAX_LEN],

        /// Length of the resolved sequence and guide
        pub resolved_len: u8,

        /// Offset from the start of the sequence
        pub offset: u8,
    }
}

pub mod aligned {
    use columnar::buffer::Schema;
    use crate::model::input::Occurrence;

    use super::*;

    /// Definition of a complete alignment
    #[derive(Debug, Columnar)]
    pub struct Alignment {

        /// Unique identifier for this particular alignment
        pub id: AlignmentId,

        /// Resolved guide
        pub rguide: [Iupac; ALIGN_RESOLVED_MAX_LEN],
        /// Resolved sequence
        pub rseq: [Iupac; ALIGN_RESOLVED_MAX_LEN],

        /// Length of the resolved sequence and guide
        pub resolved_len: u8,

        /// Offset from the start of the sequence
        pub offset: u8,

        /// Where this alignment occures
        pub occurence: Occurrence,

        /// Features
        #[columnar(group)]
        pub features: [u32; ALIGN_MAX_FEATURES],

        /// Scores
        #[columnar(group)]
        pub scores: [f32; ALIGN_MAX_SCORES],
    }
}

/// Metadata for a batch of mined alignments
#[derive(Debug)]
pub struct MinedBatchMetadata {
    pub sequences: SeqBatchMetadata,
}

/// Metadata for a batch of resolved alignments
#[derive(Debug)]
pub struct ResolvedBatchMetadata {
    /// All occurences of the resolved sequences
    pub occurences: Vec<BatchRef<SeqOccSchema, ()>>,
}


pub use mined::MinedSchema;
pub use resolved::ResolvedSchema;
pub use aligned::AlignmentSchema;