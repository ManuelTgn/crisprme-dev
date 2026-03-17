use std::ffi::CStr;

use columnar::{Columnar, Schema};
use columnar::python::PyBufferFormat;

use crate::crispr::guide::Guide;
use crate::model::cigarx::Cigarx64;
use crate::model::input::{SeqBatch, SeqFrame, SeqId, SeqOccFrame};
use crate::model::occurence::Occurence;
use crate::sequence::iupac::Iupac;

/// Max length of the resolved guides and sequence
pub const ALIGN_RESOLVED_MAX_LEN: usize = 32;

/// Definition of a mined alignment
#[derive(Debug, Columnar)]
pub struct SeqMined {
    /// Unique identifier of the source sequence
    pub seq_id: SeqId,
    /// Cigarx that represents the alignment
    pub cigarx: Cigarx64,
    /// Offset from the start of the sequence
    pub offset: u8,
}

/// Definition of a resolved alignment
#[derive(Debug, Columnar)]
pub struct SeqResolved {
    /// Unique identifier of the source sequence
    pub seq_id: SeqId,
    /// Resolved guide
    pub rguide: [u8; ALIGN_RESOLVED_MAX_LEN],
    /// Resolved sequence
    pub rseq: [u8; ALIGN_RESOLVED_MAX_LEN],
    /// Offset from the start of the sequence
    pub offset: u8,
}

/// Definition of a complete alignment
#[derive(Debug, Columnar)]
pub struct Alignment {

    /// Unique identifier of the source sequence
    pub seq_id: SeqId,

    /// Resolved guide
    pub rguide: [u8; ALIGN_RESOLVED_MAX_LEN],
    /// Resolved sequence
    pub rseq: [u8; ALIGN_RESOLVED_MAX_LEN],

    /// Where this alignment occures
    pub occurence: Occurence,
    /// Offset from the start of the sequence
    pub offset: u8,

    /// Features
    #[columnar(group)]
    pub features: [u32; 10],

    /// Scores
    #[columnar(group)]
    pub scores: [f32; 4],
}

/// Inform columnar that iupac characters are just u8 in python
impl PyBufferFormat for Iupac {
    const FORMAT: &'static CStr = unsafe { 
        CStr::from_bytes_with_nul_unchecked("B\0".as_bytes()) 
    };
}

pub struct SeqMinedBatch {

    pub guide: Guide,

    pub sequences: SeqFrame,
    pub occurences: SeqOccFrame,
    pub mined: SeqMinedFrame,
}

pub struct SeqResolvedBatch {

    pub occurences: SeqOccFrame,
    pub resolved: SeqResolvedFrame,
}

pub struct AlignmentBatch {
    content: AlignmentFrame
}
/*
/// Metadata for a batch of mined alignments
#[derive(Debug)]
pub struct MinedBatchMetadata {
    /// Source sequence batch
    pub sequences: BatchRef<sequences::SeqSchema, SeqBatchMetadata>,
}

/// Metadata for a batch of resolved alignments
#[derive(Debug)]
pub struct ResolvedBatchMetadata {
    /// All occurences of the resolved sequences
    pub occurences: Vec<BatchRef<SeqOccSchema, ()>>,
}
 */
