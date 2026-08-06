//! Assembly report finalization: sample presence sets, the canonical
//! merged-row schema, and the within/cross-sample mergers.
pub mod row;
pub mod samples;

pub use row::{ClusterKey, MergedRow, NativeLoc, RefLoc, Scores};
pub use samples::{SamplePresence, SampleSet, SampleSetId, SampleSetRegistry, SampleTable};