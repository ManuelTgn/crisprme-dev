//! Standalone, dependency-free liftover (haplotype -> reference)
pub mod apply;
pub mod chain;
pub mod mapper;

pub use apply::{lift_report_core, LiftReportStats, LIFTOVER_HEADER};
pub use chain::{Block, Chain, ChainFile, Strand};
pub use mapper::{LiftOver, LiftResult, Lifted, MapStatus};
