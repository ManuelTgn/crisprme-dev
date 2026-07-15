pub mod column;
pub mod frame;
pub mod memory;
pub mod pipeline;
pub mod python;
pub mod shared;
pub mod typed;

// Re-export derive macro
pub use columnar_derive::Columnar;

pub use column::{Column, ColumnGroup};
pub use memory::{ChunkArray, MemoryPool};
pub use shared::Share;
pub use typed::{Schema, TypedFrame};

#[cfg(test)]
extern crate self as columnar;

#[cfg(test)]
mod test {
    use crate::memory::{CHUNK_SIZE, region_min_size};

    use super::*;

    #[derive(Columnar)]
    pub struct Example {
        pub a: u32,       // Scalar
        pub b: [u32; 10], // Array
        #[columnar(group)]
        pub c: [u32; 4], // Group
    }

    #[test]
    fn new_allocates_all_columns() {
        let pool = MemoryPool::new(CHUNK_SIZE * 10, |_, _| {});
        let mut frame = ExampleFrame::alloc(&pool, 100);
        frame.with_cols(|c| {
            assert_eq!(c.a.rows(), 100);
            assert_eq!(c.b.rows(), 100);
            assert_eq!(c.c.rows(), 100);
        });
    }

    #[test]
    fn with_cols_shared_and_split() {
        let pool = MemoryPool::new(CHUNK_SIZE * 10, |_, _| {});

        let mut src = ExampleFrame::alloc(&pool, 100);
        let mut dst = ExampleFrame::empty();

        src.with_cols(|mut s| {
            dst.with_cols(|mut d| {
                d.a.shared(&mut s.a);
                d.b.shared(&mut s.b);
                d.c.alloc(&pool, s.c.rows());
            });
        });

        dst.with_cols(|c| {
            let [mut d0, ..] = c.c.split();
            for chunk in d0.chunks_mut() {
                for elem in chunk {
                    *elem = 90;
                }
            }
        });
    }

    #[derive(Columnar)]
    pub struct Sequences {
        pub id: u32,
        pub sequence: [u8; 2],
    }

    #[derive(Columnar)]
    pub struct Positions {
        pub position: u64,
    }

    #[derive(Columnar)]
    pub struct Merged {
        pub seq_id: u32,
        pub sequence: [u8; 2],
        pub position: u64,
    }

    #[test]
    fn merge() {
        let pool = MemoryPool::new(CHUNK_SIZE * 10, |_, _| {});

        let mut sequences = SequencesFrame::alloc(&pool, 100);
        sequences.with_cols(|mut c| {
            for (i, (id, sequence)) in c.id.iter_mut().zip(c.sequence.iter_mut()).enumerate() {
                *id = i as u32;
                *sequence = [0, 1];
            }
        });

        let mut positions = PositionsFrame::alloc(&pool, 100);
        positions.with_cols(|mut c| {
            for (i, p) in c.position.iter_mut().enumerate() {
                *p = (i % 100) as u64;
            }
        });

        // Merge frames
        let mut aggregated = MergedFrame::empty();
        cols!(mut s = sequences, mut p = positions, mut m = aggregated => {

            m.seq_id.shared(&mut s.id);
            m.sequence.shared(&mut s.sequence);
            m.position.shared(&mut p.position);

        });

        // Sequences and positions don't exist anymore, but memory does
        drop(sequences);
        drop(positions);

        // Data still available
        aggregated.with_cols(|c| {
            for (i, id) in c.seq_id.iter().enumerate() {
                assert_eq!(*id, i as u32);
            }
        });
    }
}
