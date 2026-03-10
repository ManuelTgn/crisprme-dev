use std::sync::Arc;

use columnar::{arena::Arena, pipeline::{Emit, Stage, StageError}, pool::Pool};

use crate::model::alignment::ResolvedSchema;

/// Resolve mined alignments using the present cigarx
pub struct AlignmentResolve {
    
    /// Pool for buffers of resolved alignment schema
    pool: Arc<Pool<ResolvedSchema>>,
    /// Temporary buffer
    arena: Arena,
}

impl Stage for AlignmentResolve {

    type Input  = ();
    type Output = ();

    fn process<E>(&mut self, input: Self::Input, emitter: &mut E) -> Result<(), StageError>
    where
        E: Emit<Self::Output> 
    {
        todo!()
    }
}