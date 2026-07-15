use crate::error::crisprme_errors::AnnotationError;

/// Compact feature annotation bitset
/// Internally stored as Vec<u8>
#[derive(Clone)]
pub struct AnnotationBits {
    bits: Vec<u8>,
}

impl AnnotationBits {
    /// Allocate bitset sized for `num_features`
    pub fn new(num_features: usize) -> Self {
        let num_bytes = (num_features + 7) >> 3;
        Self {
            bits: vec![0u8; num_bytes],
        }
    }

    /// Set feature bit
    #[inline(always)]
    pub fn set(&mut self, feature_id: usize) -> Result<(), AnnotationError> {
        // Retrieve bit for query feature
        let byte_index = feature_id >> 3;
        if byte_index >= self.bits.len() {
            return Err(AnnotationError::InvalidFeatureId(feature_id));
        }

        // Flip to true feature's bit
        let bit_offset = feature_id & 7;
        self.bits[byte_index] |= 1u8 << bit_offset;

        Ok(())
    }

    #[inline(always)]
    pub fn as_slice(&self) -> &[u8] {
        &self.bits
    }
}
