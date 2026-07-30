//! PyO3 wrapper exposing the annotation feature registry to Python.
//!
//! The heavy [`FeatureRegistry`] stays Rust-owned. Python only needs the
//! term -> bit mapping so the annotation transform can turn each pysam-fetched
//! term into a bit and OR it into an alignment's `u32` annotation slot.
//!
//! The mapping is handed over **once** as a `dict[str, int]` (see
//! [`PyRegistry::masks`]) and cached Python-side. There is deliberately no
//! per-term lookup method: the vocabulary is tiny (<= 32 terms) and a cached
//! Python dict avoids millions of boundary crossings across a run.

use crate::annotation::features::FeatureRegistry;

use pyo3::prelude::*;

use std::collections::HashMap;

/// Python handle to one annotation BED file's feature registry.
///
/// Build one per BED, then call [`masks`](Self::masks) once to obtain the full
/// term -> bitmask mapping.
#[pyclass(module = "crisprme2._crisprme2_native")]
pub struct PyRegistry {
    inner: FeatureRegistry,
}

#[pymethods]
impl PyRegistry {
    /// Build a registry from a BED file (plain or bgzip-compressed).
    ///
    /// Args:
    ///     path (str): Path to the annotation BED. Compression is detected by
    ///         magic bytes, so `.bed` and `.bed.gz` both work.
    ///
    /// Raises:
    ///     IOError: the file cannot be opened or read.
    ///     ValueError: the BED is malformed, or defines more than 32 distinct
    ///         terms (one per bit of the u32 annotation slot).
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        Ok(Self {
            inner: FeatureRegistry::from_bed(path)?,
        })
    }

    /// Number of distinct annotation terms (bits in use); in `0..=32`.
    fn num_features(&self) -> usize {
        self.inner.num_features()
    }

    /// The complete term -> bitmask mapping as a `dict[str, int]`.
    ///
    /// Each value is a **pre-shifted mask** (`1 << bit`), so the transform ORs
    /// it directly into an alignment's `u32` annotation slot with no shift::
    ///
    ///     masks = registry.masks()          # call ONCE, then cache
    ///     acc = 0
    ///     for term in fetched_terms:        # terms from pysam .fetch()
    ///         acc |= masks[term]
    ///
    /// A fresh dict is allocated on every call, so cache the result rather than
    /// calling this per batch. (If you would rather have raw bit indices than
    /// masks, swap `mask_of` for `bit_of` below — the transform then ORs
    /// `1 << masks[term]`.)
    fn masks(&self) -> HashMap<String, u32> {
        self.inner
            .terms()
            .map(|term| {
                let mask = self
                    .inner
                    .mask_of(term)
                    .expect("a term yielded by terms() is always registered");
                (term.to_string(), mask)
            })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner.num_features()
    }

    fn __repr__(&self) -> String {
        format!("<PyRegistry: {} features>", self.inner.num_features())
    }
}
