use crate::iupac::iupac_to_char;

use pyo3::prelude::*;
use pyo3::types::PyBytes;

/// Represents a potential target sequence (e.g., a protospacer) identified during a genome scan.
/// 
/// This struct is exposed to Python via `pyo3` and holds all necessary information
/// about a single matching site, including its location, orientation, and IUPAC bitmask sequence.
#[pyclass]
#[derive(Clone)]
pub struct Target {
    /// the identifier for the sequence (e.g., chromosome or contig name). Accessible via Python getter
    #[pyo3(get)]
    pub contig: String,

    /// the 0-based starting position of the target sequence within the contig. Accessible via Python getter
    #[pyo3(get)]
    pub position: usize,

    /// the orientation of the target relative to the contig. Accessible via Python getter.
    /// `true` indicates the forward strand; `false` indicates the reverse complement strand
    #[pyo3(get)]
    pub orientation: bool,  // true = foward, false = reverse
    
    // the core sequence data stored as IUPAC bitmasks. Accessed via the custom `target()` getter method
    pub target: Vec<u8>,  // IUPAC bit masks
}

#[pymethods]
impl Target {
    /// Creates a new `Target` instance.
    /// 
    /// This constructor is exposed as the standard way to create the object in Python.
    /// 
    /// # Arguments
    /// * `contig` - The sequence identifier.
    /// * `pos` - The start position.
    /// * `ori` - The orientation (`true` for forward, `false` for reverse).
    /// * `target_seq` - The sequence of IUPAC bitmasks (`Vec<u8>`).
    #[new]
    pub fn new(contig: &str, pos: usize, ori: bool, target_seq: Vec<u8>) -> Self {
        Self {
            contig: contig.to_string(),
            position: pos,
            orientation: ori,
            target: target_seq,
        }
    }

    /// Custom getter method to expose the internal IUPAC bitmask data to Python.
    /// 
    /// This converts the Rust `Vec<u8>` into a Python `bytes` object, which is the 
    /// idiomatic way to handle raw byte sequences in Python.
    /// 
    /// **Access in Python:** `my_target.target()` (or `my_target.target` if `#[getter]` is used).
    /// 
    /// # Arguments
    /// * `py` - The Python interpreter token.
    /// 
    /// # Returns
    /// A Python `bytes` object containing the raw IUPAC bitmasks.
    #[getter]
    pub fn target<'py>(&self, py: Python<'py>) -> Py<PyBytes> {
        // convert the internal Vec<u8> slice into a PyBytes object
        PyBytes::new(py, &self.target).into()
    }

    fn __repr__(&self) -> String {
        // Decode full target (Vec<u8> bitmasks) into IUPAC string
        let mut decoded = String::with_capacity(self.target.len());
        for &b in &self.target {
            decoded.push(iupac_to_char(b));
        }

        format!(
            "<Target object; contig={}, position={}, strand={}, length={}, target={}>",
            self.contig,
            self.position,
            if self.orientation { "+" } else { "-" },
            self.target.len(),
            decoded
        )
    }
}