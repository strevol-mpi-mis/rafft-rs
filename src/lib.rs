#![warn(missing_docs)]
#![warn(rustdoc::missing_crate_level_docs)]

//! Rust implementation of [`RAFFT`](https://www.biorxiv.org/content/10.1101/2021.07.02.450908v1.full)
//! TODO: crate level documentation

/// Autocorrelation of an encoded RNA sequence using FFT
#[allow(dead_code)]
pub mod autocorrelation;
/// Encoding of RNA sequences using nucleotide representations suitable for FFT
#[allow(dead_code)]
pub mod encoding;
//#[cfg(feature = "rna")]
//#[allow(dead_code)]
//pub mod optionalmodule;

#[cfg(feature = "bindings")]
use pyo3::prelude::*;

/// TODO
#[cfg(feature = "bindings")]
#[pymodule]
fn rafft(_py: Python, m: &PyModule) -> PyResult<()> {
    /// TODO documentation and API
    #[pyfn(m)]
    fn fold(_py: Python, sequence: &str) -> PyResult<(f64, String)> {
        unimplemented!()
        //Ok(crate_fold(s1, s2).unwrap())
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
