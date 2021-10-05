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
/// Implementation of the RAFFT fast-folding algorithm.
#[allow(dead_code)]
pub mod fast_folding;
//#[cfg(feature = "rna")]
/// Purpose-specific bindings for ViennaRNA
#[allow(dead_code)]
mod vienna;

pub use vienna::{set_global_energy_parameters, set_global_temperature};

#[cfg(feature = "bindings")]
use pyo3::prelude::*;

/// TODO
#[cfg(feature = "bindings")]
#[pymodule]
fn rafft(_py: Python, m: &PyModule) -> PyResult<()> {
    /// TODO documentation and API
    #[pyfn(m)]
    fn set_temperature(_py: Python, temperature: f64) -> PyResult<()> {
        set_global_temperature(temperature);
        Ok(())
    }

    #[pyfn(m)]
    fn set_energy_parameters(_py: Python, _parameters: &str) -> PyResult<()> {
        set_global_energy_parameters(std::path::PathBuf::from(_parameters));
        Ok(())
    }

    #[pyfn(m)]
    fn fold(_py: Python, _sequence: &str) -> PyResult<(f64, String)> {
        unimplemented!()
        //Ok(crate_fold(s1, s2).unwrap())
    }

    Ok(())
}
