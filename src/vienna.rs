//! This module provides some limited functionality of ViennaRNA for use in RAFFT.
use librna_sys::{
    vrna_eval_structure_pt, vrna_fold_compound, vrna_fold_compound_free, vrna_fold_compound_t,
    vrna_md_defaults_temperature, vrna_md_t, vrna_params_load, VRNA_OPTION_EVAL_ONLY,
    VRNA_PARAMETER_FORMAT_DEFAULT, VRNA_VERSION,
};
use ndarray::ArrayView1;
use std::ffi::CString;
use std::path::PathBuf;

/// A wrapper struct around `vrna_fold_compound_t` from ViennaRNA with limited functionality,
pub struct VCompound {
    fc: *mut vrna_fold_compound_t,
    length: usize,
}

impl VCompound {
    /// Create a new `VCompound` wrapper object from a string representing an RNA sequence.
    pub fn new(sequence: &str) -> Self {
        let csequence = CString::new(sequence).expect("CString::new failed");
        let fc = unsafe {
            let md = std::ptr::null::<vrna_md_t>();

            vrna_fold_compound(csequence.as_ptr(), md, VRNA_OPTION_EVAL_ONLY)
        };

        Self {
            fc,
            length: sequence.len(),
        }
    }

    /// Compute the minimum free energy of an RNA secondary structure provided as a pair table.
    /// Pair tables are 1-indexed and contain the structure's length at position 0.
    pub fn evaluate_structure(&self, pairtable: ArrayView1<i16>) -> i32 {
        assert_eq!(pairtable.len(), self.length + 1);
        unsafe { vrna_eval_structure_pt(self.fc, pairtable.as_ptr()) }
    }

    /// Compute the minimum free energy of an RNA secondary structure provided as a pair table.
    /// Pair tables are 1-indexed and contain the structure's length at position 0.
    /// **This returns `[evaluate_structure()] * 0.01` (`kcal/mol`)**
    pub fn evaluate_structure_f64(&self, pairtable: ArrayView1<i16>) -> f64 {
        self.evaluate_structure(pairtable) as f64 * 0.01
    }

    /// Return the length of the underlying sequence of the fold compound.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.length
    }
}

impl Drop for VCompound {
    fn drop(&mut self) {
        unsafe {
            vrna_fold_compound_free(self.fc);
        }
    }
}

// TODO: check safety for this
// https://medium.com/dwelo-r-d/wrapping-unsafe-c-libraries-in-rust-d75aeb283c65
// see bindings.rs #[pyclass(unsendable)]
//unsafe impl Send for VCompound {}
//unsafe impl Sync for VCompound {}

/// Set the temperature of the Nearest-Neighbor model in `ViennaRNA` globally.
/// Refer to the [upstream API](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__model__details.html#gaf9e527e9a2f7e6fd6e42bc6e602f5445) for details.
pub fn set_global_temperature(temperature: f64) {
    unsafe {
        vrna_md_defaults_temperature(temperature);
    }
}

/// Read the parameters of the Nearest-Neighbor model from a file and sets them globally.
/// Refer to the [upstream API](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__energy__parameters__rw.html#gabb0583595c67094986ef90cb4f1c7555) for details.
pub fn set_global_energy_parameters(parameters: PathBuf) {
    let cparams =
        std::ffi::CString::new(parameters.to_str().unwrap()).expect("CString::new failed");
    unsafe {
        vrna_params_load(cparams.as_ptr(), VRNA_PARAMETER_FORMAT_DEFAULT);
    }
}

/// Return the version string of the statically linked `ViennaRNA` library.
// Safety: the version string of ViennaRNA should always be valid Unicode
pub const VIENNA_VERSION: &str = unsafe { std::str::from_utf8_unchecked(VRNA_VERSION) };

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array1;

    #[test]
    fn test_vrna_eval() {
        let sequence =
            "GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU";
        let vc = VCompound::new(sequence);

        // MFE structure:
        //(((((.(((.......))))))))(((.....(((((((((((((.........))).))))))))...))......)))..
        // minimum free energy = -25.80 kcal/mol
        let pt = Array1::from_vec(vec![
            82, 24, 23, 22, 21, 20, 0, 19, 18, 17, 0, 0, 0, 0, 0, 0, 0, 9, 8, 7, 5, 4, 3, 2, 1, 80,
            79, 78, 0, 0, 0, 0, 0, 71, 70, 66, 65, 64, 63, 62, 61, 60, 59, 57, 56, 55, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 45, 44, 43, 0, 42, 41, 40, 39, 38, 37, 36, 35, 0, 0, 0, 34, 33, 0, 0, 0,
            0, 0, 0, 27, 26, 25, 0, 0,
        ]);

        assert_eq!(-25.8f64, vc.evaluate_structure_f64(pt.view()));
    }
}
