use librna_sys::{
    vrna_alloc, vrna_eval_structure_pt, vrna_fold_compound,
    vrna_fold_compound_free, vrna_fold_compound_t, vrna_md_defaults_temperature,
    vrna_md_set_default, vrna_md_t, vrna_params_load, VRNA_OPTION_EVAL_ONLY,
    VRNA_PARAMETER_FORMAT_DEFAULT,
};
use ndarray::ArrayView1;
use std::convert::TryInto;
use std::ffi::CString;
use std::path::PathBuf;

pub struct VCompound {
    fc: *mut vrna_fold_compound_t,
    length: usize,
}

impl VCompound {
    /// Create a new `VCompound` wrapper object from a string representing an RNA sequence.
    pub fn new(sequence: &str) -> Self {
        let csequence = CString::new(sequence).expect("CString::new failed");
        let fc = unsafe {
            let md: *mut vrna_md_t = vrna_alloc(
                std::mem::size_of::<vrna_md_t>()
                    .try_into()
                    .expect("vrna_md_t allocation failed."),
            ) as *mut vrna_md_t; // TODO: I think vrna_fold_compound_free() frees this
            vrna_md_set_default(md);

            // https://doc.rust-lang.org/std/primitive.pointer.html#method.as_ref
            let _fc = vrna_fold_compound(
                csequence.as_ptr(),
                md.as_ref().unwrap(),
                VRNA_OPTION_EVAL_ONLY,
            );
            _fc
        };

        Self {
            fc,
            length: sequence.len(),
        }
    }

    /// Compute the minimum free energy of an RNA secondary structure provided as a pair table.
    /// Pair tables are 1-indexed and contain the structure's length at position 0.
    pub fn evaluate_structure(&self, pairtable: ArrayView1<i16>) -> f64 {
        assert_eq!(pairtable.len(), self.length + 1);
        let mfe = unsafe { vrna_eval_structure_pt(self.fc, pairtable.as_ptr()) };
        mfe as f64 * 0.01
    }
}

impl Drop for VCompound {
    fn drop(&mut self) {
        unsafe {
            vrna_fold_compound_free(self.fc);
        }
    }
}

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

        assert_eq!(-25.8f64, vc.evaluate_structure(pt.view()));
    }
}
