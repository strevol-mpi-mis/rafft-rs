//! TODO: Brief description of RAFFT and data structures
//! TODO: documentation
//! TODO: Sanity checks on associated methods in impl RafftConfig
//! Note that energy parameters and temperature are set globally (available via CLI, crate root and python bindings)

use crate::autocorrelation::*;
use crate::encoding::{BasePairWeights, EncodedSequence};
use crate::vienna::VCompound;

pub struct RafftConfig {
    basepair_weights: BasePairWeights,
    min_unpaired: usize,
    min_loop_energy: f64,
    number_of_lags: usize,
    number_of_branches: usize,
    saved_trajectories: usize,
    stored_trajectories: Vec<String>, // TODO: actual tree structure
}

impl Default for RafftConfig {
    fn default() -> Self {
        Self {
            basepair_weights: BasePairWeights {
                AU: 2.0,
                GC: 3.0,
                GU: 1.0,
            },
            min_unpaired: 3,
            min_loop_energy: 0.0,
            number_of_lags: 100,
            number_of_branches: 1000,
            saved_trajectories: 1,
            stored_trajectories: vec![],
        }
    }
}

impl RafftConfig {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn basepair_weights(&mut self, au: f64, gc: f64, gu: f64) -> &mut Self {
        self.basepair_weights = BasePairWeights {
            AU: au,
            GC: gc,
            GU: gu,
        };
        self
    }

    pub fn minimum_unpaired_in_hairpins(&mut self, min_unpaired: usize) -> &mut Self {
        self.min_unpaired = min_unpaired;
        self
    }

    pub fn minimum_loop_energy(&mut self, min_loop_energy: f64) -> &mut Self {
        self.min_loop_energy = min_loop_energy;
        self
    }

    pub fn positional_lags(&mut self, number_of_lags: usize) -> &mut Self {
        self.number_of_lags = number_of_lags;
        self
    }

    pub fn maximum_branches(&mut self, number_of_branches: usize) -> &mut Self {
        self.number_of_branches = number_of_branches;
        self
    }

    pub fn maximum_trajectories(&mut self, saved_trajectories: usize) -> &mut Self {
        self.saved_trajectories = saved_trajectories;
        self
    }

    pub fn fold(&mut self, sequence: &str) -> (f64, String) {
        todo!()
    }

    pub fn trajectories(&mut self) -> Vec<String> {
        if self.stored_trajectories.is_empty() {
            eprintln!("WARNING: No trajectories stored yet. Consider calling fold() first.");
        }
        todo!()
    }
}
