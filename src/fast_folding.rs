//! TODO: Brief description of RAFFT and data structures
//! TODO: documentation
//! TODO: Sanity checks on associated methods in impl RafftConfig
//! Note that energy parameters and temperature are set globally (available via CLI, crate root and python bindings)

use crate::autocorrelation::*;
use crate::encoding::{BasePairWeights, EncodedSequence, PairTable};
use crate::folding_graph::*;
use crate::vienna::VCompound;

pub struct RafftConfig {
    basepair_weights: BasePairWeights,
    min_unpaired: usize,
    min_loop_energy: f64,
    number_of_lags: usize,
    number_of_branches: usize,
    saved_trajectories: usize,
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
        }
    }
}

impl RafftConfig {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn basepair_weights(mut self, au: f64, gc: f64, gu: f64) -> Self {
        self.basepair_weights = BasePairWeights {
            AU: au,
            GC: gc,
            GU: gu,
        };
        self
    }

    pub fn minimum_unpaired_in_hairpins(mut self, min_unpaired: usize) -> Self {
        self.min_unpaired = min_unpaired;
        self
    }

    pub fn minimum_loop_energy(mut self, min_loop_energy: f64) -> Self {
        self.min_loop_energy = min_loop_energy;
        self
    }

    pub fn positional_lags(mut self, number_of_lags: usize) -> Self {
        self.number_of_lags = number_of_lags;
        self
    }

    pub fn maximum_branches(mut self, number_of_branches: usize) -> Self {
        self.number_of_branches = number_of_branches;
        self
    }

    pub fn maximum_trajectories(mut self, saved_trajectories: usize) -> Self {
        self.saved_trajectories = saved_trajectories;
        self
    }

    /// Return an empty [`RafftGraph`] that can be used to construct fast folding trajectories.
    pub fn folding_graph(&self, sequence: &str) -> RafftGraph {
        let fc = VCompound::new(sequence);

        let encoded = EncodedSequence::with_basepair_weights(sequence, &self.basepair_weights)
            .expect("Not a valid RNA Sequence!");

        RafftGraph::new(
            encoded,
            fc,
            self.min_unpaired,
            self.min_loop_energy,
            self.number_of_lags,
            self.number_of_branches,
            self.saved_trajectories,
        )
    }
}

mod tests {
    use super::*;

    #[test]
    fn test_folding() {
        let sequence =
            "GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU";

        let config = RafftConfig::new().maximum_trajectories(5);

        let mut ffgraph = config.folding_graph(sequence);

        ffgraph.construct_trajectories();

        ffgraph.iter().for_each(|node| {
            println!(
                "[{}] {} {:.2}",
                node.depth,
                node.structure.to_string(),
                node.energy as f64 * 0.01
            );
        });
    }
}
