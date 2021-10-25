//! This module provides `RafftConfig`, a convenient wrapper type to construct [`RafftGraph`]s.
//! Note that energy parameters and temperature are set globally (available via CLI, crate root and python bindings)

use crate::encoding::{BasePairWeights, EncodedSequence};
use crate::folding_graph::*;
use crate::vienna::VCompound;

/// A builder type for [`RafftGraph`] allowing to adjust parameters as necessary and to finally construct
/// the graph type per individual RNA sequence.
/// A single `RafftConfig` can be re-used to construct `RafftGraph`s for different sequences.
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
    /// Create a new `RafftConfig` instance with default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set weights of the different legal base pairs.
    /// This affects the autocorrelation computed using FFT.
    pub fn basepair_weights(mut self, au: f64, gc: f64, gu: f64) -> Self {
        self.basepair_weights = BasePairWeights {
            AU: au,
            GC: gc,
            GU: gu,
        };
        self
    }

    /// Set the minimum amount of unpaired positions enclosed by a hairpin loop.
    /// Usually, the default value of `3` is okay.
    pub fn minimum_unpaired_in_hairpins(mut self, min_unpaired: usize) -> Self {
        self.min_unpaired = min_unpaired;
        self
    }

    /// Set the minimum energy value new loops have to contribute in order to be formed.
    pub fn minimum_loop_energy(mut self, min_loop_energy: f64) -> Self {
        self.min_loop_energy = min_loop_energy;
        self
    }

    /// Set the number of positional lags between forward and mirrored encoded RNA sequence strands
    /// that should be searched for base pair stacks.
    pub fn positional_lags(mut self, number_of_lags: usize) -> Self {
        self.number_of_lags = number_of_lags;
        self
    }

    /// Set the number of branches to be explored during construction of the fast folding graph.
    pub fn maximum_branches(mut self, number_of_branches: usize) -> Self {
        self.number_of_branches = number_of_branches;
        self
    }

    /// Set the number of trajectories, or structures per steps, to be saved.
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
    #[test]
    fn test_folding() {
        use super::RafftConfig;
        let sequence =
            "GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU";

        let config = RafftConfig::new().maximum_trajectories(1);

        let mut ffgraph = config.folding_graph(sequence);

        ffgraph.construct_trajectories();

        let trajectory: Vec<_> = ffgraph
            .iter()
            .map(|node| node.structure.to_string())
            .collect();

        assert_eq!(
            trajectory[trajectory.len() - 1],
            "..((((((((((((((.((.....))))))))))))).))).(((.........)))((((((.............))))))"
        );
    }
}
