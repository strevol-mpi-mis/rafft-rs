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
        let fc = VCompound::new(sequence);

        let encoded = EncodedSequence::with_basepair_weights(sequence, &self.basepair_weights)
            .expect("Not a valid RNA Sequence!");

        todo!()
    }

    pub fn trajectories(&mut self) -> Vec<String> {
        if self.stored_trajectories.is_empty() {
            eprintln!("WARNING: No trajectories stored yet. Consider calling fold() first.");
        }
        todo!()
    }
}

impl<'a> EncodedSequence<'a> {
    /// Find the `n` best consecutive stacks for this `EncodedSequence` that can be formed
    /// by computing first the autocorrelation and then computing the energy reduction
    /// for the `n` best positional lags.
    /// TODO: I need a reference structure that I can augment here
    /// TODO: Ideally I'd already store this in a tree-like structure?
    pub fn best_consecutive_stacks(&self, n: usize) -> bool {

        todo!()
    }
}

mod tests {
    use super::*;
    use crate::encoding::*;

    //#[test]
    fn test_folding() {
        let sequence =
            "GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU";
        let bpw = BasePairWeights {
            AU: 2.0,
            GC: 3.0,
            GU: 1.0,
        };
        let encoded = EncodedSequence::with_basepair_weights(sequence, &bpw).unwrap();
        let fc = VCompound::new(sequence);

        let mut structure = PairTable::new(sequence.len());

        let corr = encoded.autocorrelation(1.0);
        let mut corr = corr.iter().enumerate().collect::<Vec<_>>();
        corr.sort_by(|a, b| a.1.partial_cmp(b.1).unwrap());
        corr.reverse();

        let (bp, mi, mj, score) = encoded.consecutive_pairs_at_lag(corr[1].0, 3);

        (0..bp).for_each(|i| {
            structure.insert((mi - i + 1) as i16, (mj + i + 1) as i16);
        });

        println!("{}", structure.to_string());

        let outer = encoded.subsequence(mj + bp, mi - bp);
        let inner = encoded.subsequence(mi + 1, mj + 1);

        let corr = outer.autocorrelation(1.0);
        let mut corr = corr.iter().enumerate().collect::<Vec<_>>();
        corr.sort_by(|a, b| a.1.partial_cmp(b.1).unwrap());
        corr.reverse();

        println!("{:?}", corr[..10].to_vec());

        let (obp, omi, omj, oscore) = outer.consecutive_pairs_at_lag(corr[3].0, 3);

        (0..obp).for_each(|i| {
            structure.insert((omj + i + 1 - outer.concatenation_site().unwrap()) as i16, (mj + bp + omi + i) as i16);
        });

        println!("{}", structure.to_string());

        let corr = inner.autocorrelation(1.0);
        let mut corr = corr.iter().enumerate().collect::<Vec<_>>();
        corr.sort_by(|a, b| a.1.partial_cmp(b.1).unwrap());
        corr.reverse();

        println!("{:?}", corr[..10].to_vec());

        let (ibp, imi, imj, iscore) = inner.consecutive_pairs_at_lag(corr[1].0, 3);

        (0..ibp).for_each(|i| {
            structure.insert((mi + 1 + imi - i + 1) as i16, (mi + 1 + imj + i + 1) as i16);
        });

        println!("{}", structure.to_string());
    }
}
