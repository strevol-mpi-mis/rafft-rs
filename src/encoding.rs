//! This module provides the central representation of an RNA sequence enabling computation
//! of an autocorrelation using FFT-based convolution.
//! Nucleotides of a sequence are encoded as tuples:
//!
//! - `A = (1, 0, 0, 0)`
//! - `C = (0, 1, 0, 0)`
//! - `G = (0, 0, 1, 0)`
//! - `U = (0, 0, 0, 1)`
//!
//! Additionally, a _mirrored_ copy of the sequence is encoded in reverse using a complementary (in a sense)
//! alphabet, effectively carrying information about the strength of legal base pairs:
//!
//! - `a = (0, 0, 0, AU)`
//! - `c = (0, 0, GC, 0)`
//! - `g = (0, GC, 0, GU)`
//! - `u = (AU, 0, GU, 0)`
//!
//! where `AU`, `GC`, `GU` are weights of the base pairs.

use ndarray::{arr1, s, Array1, Array2, ArrayView1, Axis, CowArray, Ix2};
use std::convert::TryInto;
use thiserror::Error;

/// Error type representing errors that may arise during sequence parsing or encoding.
#[derive(Error, Debug)]
pub enum Error {
    /// Error variant corresponding to invalid nucleotides in the supplied sequence string.
    #[error("invalid nucleotide (expected one of [A, C, G, U], found {0:?})")]
    InvalidNucleotide(char),
}

// emulating an enum with array variants
#[allow(non_snake_case)]
mod Alphabet {
    pub(crate) const A: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
    pub(crate) const C: [f64; 4] = [0.0, 1.0, 0.0, 0.0];
    pub(crate) const G: [f64; 4] = [0.0, 0.0, 1.0, 0.0];
    pub(crate) const U: [f64; 4] = [0.0, 0.0, 0.0, 1.0];
}

/// See the [module-level description](crate::encoding).
#[allow(missing_docs)]
#[allow(non_snake_case)]
pub struct BasePairWeights {
    pub AU: f64,
    pub GC: f64,
    pub GU: f64,
}

#[allow(non_snake_case)]
struct MirrorAlphabet {
    A: Array1<f64>,
    C: Array1<f64>,
    G: Array1<f64>,
    U: Array1<f64>,
}

impl MirrorAlphabet {
    pub fn new(weights: &BasePairWeights) -> Self {
        Self {
            A: arr1(&[0.0, 0.0, 0.0, weights.AU]),
            C: arr1(&[0.0, 0.0, weights.GC, 0.0]),
            G: arr1(&[0.0, weights.GC, 0.0, weights.GU]),
            U: arr1(&[weights.AU, 0.0, weights.GU, 0.0]),
        }
    }
}

impl Default for MirrorAlphabet {
    fn default() -> Self {
        Self {
            A: arr1(&[0.0, 0.0, 0.0, 1.0]),
            C: arr1(&[0.0, 0.0, 1.0, 0.0]),
            G: arr1(&[0.0, 1.0, 0.0, 1.0]),
            U: arr1(&[1.0, 0.0, 1.0, 0.0]),
        }
    }
}

/// An [`EncodedSequence`] consists of a _forward_ encoding and a _mirrored_ encoding.
/// See the [module-level description](crate::encoding) for details.
#[derive(Debug)]
pub struct EncodedSequence<'a> {
    pub(crate) forward: CowArray<'a, f64, Ix2>,
    pub(crate) mirrored: CowArray<'a, f64, Ix2>,
    concatenation_site: Option<usize>,
}

impl<'a> EncodedSequence<'a> {
    /// Encode an RNA sequence with given [`BasePairWeights`] being stored in the mirrored encoded sequence.
    pub fn with_basepair_weights(sequence: &str, weights: &BasePairWeights) -> Result<Self, Error> {
        let mirrored_alphabet = MirrorAlphabet::new(weights);

        let length = sequence.len();

        let mut forward = Array2::default((4, length));
        let mut mirrored = Array2::default((4, length));

        match sequence.chars().enumerate().try_for_each(|(i, c)| match c {
            'A' => {
                forward
                    .column_mut(i)
                    .zip_mut_with(&arr1(&Alphabet::A), |ci, ni| *ci = *ni);
                mirrored
                    .column_mut(i)
                    .zip_mut_with(&mirrored_alphabet.A.view(), |ci, ni| *ci = *ni);

                Ok(())
            }
            'C' => {
                forward
                    .column_mut(i)
                    .zip_mut_with(&arr1(&Alphabet::C), |ci, ni| *ci = *ni);
                mirrored
                    .column_mut(i)
                    .zip_mut_with(&mirrored_alphabet.C.view(), |ci, ni| *ci = *ni);

                Ok(())
            }
            'G' => {
                forward
                    .column_mut(i)
                    .zip_mut_with(&arr1(&Alphabet::G), |ci, ni| *ci = *ni);
                mirrored
                    .column_mut(i)
                    .zip_mut_with(&mirrored_alphabet.G.view(), |ci, ni| *ci = *ni);

                Ok(())
            }
            'U' => {
                forward
                    .column_mut(i)
                    .zip_mut_with(&arr1(&Alphabet::U), |ci, ni| *ci = *ni);
                mirrored
                    .column_mut(i)
                    .zip_mut_with(&mirrored_alphabet.U.view(), |ci, ni| *ci = *ni);

                Ok(())
            }
            _ => Err(Error::InvalidNucleotide(c)),
        }) {
            Err(e) => Err(e),
            _ => Ok(Self {
                forward: CowArray::from(forward),
                mirrored: CowArray::from(mirrored),
                concatenation_site: None,
            }),
        }
    }

    /// Encode an RNA sequence with equal [`BasePairWeights`].
    pub fn new(sequence: &str) -> Result<Self, Error> {
        Self::with_basepair_weights(
            sequence,
            &BasePairWeights {
                AU: 1.0,
                GC: 1.0,
                GU: 1.0,
            },
        )
    }

    /// Return the length of the encoded sequence.
    pub fn len(&self) -> usize {
        self.forward.len_of(Axis(1))
    }

    /// Return whether the encoded sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.forward.is_empty()
    }

    /// Return the position of concatenation if `&self` was created
    /// using [`subsequence()`] from two non-contiguous fragments and `None` otherwise.
    pub fn concatenation_site(&self) -> Option<usize> {
        self.concatenation_site
    }

    /// Return `true` if `&self` represents a sequence concatenated from two fragments
    /// using [`subsequence()`] and `false` otherwise.
    pub fn is_concatenated(&self) -> bool {
        self.concatenation_site.is_some()
    }

    /// Get an copy-on-write slice of a subsequence (0-indexed).
    /// The range defined by `start` and `end` is exclusive.
    /// TODO: If `start >= end`, a contiguous [`EncodedSequence`] is newly created, with `start` as `5'` and `end-1` as `3'`.
    // TODO: vaitea's RAFFT splits without wrapping around (which does not make a difference for pairing per se)
    // TODO I think my approach is more intuitive but vaitea's might make it more easy later on?
    // TODO Also: maybe I don't need to store the actual concatenation site
    // TODO but the relative 5' site (could be negative) in my approach
    pub fn subsequence(&'a self, start: usize, end: usize) -> Self {
        if start < end {
            let sub_fwd = self.forward.slice(s![.., start..end]);
            let sub_mrrd = self.mirrored.slice(s![.., start..end]);

            Self {
                forward: CowArray::from(sub_fwd),
                mirrored: CowArray::from(sub_mrrd),
                concatenation_site: None,
            }
        } else {
            // let indices: Vec<usize> = (0..end).chain(start..self.len())
            // should work as well since it does not change pairing
            // in which case `end` should be stored as concatenation site
            //let indices: Vec<usize> = (0..end).chain(start..self.len()).collect();
            let indices: Vec<usize> = (start..self.len()).chain(0..end).collect();

            // double-select to force C standard layout
            // this is hacky and not as efficient as possible but should suffice for now
            let sub_fwd = self
                .forward
                .select(Axis(1), &indices)
                .select(Axis(0), &[0, 1, 2, 3]);
            let sub_mrrd = self
                .mirrored
                .select(Axis(1), &indices)
                .select(Axis(0), &[0, 1, 2, 3]);

            Self {
                forward: CowArray::from(sub_fwd),
                mirrored: CowArray::from(sub_mrrd),
                //concatenation_site: Some(end),
                //concatenation_site: Some(start),
                concatenation_site: Some(self.len() - start), //TODO: + 1?
            }
        }
    }
}

impl<'a> EncodedSequence<'a> {
    /// Search for the longest sequence of consecutive pairs of the encoded sequence and its (reversed) mirror
    /// offset-aligned by `positional_lag` using a sliding-window approach.
    ///
    /// Returns a quadruple containing the number of pairs in the sequence,
    /// the first paired positions of both strands, and a score based on the underlying [`BasePairWeights`]
    pub fn consecutive_pairs_at_lag(&self, positional_lag: usize) -> (usize, usize, usize, f64) {
        let min_hairpin = 3;

        // Slicing this way since self.mirrored is stored in the same direction as self.forward
        // Maybe this would be simpler using `%`?
        let (fwd_sliceinfo, mrrd_sliceinfo, concat_site) = if positional_lag < self.len() {
            let concat_site = match self.concatenation_site() {
                Some(site) => {
                    // This is not directly apparent, so here are some notes:
                    // If site < positional_lag the window includes the concatenation site.
                    // In this case, we can't accumulate the scores over this site
                    // Since we're only sliding over the first half of the window and the mirrored (reversed) sequence
                    // is, well, mirrored, we have a symmetrical situation
                    // and choose the minimum of site and positional_lag - site to cut off accumulation
                    // Otherwise our offset-aligned sequences do not include the concatenation site and we're fine.
                    if site < positional_lag {
                        Some(site.min(positional_lag - site) + 1)
                    } else {
                        None
                    }
                }
                _ => None,
            };

            (
                s![.., ..=positional_lag],
                s![.., ..=positional_lag;-1],
                concat_site,
            )
        } else {
            let concat_site = match self.concatenation_site() {
                Some(site) => {
                    if site > positional_lag - self.len() + 1 {
                        Some(site.min(site - positional_lag + self.len() - 1))
                    } else {
                        None
                    }
                }
                _ => None,
            };

            (
                s![.., positional_lag - self.len() + 1..],
                s![.., positional_lag - self.len() + 1..;-1],
                concat_site,
            )
        };

        //println!("{}", self.forward);
        //println!("{}", self.mirrored);

        let fwd_slice = self.forward.slice(fwd_sliceinfo);
        let mrrd_slice = self.mirrored.slice(mrrd_sliceinfo);

        //println!("{}", fwd_slice);
        //println!("{}", mrrd_slice);

        // Slide over half of the offset-aligned sequences since they are complementary
        let halved_length = fwd_slice.len_of(Axis(1)) / 2 + fwd_slice.len_of(Axis(1)) % 2;

        // This is not pretty but we can save us the work if that's the case
        if halved_length < min_hairpin {
            return (0, 0, 0, 0.0);
        }

        // The total pairing score per position is computed as the pairwise product
        // of the offset-aligned sequences (actually, only their first halves)
        // and then summed over all four nucleotides.
        let mut total_pairing_scores = (fwd_slice.slice(s![.., ..halved_length]).to_owned()
            * mrrd_slice.slice(s![.., ..halved_length]))
        .sum_axis(Axis(0));

        //println!("{}", total_pairing_scores);
        //println!("{:?} {:?}", concat_site, self.concatenation_site());

        // not very idiomatic but I'm trying to stay close to the reference implementation
        let mut max_score = total_pairing_scores[0];
        let mut acc_pairs = if max_score > 0.0 { 1 } else { 0 };
        let mut max_pairs = acc_pairs;
        let mut max_i = 0;
        let mut i = 0;
        let (mut max_lower, mut max_upper) = if positional_lag < self.len() {
            (0, positional_lag)
        } else {
            (positional_lag - self.len() + 1, self.len() - 1)
        };

        let accumulate_scores = |&prev: &f64, curr: &mut f64| {
            i += 1;

            if concat_site != Some(i) {
                *curr *= prev + *curr;
            } else {
                // if Some(i) == concat_site, reset accumulated pairs
                acc_pairs = 0;
            }

            if *curr > 0.0 {
                acc_pairs += 1;
            } else {
                acc_pairs = 0;
            }

            let (lower_position, upper_position) = if positional_lag < self.len() {
                (i, positional_lag - i)
            } else {
                (positional_lag - self.len() + 1 + i, self.len() - i - 1)
            };

            let distance = match self.concatenation_site() {
                Some(site) => {
                    if lower_position >= site || upper_position < site {
                        upper_position - lower_position
                    } else {
                        // We don't have access to the actual distance but both positions
                        // flank the concatenation site
                        // so there's enough space
                        upper_position - lower_position + min_hairpin
                    }
                }
                _ => upper_position - lower_position,
            };

            if *curr >= max_score
            // check if there are at least 3 unpaired positions between paired positions of stack
            && distance > min_hairpin
            {
                max_score = *curr;
                max_i = i;
                max_upper = upper_position;
                max_lower = lower_position;
                max_pairs = acc_pairs;
            }
        };

        total_pairing_scores.accumulate_axis_inplace(Axis(0), accumulate_scores);

        //println!("{}", total_pairing_scores);
        //println!("{} {} {} {}", max_pairs, max_lower, max_upper, max_score);

        (max_pairs, max_lower, max_upper, max_score)
    }
}

/// A wrapper type for pair tables in `ViennaRNA`.
/// This struct stores `i16` internally and is `1`-indexed.
///
/// Refer to the [upstream API](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__struct__utils__pair__table.html) for details.
pub struct PairTable(Array1<i16>);

impl PairTable {
    /// Create a new [`PairTable`].
    pub fn new(length: usize) -> Self {
        let mut inner = Array1::zeros(length + 1);
        inner[0] = length.try_into().unwrap();
        PairTable(inner)
    }

    /// Return the `length` of the represented structure.
    /// The internal representation has `length + 1` elements for compatibility with ViennaRNA.
    pub fn len(&self) -> usize {
        self.0[0] as usize
    }

    /// Return whether the represented structure is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Return an iterator over all unpaired positions (`1`-indexed).
    pub fn unpaired(&self) -> impl Iterator<Item = usize> + '_ {
        self.0
            .indexed_iter()
            .filter(|(_, &u)| u == 0)
            .map(|(i, _)| i as usize)
    }

    /// Return an iterater over all ordered tuples of paired positions (`1-indexed`).
    pub fn paired(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.0
            .indexed_iter()
            .skip(1)
            .filter(|(i, &u)| i < &(u as usize))
            .map(|(i, &u)| (i, u as usize))
    }

    /// Insert a new pair into the [`PairTable`].
    /// Does not check for crossing pairs.
    /// Panics if supplied positions are out of range or already paired.
    pub fn insert(&mut self, i: i16, j: i16) {
        assert!(0 < i && i <= self.len().try_into().unwrap());
        assert!(0 < j && j <= self.len().try_into().unwrap());

        assert_ne!(i, j);

        assert_eq!(self.0[i as usize], 0);
        assert_eq!(self.0[j as usize], 0);

        self.0[i as usize] = j;
        self.0[j as usize] = i;
    }

    /// Return a view of the inner array.
    pub fn view(&self) -> ArrayView1<i16> {
        self.0.view()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array;

    #[test]
    fn test_encoding() {
        let sequence =
            "GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU";
        let bpw = BasePairWeights {
            AU: 2.0,
            GC: 3.0,
            GU: 1.0,
        };
        let encoded = EncodedSequence::with_basepair_weights(sequence, &bpw).unwrap();

        let fwd = Array::from_vec(vec![
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0., 1., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 1.,
            0., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0.,
            1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0.,
            0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1., 1., 0., 0., 1., 0.,
            0., 0., 1., 0., 1., 1., 0., 0., 0., 1., 0., 0., 1., 0., 1., 0., 0., 0., 1., 0., 1., 0.,
            0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
            1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1., 1., 0., 0., 0., 1., 0., 1., 1., 0., 1.,
            0., 0., 0., 1., 0., 1., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            1., 0., 1., 0., 1., 1., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
            1., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1., 0., 0., 0., 0., 0., 1.,
            0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 1., 0.,
            0., 0., 0., 0., 0., 0., 0., 1., 0., 1., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 1., 0.,
            0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 0., 0., 0., 1.,
        ])
        .into_shape((4, 82))
        .unwrap();

        /*let mrrd = Array::from_vec(vec![
            2., 0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 2., 0., 2., 0., 0., 2.,
            0., 2., 0., 0., 2., 0., 0., 2., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 2., 0., 0., 0., 0., 0., 0., 2., 2., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 2.,
            0., 0., 0., 2., 0., 2., 0., 0., 0., 0., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 3., 0.,
            0., 0., 0., 0., 3., 3., 3., 0., 0., 0., 0., 0., 0., 3., 0., 0., 3., 0., 0., 3., 0., 0.,
            0., 0., 3., 0., 0., 0., 0., 0., 3., 3., 0., 0., 0., 0., 3., 3., 0., 3., 0., 3., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 3., 0., 0., 0., 3., 0., 0., 3., 0., 3., 0., 0., 0., 3., 0.,
            3., 3., 0., 3., 0., 0., 0., 3., 3., 3., 1., 0., 3., 0., 0., 1., 1., 1., 1., 3., 0., 0.,
            0., 0., 3., 0., 1., 0., 1., 0., 3., 1., 0., 1., 0., 0., 1., 3., 0., 1., 0., 0., 1., 3.,
            0., 3., 0., 0., 0., 3., 0., 3., 0., 0., 3., 0., 1., 0., 3., 3., 0., 3., 0., 1., 1., 3.,
            1., 0., 3., 3., 3., 0., 0., 3., 0., 1., 0., 0., 0., 1., 0., 1., 0., 0., 3., 0., 1., 1.,
            1., 0., 0., 0., 0., 2., 0., 2., 1., 0., 0., 0., 0., 0., 1., 1., 1., 2., 0., 2., 0., 2.,
            0., 1., 0., 0., 1., 0., 2., 1., 0., 0., 2., 0., 1., 2., 0., 0., 2., 0., 1., 1., 2., 0.,
            2., 0., 1., 1., 0., 1., 0., 1., 0., 0., 2., 0., 2., 0., 0., 0., 0., 1., 0., 0., 0., 1.,
            2., 0., 1., 0., 1., 2., 2., 0., 1., 0., 1., 1., 0., 1., 0., 0., 0., 1., 1., 1.,
        ])
        .into_shape((4, 82))
        .unwrap();*/
        let mrrd = Array::from_vec(vec![
            0., 0., 0., 2., 2., 2., 0., 0., 0., 0., 2., 0., 2., 0., 0., 0., 2., 0., 0., 0., 0., 0.,
            0., 0., 0., 2., 0., 2., 2., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 2., 0., 0., 2., 0., 0., 2., 0., 0., 2., 0., 2., 0., 0., 2., 0., 2.,
            0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0., 2., 3., 3., 3., 0., 0., 0.,
            3., 0., 3., 3., 0., 3., 0., 0., 0., 3., 0., 3., 0., 0., 3., 0., 0., 0., 3., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 3., 0., 3., 0., 3., 3., 0., 0., 0., 0., 3., 3., 0., 0., 0., 0.,
            0., 3., 0., 0., 0., 0., 3., 0., 0., 3., 0., 0., 3., 0., 0., 0., 0., 0., 0., 3., 3., 3.,
            0., 0., 0., 0., 0., 3., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 0., 3., 0., 0., 1., 0.,
            1., 0., 0., 0., 1., 0., 3., 0., 0., 3., 3., 3., 0., 1., 3., 1., 1., 0., 3., 0., 3., 3.,
            0., 1., 0., 3., 0., 0., 3., 0., 3., 0., 0., 0., 3., 0., 3., 1., 0., 0., 1., 0., 3., 1.,
            0., 0., 1., 0., 1., 3., 0., 1., 0., 1., 0., 3., 0., 0., 0., 0., 3., 1., 1., 1., 1., 0.,
            0., 3., 0., 1., 1., 1., 1., 0., 0., 0., 1., 0., 1., 1., 0., 1., 0., 2., 2., 1., 0., 1.,
            0., 2., 1., 0., 0., 0., 1., 0., 0., 0., 0., 2., 0., 2., 0., 0., 1., 0., 1., 0., 1., 1.,
            0., 2., 0., 2., 1., 1., 0., 2., 0., 0., 2., 1., 0., 2., 0., 0., 1., 2., 0., 1., 0., 0.,
            1., 0., 2., 0., 2., 0., 2., 1., 1., 1., 0., 0., 0., 0., 0., 1., 2., 0., 2., 0.,
        ])
        .into_shape((4, 82))
        .unwrap();

        assert_eq!(encoded.forward, fwd);
        assert_eq!(encoded.mirrored, mrrd);
    }

    #[test]
    fn test_subsequence() {
        let sequence =
            "GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU";
        let bpw = BasePairWeights {
            AU: 2.0,
            GC: 3.0,
            GU: 1.0,
        };
        let encoded = EncodedSequence::with_basepair_weights(sequence, &bpw).unwrap();

        let sub = encoded.subsequence(0, 5);

        // TODO: These assertions are rather implementation-specific.
        assert!(sub.forward.is_view());
        assert!(sub.mirrored.is_view());
        assert_eq!(
            sub.forward,
            CowArray::from(encoded.forward.slice(s![.., 0..5]))
        );
        assert_eq!(
            sub.mirrored,
            CowArray::from(encoded.mirrored.slice(s![.., 0..5]))
        );

        let oligo = "AUGGG";
        let encoded_oligo = EncodedSequence::with_basepair_weights(oligo, &bpw).unwrap();
        let concat_oligo = encoded.subsequence(80, 3);

        assert_eq!(concat_oligo.forward, encoded_oligo.forward);
        assert_eq!(concat_oligo.mirrored, encoded_oligo.mirrored);
    }

    #[test]
    fn test_consecutivepairs() {
        let sequence = "UGCGGUGUAAGUGC";
        let bpw = BasePairWeights {
            AU: 2.0,
            GC: 3.0,
            GU: 1.0,
        };
        let encoded = EncodedSequence::with_basepair_weights(sequence, &bpw).unwrap();

        assert_eq!(encoded.consecutive_pairs_at_lag(25), (0, 0, 0, 0.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(23), (0, 0, 0, 0.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(21), (0, 8, 13, 0.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(16), (1, 3, 13, 3.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(15), (2, 5, 10, 2.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(12), (3, 2, 10, 15.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(9), (1, 0, 9, 2.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(5), (0, 0, 5, 0.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(4), (1, 0, 4, 1.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(3), (0, 0, 0, 0.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(2), (0, 0, 0, 0.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(1), (0, 0, 0, 0.0));
        assert_eq!(encoded.consecutive_pairs_at_lag(0), (0, 0, 0, 0.0));

        // CGGCA ACGUAG GGGUU
        let tobesplit = "CGGCAACGUAGGGGUU";
        let tobesplitenc = EncodedSequence::with_basepair_weights(tobesplit, &bpw).unwrap();

        let splitenc = tobesplitenc.subsequence(11, 5);
        assert_eq!(splitenc.consecutive_pairs_at_lag(6), (1, 1, 5, 9.0));
        assert_eq!(splitenc.consecutive_pairs_at_lag(11), (1, 4, 7, 1.0));
    }
}
