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

    /// Returns the length of the encoded sequence.
    pub fn len(&self) -> usize {
        self.forward.len_of(Axis(1))
    }

    /// Get an copy-on-write slice of a subsequence (0-indexed).
    /// The range defined by `start` and `end` is exclusive.
    /// If `start >= end`, a contiguous EncodedSequence is newly created, with `end` as `5'` and `start-1` as `3'`.
    pub fn subsequence(&'a self, start: usize, end: usize) -> Self {
        if start < end {
            let sub_fwd = self.forward.slice(s![.., start..end]);
            let sub_mrrd = self.mirrored.slice(s![.., start..end]);

            Self {
                forward: CowArray::from(sub_fwd),
                mirrored: CowArray::from(sub_mrrd),
            }
        } else {
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
            }
        }
    }
}

impl<'a> EncodedSequence<'a> {
    /// Search for the longest sequence of consecutive pairs of the encoded sequence and its (reversed) mirror
    /// offset-aligned by `positional_lag` using a sliding-window approach.
    ///
    /// Returns a quadruple containing the number of pairs in the sequence,
    /// the first and last position of that sequence, and a score based on the underlying [`BasePairWeights`]
    pub fn consecutive_pairs_at_lag(&self, positional_lag: usize) -> (usize, usize, usize, usize) {
        // Slicing this way since self.mirrored is stored in the same direction as self.forward
        // Maybe this would be simpler using `%`?
        let (fwd_sliceinfo, mrrd_sliceinfo) = if positional_lag < self.len() {
            (s![.., ..=positional_lag], s![.., ..=positional_lag;-1])
        } else {
            (
                s![.., self.len() - positional_lag + 1..],
                s![.., self.len()-positional_lag+1..;-1],
            )
        };

        let fwd_slice = self.forward.slice(fwd_sliceinfo);
        let mrrd_slice = self.mirrored.slice(mrrd_sliceinfo);

        // window only needs to slide over half of the offset-aligned sequences
        // I don't think I need pos_list to check for conti
        todo!()
    }
}

/// A wrapper type for pair tables in `ViennaRNA`.
/// This struct stores `i16` internally and is `1`-indexed.
///
/// Refer to the [upstream API](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__struct__utils__pair__table.html) for details.
pub struct PairTable(Array1<i16>);

impl PairTable {
    /// Creates a new [`PairTable`].
    pub fn new(length: usize) -> Self {
        let mut inner = Array1::zeros(length + 1);
        inner[0] = length.try_into().unwrap();
        PairTable(inner)
    }

    /// Returns the `length` of the represented structure.
    /// The internal representation has `length + 1` elements for compatibility with ViennaRNA.
    pub fn len(&self) -> usize {
        self.0[0] as usize
    }

    /// Returns `true` if the represented structure is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns an iterator over all unpaired positions (`1`-indexed).
    pub fn unpaired(&self) -> impl Iterator<Item = usize> + '_ {
        self.0
            .indexed_iter()
            .filter(|(_, &u)| u == 0)
            .map(|(i, _)| i as usize)
    }

    /// Returns an iterater over all ordered tuples of paired positions (`1-indexed`).
    pub fn paired(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.0
            .indexed_iter()
            .skip(1)
            .filter(|(i, &u)| i < &(u as usize))
            .map(|(i, &u)| (i, u as usize))
    }

    /// Inserts a new pair into the [`PairTable`].
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

    /// Returns a view of the inner array.
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
}
