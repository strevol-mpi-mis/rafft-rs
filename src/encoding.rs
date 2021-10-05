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

use ndarray::{arr1, s, Array1, Array2, Axis, CowArray, Ix2};
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

/// An [EncodedSequence] consists of a _forward_ encoding and a _mirrored_ encoding.
/// See the [module-level description](crate::encoding) for details.
#[derive(Debug)]
pub struct EncodedSequence<'a> {
    pub(crate) forward: CowArray<'a, f64, Ix2>,
    pub(crate) mirrored: CowArray<'a, f64, Ix2>,
}

impl<'a> EncodedSequence<'a> {
    /// Encode an RNA sequence with given [BasePairWeights] being stored in the mirrored encoded sequence.
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
            _ => {
                // I think I don't need to reverse the encoded mirrored sequence because then I'd need to reverse again for autocorrelation::convolution()
                //mirrored.invert_axis(Axis(1));
                Ok(Self {
                    forward: CowArray::from(forward),
                    mirrored: CowArray::from(mirrored),
                })
            }
        }
    }

    /// Encode an RNA sequence with equal [BasePairWeights].
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
            let indices: Vec<usize> = (start..self.forward.len_of(Axis(1)))
                .chain(0..end)
                .collect();

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

    // TODO: method that returns an [EncodedSequence] built from a sub sequence
    // or perform convolution/autocorrelation directly on slices by providing indices and lengths?
    // I think in the latter case it might be better to store `mirrored` _not_ in reverse
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
