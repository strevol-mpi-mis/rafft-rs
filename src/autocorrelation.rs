//! TODO: this extends the module encoding, so not sure if I need module-level documentation here

use crate::encoding::EncodedSequence;
use ndarray::Array1;
use realfft::RealFftPlanner;

fn convolution(a: &[f64], b: &[f64]) -> Array1<f64> {
    let length = a.len() + b.len() - 1;

    let mut planner = RealFftPlanner::<f64>::new();

    let mut in_a = a
        .iter()
        .cloned()
        .chain(std::iter::repeat(0f64))
        .take(length)
        .collect::<Vec<_>>();
    let mut in_b = b
        .iter()
        .cloned()
        .chain(std::iter::repeat(0f64))
        .take(length)
        .collect::<Vec<_>>();

    let fft = planner.plan_fft_forward(length);

    let mut out_a = fft.make_output_vec();
    let mut out_b = fft.make_output_vec();

    fft.process(&mut in_a, &mut out_a).unwrap();
    fft.process(&mut in_b, &mut out_b).unwrap();

    let mut in_ab = out_a
        .iter()
        .zip(out_b.iter())
        .map(|(ai, bi)| ai * bi / length as f64)
        .collect::<Vec<_>>();
    let ifft = planner.plan_fft_inverse(length);
    let mut out_ab = ifft.make_output_vec();

    ifft.process(&mut in_ab, &mut out_ab).unwrap();

    Array1::from_vec(out_ab)
}

impl<'a> EncodedSequence<'a> {
    /// TODO find a simpler way?
    /// explain `padding`
    pub fn autocorrelation(&self, padding: f64) -> Array1<f64> {
        let correlates = self
            .forward
            .rows()
            .into_iter()
            .zip(self.mirrored.rows().into_iter())
            .map(|(f, m)| convolution(f.as_slice().unwrap(), m.as_slice().unwrap()))
            .collect::<Vec<_>>();

        let shape = correlates[0].dim();

        let norm = Array1::from_iter(
            (0..(shape + 1) / 2)
                .chain((0..(shape + 1) / 2 - 1).rev())
                .map(|norm| norm as f64 + padding),
        );

        let mut correlates = correlates
            .iter()
            .fold(Array1::zeros(shape), |acc, c| acc + c);

        correlates.zip_mut_with(&norm, |c, n| *c /= *n);
        correlates
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::*;

    #[test]
    fn test_convolution() {
        let a = vec![
            0., 0., 0., 1., 1., 1., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
            0., 0., 0., 1., 0., 1., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 1., 0., 0., 1., 0., 1.,
            0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 0., 0., 0., 1.,
        ];

        let bb = vec![
            1., 1., 1., 0., 0., 0., 1., 0., 1., 1., 0., 1., 0., 2., 2., 1., 0., 1., 0., 2., 1., 0.,
            0., 0., 1., 0., 0., 0., 0., 2., 0., 2., 0., 0., 1., 0., 1., 0., 1., 1., 0., 2., 0., 2.,
            1., 1., 0., 2., 0., 0., 2., 1., 0., 2., 0., 0., 1., 2., 0., 1., 0., 0., 1., 0., 2., 0.,
            2., 0., 2., 1., 1., 1., 0., 0., 0., 0., 0., 1., 2., 0., 2., 0.,
        ];

        /*let b = vec![
            0., 2., 0., 2., 1., 0., 0., 0., 0., 0., 1., 1., 1., 2., 0., 2., 0., 2., 0., 1., 0., 0.,
            1., 0., 2., 1., 0., 0., 2., 0., 1., 2., 0., 0., 2., 0., 1., 1., 2., 0., 2., 0., 1., 1.,
            0., 1., 0., 1., 0., 0., 2., 0., 2., 0., 0., 0., 0., 1., 0., 0., 0., 1., 2., 0., 1., 0.,
            1., 2., 2., 0., 1., 0., 1., 1., 0., 1., 0., 0., 0., 1., 1., 1.,
        ];*/

        let mut ab_conv = convolution(&a, &bb);
        ab_conv.mapv_inplace(|f| f.round());

        let ab = Array1::from_vec(vec![
            0., 0., 0., 1., 2., 3., 2., 1., 0., 1., 2., 3., 4., 3., 3., 1., 5., 5., 8., 4., 3., 3.,
            4., 6., 6., 6., 3., 6., 3., 8., 4., 4., 3., 5., 7., 6., 8., 3., 4., 7., 5., 11., 6.,
            6., 8., 6., 9., 8., 9., 7., 7., 7., 6., 7., 10., 10., 9., 11., 7., 10., 10., 10., 9.,
            13., 9., 8., 17., 10., 12., 20., 6., 16., 14., 14., 17., 13., 15., 9., 19., 9., 15.,
            13., 18., 13., 16., 11., 9., 14., 11., 16., 13., 11., 15., 13., 19., 10., 17., 7., 7.,
            14., 5., 9., 11., 11., 10., 19., 12., 8., 12., 6., 9., 8., 10., 10., 8., 12., 10., 12.,
            9., 10., 11., 7., 11., 7., 11., 6., 12., 9., 7., 12., 6., 10., 9., 8., 5., 6., 4., 5.,
            9., 4., 8., 6., 6., 9., 5., 7., 2., 3., 0., 2., 2., 4., 4., 5., 4., 2., 2., 0., 1., 2.,
            0., 2., 0.,
        ]);

        assert_eq!(ab, ab_conv);
    }

    #[test]
    fn test_autocorrelation() {
        let sequence =
            "GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU";
        let bpw = BasePairWeights {
            AU: 2.0,
            GC: 3.0,
            GU: 1.0,
        };
        let encoded = EncodedSequence::with_basepair_weights(sequence, &bpw).unwrap();

        let ac = encoded.autocorrelation(1.0);

        //TODO finish test
        assert!(false);
    }
}
