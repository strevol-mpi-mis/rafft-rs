//! This extends the module [`encoding`].
use crate::encoding::EncodedSequence;
use ndarray::Array1;
use realfft::RealFftPlanner;

fn convolution(a: &[f64], b: &[f64]) -> Array1<f64> {
    assert_ne!(a.len(), 0);
    assert_ne!(b.len(), 0);

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

impl EncodedSequence {
    /// Compute the (auto)correlation of an [`EncodedSequence`] with its complementary strand representation using FFT.
    /// A sane value for `padding` is `1.0`.
    /// The `padding` parameter might be removed in the future.
    pub fn autocorrelation(&self, padding: f64) -> Array1<f64> {
        // TODO: remove padding parameter
        assert!(padding > 0.0);

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
    use approx::assert_relative_eq;

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

        let _ac = Array1::from_vec(vec![
            0., 0., 0., 0.5, 0.8, 1., 0.57142857, 1., 0.66666667, 0.8, 0.36363636, 0.5, 0.61538462,
            0.85714286, 0.4, 0.5, 0.94117647, 0.55555556, 1.47368421, 0.7, 0.57142857, 0.54545455,
            1.13043478, 1.25, 1.44, 0.69230769, 0.66666667, 1.28571429, 0.62068966, 1.13333333,
            0.83870968, 1., 1.09090909, 1., 1.08571429, 0.83333333, 0.75675676, 0.63157895,
            1.12820513, 0.95, 0.68292683, 1.52380952, 0.97674419, 1.09090909, 0.62222222,
            0.7826087, 1.14893617, 0.83333333, 1.10204082, 0.64, 0.8627451, 0.61538462, 0.79245283,
            0.48148148, 1.01818182, 0.78571429, 0.84210526, 1.31034483, 0.44067797, 0.73333333,
            0.81967213, 0.80645161, 0.85714286, 0.96875, 0.55384615, 0.60606061, 1.13432836,
            0.73529412, 0.7826087, 1.17142857, 0.50704225, 0.86111111, 0.79452055, 0.7027027,
            0.93333333, 0.57894737, 0.93506494, 0.69230769, 1.16455696, 0.45, 0.81481481,
            0.82926829, 0.88888889, 0.55, 0.93670886, 0.74358974, 0.46753247, 0.84210526,
            0.77333333, 0.75675676, 0.68493151, 0.80555556, 1.09859155, 0.88571429, 0.8115942,
            0.64705882, 0.95522388, 0.48484848, 0.49230769, 1.1875, 0.44444444, 0.67741935,
            0.75409836, 0.86666667, 0.6440678, 1.06896552, 0.84210526, 0.60714286, 0.76363636,
            0.44444444, 0.90566038, 0.76923077, 0.74509804, 0.76, 0.44897959, 0.75, 0.80851064,
            1.30434783, 0.8, 0.72727273, 0.65116279, 0.33333333, 0.53658537, 1.25, 0.87179487,
            0.63157895, 0.81081081, 0.5, 0.57142857, 0.88235294, 0.72727273, 1.1875, 0.77419355,
            0.53333333, 0.55172414, 0.64285714, 0.51851852, 0.61538462, 1.44, 0.33333333,
            0.69565217, 1.09090909, 0.85714286, 1.2, 0.84210526, 0.77777778, 0.23529412, 0.375,
            0.4, 1.14285714, 0.76923077, 0.66666667, 0.72727273, 1., 0.88888889, 0.5, 1.42857143,
            0., 0.4, 1., 0., 2., 0.,
        ]);

        let ac = encoded.autocorrelation(1.0);

        assert_relative_eq!(ac, _ac, epsilon = std::f32::EPSILON as f64);
    }
}
