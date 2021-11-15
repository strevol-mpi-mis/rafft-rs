use crate::fast_folding::RafftConfig;
use crate::folding_graph::RafftGraph;
use pyo3::prelude::*;

pub(crate) fn register(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<FastFoldingGraph>()?;
    Ok(())
}

// TODO: see if VCompound in vienna.rs is safe to send before removing `unsendable`
#[pyclass(module = "rafft", unsendable)]
struct FastFoldingGraph {
    inner: RafftGraph,
}

#[pymethods]
impl FastFoldingGraph {
    #[new]
    #[args(
        number_of_lags = "100",
        number_of_branches = "1000",
        saved_trajectories = "1",
        au = "2.0",
        gc = "3.0",
        gu = "1.0",
        min_unpaired = "3",
        min_loop_energy = "0.0"
    )]
    #[allow(clippy::too_many_arguments)]
    fn new(
        sequence: &str,
        number_of_lags: usize,
        number_of_branches: usize,
        saved_trajectories: usize,
        au: f64,
        gc: f64,
        gu: f64,
        min_unpaired: usize,
        min_loop_energy: f64,
    ) -> Self {
        let config = RafftConfig::new()
            .maximum_trajectories(saved_trajectories)
            .basepair_weights(au, gc, gu)
            .minimum_unpaired_in_hairpins(min_unpaired)
            .minimum_loop_energy(min_loop_energy)
            .maximum_branches(number_of_branches)
            .positional_lags(number_of_lags);

        FastFoldingGraph {
            inner: config.folding_graph(sequence),
        }
    }

    fn trajectories(&mut self) -> PyResult<Vec<(usize, String, f64)>> {
        self.inner.construct_trajectories();

        let trajectories = self
            .inner
            .iter()
            .map(|node| {
                (
                    node.depth,
                    node.structure.to_string(),
                    node.energy as f64 * 0.01,
                )
            })
            .collect();

        Ok(trajectories)
    }

    #[args(beta = "0.61")]
    fn transition_rates(&self, beta: f64) -> PyResult<(Vec<f64>, Vec<usize>, Vec<usize>)> {
        Ok(self.inner.transition_rates(beta))
    }

    fn directed_edges(&self) -> PyResult<(Vec<usize>, Vec<usize>)> {
        let (is, js): (Vec<_>, Vec<_>) = self.inner.adjacent_indices().unzip();
        Ok((is, js))
    }
}
