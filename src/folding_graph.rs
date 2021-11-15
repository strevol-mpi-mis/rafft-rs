//! This module provides the core algorithm of RAFFT to construct fast folding graphs.
use crate::encoding::{EncodedSequence, PairTable};
use crate::vienna::VCompound;
use itertools::Itertools;
use petgraph::graph::DiGraph;
use std::collections::{HashMap, HashSet};

pub use petgraph::graph::NodeIndex;

/// Information stored per Node in a `RafftGraph`
//#[derive(Clone, Eq, Hash, PartialEq)]
pub struct RafftNodeInfo {
    /// Encoded subsequences for this structure,
    pub sub_nodes: Vec<EncodedSequence>,
    /// structure of this node, corresponds to its parent's structure + stack gained through the corresponding edge
    pub structure: PairTable,
    /// cached free energy of the structure
    pub energy: i32,
    /// depth of the node, i.e. number of edges starting from the root node
    pub depth: usize,
}

/// Fast-folding graph containing the folding trajectories and associated information.
pub struct RafftGraph {
    pub(crate) inner: DiGraph<RafftNodeInfo, ()>,
    node_table: HashMap<String, NodeIndex>,
    root: NodeIndex,
    fc: VCompound,
    min_unpaired: usize,
    min_loop_energy: f64,
    number_of_lags: usize,
    number_of_branches: usize,
    saved_trajectories: usize,
}

impl RafftGraph {
    /// Construct new graph containing only the root node
    pub fn new(
        root: EncodedSequence,
        fold_compound: VCompound,
        min_unpaired: usize,
        min_loop_energy: f64,
        number_of_lags: usize,
        number_of_branches: usize,
        saved_trajectories: usize,
    ) -> Self {
        let mut inner = DiGraph::new();
        let mut node_table = HashMap::new();

        let root_structure = PairTable::new(root.len());
        let root_string = root_structure.to_string();

        let root_info = RafftNodeInfo {
            sub_nodes: vec![root],
            structure: root_structure,
            energy: 0,
            depth: 0,
        };

        let _root = inner.add_node(root_info);
        node_table.insert(root_string, _root);

        Self {
            inner,
            node_table,
            root: _root,
            fc: fold_compound,
            min_unpaired,
            min_loop_energy,
            number_of_lags,
            number_of_branches,
            saved_trajectories,
        }
    }

    /// Return the `NodeIndex` of the root node.
    pub fn root(&self) -> NodeIndex {
        self.root
    }

    /// Insert a new structure as child of `parent`.
    /// If the structure is already present, the `NodeIndex` of the existing node is returned.
    /// A new edge is added anyway if there was not already an edge starting from `parent`.
    /// Therefore, a `RafftGraph` is usually not a tree.
    pub fn insert(
        &mut self,
        parent: NodeIndex,
        sub_nodes: Vec<EncodedSequence>,
        structure: PairTable,
        energy: i32,
    ) -> NodeIndex {
        let depth = self.inner[parent].depth + 1;

        let structure_string = structure.to_string();

        let info = RafftNodeInfo {
            sub_nodes,
            structure,
            energy,
            depth,
        };

        let node_index = if let Some(index) = self.node_table.get(&structure_string) {
            *index
        } else {
            let index = self.inner.add_node(info);
            self.node_table.insert(structure_string, index);

            index
        };

        self.inner.update_edge(parent, node_index, ());
        node_index
    }

    /// Return whether the fast folding graph already contains a structure with the provided dot-bracket notation.
    pub fn contains(&self, structure: &str) -> bool {
        self.node_table.get(structure).is_some()
    }
}

impl RafftGraph {
    /// Construct folding trajectories recursively in a breadth-first fashion, starting from the root.
    pub fn construct_trajectories(&mut self) {
        let current_nodes = vec![self.root()];
        self.breadth_first_search(&current_nodes);
    }

    /// Return an iterator over all structures represented as [`&RafftNodeInfo`] and
    /// in insertion order (i.e. breadth-first and sorted by energy).
    pub fn iter(&self) -> impl Iterator<Item = &RafftNodeInfo> + '_ {
        self.inner.node_weights()
    }

    /// Return the (Metropolis) transition rates `r(i->j) = min(1, exp(-beta * (dGj - dGi)))` between each pair of connected structures,
    /// with `beta = k * T`, where `k` is the Boltzmann constant and `T` the absolute temperature.
    /// The output is in COO format (`(r, i, j)`) and should work nicely with [`scipy.sparse.coo_matrix()`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html#scipy.sparse.coo_matrix).
    /// The indices correspond to the order of [`iter()`] and the output of `rufft` as well as the python bindings.
    ///
    /// Try `beta = 0.61`.
    #[allow(non_snake_case)]
    pub fn transition_rates(&self, beta: f64) -> (Vec<f64>, Vec<usize>, Vec<usize>) {
        // TODO need to store diaganol entries such that rows sum to 0
        let nnz = self.inner.edge_count() + self.inner.node_count();
        let mut data: Vec<f64> = Vec::with_capacity(nnz);
        let mut is: Vec<usize> = Vec::with_capacity(nnz);
        let mut js: Vec<usize> = Vec::with_capacity(nnz);

        // diagonal entries
        for i in 0..self.inner.node_count() {
            data.push(0.0);
            is.push(i);
            js.push(i);
        }

        for edge in self.inner.raw_edges() {
            let dGi = self.inner[edge.source()].energy;
            let dGj = self.inner[edge.target()].energy;

            let ddG = (dGj - dGi) as f64;

            let rate = 1.0f64.min((-beta * ddG).exp());
            let rev_rate = 1.0f64.min((beta * ddG).exp());

            data.push(rate);
            is.push(edge.source().index());
            js.push(edge.target().index());

            data.push(rev_rate);
            is.push(edge.target().index());
            js.push(edge.source().index());

            // diagonal entries
            data[edge.source().index()] -= rate;
            data[edge.target().index()] -= rev_rate;
        }

        (data, is, js)
    }

    /// Return the directed edges `(i, j)` of the fast folding graph, where `i`, `j` are
    /// the indices of the participating structures.
    /// This is the sparse COO format with empty weights.
    /// The indices correspond to the order of [`iter()`] and the output of `rufft` as well as the python bindings.
    pub fn adjacent_indices(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.inner
            .raw_edges()
            .iter()
            .map(|edge| (edge.source().index(), edge.target().index()))
    }

    /// Recursively construct the fast folding graph layer-per-layer.
    ///
    /// _Implementation Detail_: This does not need to be done recursively. In fact, this
    /// implementation is easily translated into an iterative style as the underlying `Vec`-backed
    /// graph structure is not closely tied to the algorithm.
    #[allow(clippy::type_complexity)]
    fn breadth_first_search(&mut self, nodes: &[NodeIndex]) {
        // Using iterators nested in a for-loop because
        // nested iterators and borrowing still is elusive to me.
        // Also, triple-nested Vec is probably not very efficient
        // `all_children` is a Vec containing all children of the current step
        // (i.e. for all `structure_id`s' sub-nodes all the inner and outer fragments
        // that were generated based on global parameters
        let mut all_children: Vec<
            Vec<
                Vec<(
                    Option<EncodedSequence>,
                    Option<EncodedSequence>,
                    PairTable,
                    i32,
                )>,
            >,
        > = Vec::with_capacity(nodes.len());

        for structure_id in nodes {
            let energy = self.inner[*structure_id].energy;
            let pt = self.inner[*structure_id].structure.clone();

            all_children.push(
                self.inner[*structure_id]
                    .sub_nodes
                    .iter()
                    .filter_map(|encoded| {
                        let children = self.create_children(encoded, energy, &pt);
                        match children.len() {
                            0 => None,
                            _ => Some(children),
                        }
                    })
                    .collect::<Vec<_>>(),
            );
        }

        let mut i_branch = 0;

        // unfortunately I seem to need this because I don't want to insert first and then remove unnecessary nodes?
        // in the reference implementation this gets passed down during recursion
        // but I think I can leave it locally for now
        let mut seen: HashSet<String> = HashSet::with_capacity(self.number_of_branches);

        // parent, sub_nodes, structure, energy
        let mut new_children: Vec<(NodeIndex, Vec<EncodedSequence>, PairTable, i32)> =
            Vec::with_capacity(self.number_of_branches);

        for (structure_id, node_children) in nodes.iter().zip(all_children.iter()) {
            for combined_helix in node_children
                .iter()
                .map(|inner| inner.iter())
                .multi_cartesian_product()
            {
                let mut sub_nodes: Vec<EncodedSequence> = vec![];
                let mut pt = PairTable::new(self.fc.len());

                for helix_part in combined_helix {
                    helix_part
                        .2
                        .paired()
                        .for_each(|(i, j)| pt.insert(i as i16, j as i16));

                    if let Some(inner) = &helix_part.0 {
                        sub_nodes.push(inner.clone());
                    }

                    if let Some(outer) = &helix_part.1 {
                        sub_nodes.push(outer.clone());
                    }
                }

                let structure_string = pt.to_string();

                if !self.contains(&structure_string) && seen.insert(structure_string) {
                    i_branch += 1;

                    let energy = self.fc.evaluate_structure(pt.view());
                    new_children.push((*structure_id, sub_nodes, pt, energy));
                }

                if i_branch >= self.number_of_branches {
                    break;
                }
            }
        }

        // sort by energy
        new_children.sort_by_key(|child| child.3);
        new_children = new_children[..self.saved_trajectories.min(new_children.len())].to_vec();

        let new_nodes: Vec<NodeIndex> = new_children
            .into_iter()
            .map(|(parent, sub_nodes, pt, energy)| self.insert(parent, sub_nodes, pt, energy))
            .collect();

        if !new_nodes.is_empty() {
            self.breadth_first_search(&new_nodes);
        }
    }

    fn create_children(
        &self,
        parent_fragment: &EncodedSequence,
        reference_energy: i32,
        parent_structure: &PairTable,
    ) -> Vec<(
        Option<EncodedSequence>,
        Option<EncodedSequence>,
        PairTable,
        i32,
    )> {
        let corr = parent_fragment.autocorrelation(1.0);
        let mut corr = corr.iter().enumerate().collect::<Vec<_>>();
        corr.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap()); // swapping a and b saves me from using `corr.reverse();`

        let mut children: Vec<_> = corr
            .iter()
            .take(self.number_of_lags)
            .filter_map(|(lag, _)| {
                let (bp, mi, mj, _score) =
                    parent_fragment.consecutive_pairs_at_lag(*lag, self.min_unpaired);

                if bp > 0 {
                    let mut pt = parent_structure.clone();

                    (0..bp).for_each(|i| {
                        pt.insert(
                            parent_fragment.parent_indices[mi - i] as i16,
                            parent_fragment.parent_indices[mj + i] as i16,
                        );
                    });

                    let energy = self.fc.evaluate_structure(pt.view());

                    if (energy - reference_energy) as f64 * 0.01 < self.min_loop_energy {
                        let inner = if mj - mi > 1 {
                            Some(parent_fragment.subsequence(mi + 1, mj))
                        } else {
                            None
                        };

                        let outer = if mi + 1 > bp || mj + bp < parent_fragment.len() {
                            Some(parent_fragment.subsequence(mj + bp, mi - bp + 1))
                        } else {
                            None
                        };

                        Some((inner, outer, pt, energy))
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect();
        // Usually, sorting children seems to be unnecessary
        // but there are cases where it makes a difference.
        children.sort_by_key(|child| child.3);

        children
    }
}
