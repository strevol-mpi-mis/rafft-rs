use crate::encoding::{EncodedSequence, PairTable};
use crate::vienna::VCompound;
use petgraph::graph::DiGraph;
use std::collections::HashMap;

pub use petgraph::graph::NodeIndex;

/// Information stored per Node in a `RafftGraph`
//#[derive(Clone, Eq, Hash, PartialEq)]
pub struct RafftNodeInfo<'a> {
    /// Encoded subsequences for this structure,
    pub sub_nodes: Vec<EncodedSequence<'a>>,
    /// structure of this node, corresponds to its parent's structure + stack gained through the corresponding edge
    pub structure: PairTable,
    /// cached free energy of the structure
    pub energy: i32,
    /// depth of the node, i.e. number of edges starting from the root node
    pub depth: usize,
}

/// Fast-folding graph containing the folding trajectories and associated information.
// TODO: do I need to handle inner  and outer fragments in the graph?
pub struct RafftGraph<'a> {
    pub(crate) inner: DiGraph<RafftNodeInfo<'a>, ()>,
    node_table: HashMap<String, NodeIndex>,
    root: NodeIndex,
    fc: VCompound,
    min_unpaired: usize,
    min_loop_energy: f64,
    number_of_lags: usize,
    number_of_branches: usize,
    saved_trajectories: usize,
}

impl<'a> RafftGraph<'a> {
    /// Construct new graph containing only the root node
    /// with `n` being the sequence length.
    pub fn new(
        root: EncodedSequence<'a>,
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

    /// Insert a new structure as child of `parent` with associated `NodeEdgeInfo` `e`.
    /// If the structure is already present, the `NodeIndex` of the existing node is returned.
    /// A new edge is added anyway if there was not already an edge starting from `parent`.
    /// Therefore, a `RafftGraph` is usually not a tree.
    /// If the edge is already present, its `RafftEdgeInfo` is updated.
    pub fn insert(
        &mut self,
        parent: NodeIndex,
        sub_nodes: Vec<EncodedSequence<'a>>,
        structure: PairTable,
    ) -> NodeIndex {
        let depth = self.inner[parent].depth + 1;

        let structure_string = structure.to_string();

        let info = RafftNodeInfo {
            sub_nodes,
            structure,
            energy: 0,
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
}

impl<'a> RafftGraph<'a> {
    /// Construct folding trajectories recursively in a breadth-first fashion, starting from the root.
    pub fn construct_trajectories(&'a mut self) {
        let current_nodes = vec![self.root()];
        self.breadth_first_search(&current_nodes);
    }

    /// Recursively construct the fast folding graph layer-per-layer.
    fn breadth_first_search(&'a mut self, nodes: &[NodeIndex]) {
        let mut new_nodes: Vec<NodeIndex> = vec![]; // Vec::with_capacity() would be better as soon as I know a good guess for the capacity

        /*let children: Vec<_> = nodes
        .iter()
        .map(|structure_id| {
            let energy = self.inner[*structure_id].energy;
            let pt = self.inner[*structure_id].structure.clone();

            self.inner[*structure_id]
                .sub_nodes
                .iter()
                .map(|encoded| {
                    self.create_children(
                        encoded,
                        energy,
                        &pt,
                    )
                })
                .collect::<Vec<_>>();

            (structure_id)
            /*(
                structure_id,
                ffgraph.inner[*structure_id]
                    .sub_nodes
                    .iter()
                    .map(|encoded| create_children(self.number_of_lags, self.min_unpaired, self.min_loop_energy, fc, encoded, energy, &pt))
                    .collect::<Vec<_>>(),
            )*/
        })
        .collect();*/
        for structure_id in nodes {
            let energy = self.inner[*structure_id].energy;
            let pt = self.inner[*structure_id].structure.clone();

            for encoded in &self.inner[*structure_id].sub_nodes {
                self.create_children(encoded, energy, &pt);

                self.inner[*structure_id]
                    .sub_nodes
                    .iter()
                    .map(|encoded| self.create_children(encoded, energy, &pt))
                    .collect::<Vec<_>>();
            }
        }

        // self.breadth_first_search(fc, ffgraph);
        todo!()
    }

    fn create_children(
        &self,
        parent_fragment: &'a EncodedSequence<'a>,
        reference_energy: i32,
        parent_structure: &PairTable,
    ) -> Vec<(
        Option<EncodedSequence<'a>>,
        Option<EncodedSequence<'a>>,
        PairTable,
        i32,
    )> {
        let corr = parent_fragment.autocorrelation(1.0);
        let mut corr = corr.iter().enumerate().collect::<Vec<_>>();
        corr.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap()); // swapping a and b saves me from using `corr.reverse();`

        corr.iter()
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

                        let outer = if mj - mi > 1 {
                            None //Some(parent_fragment.subsequence(mj + bp, mi - bp))
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
            .collect()
    }
}
