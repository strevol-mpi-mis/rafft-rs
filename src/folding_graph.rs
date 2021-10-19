use crate::encoding::{EncodedSequence, PairTable};
use petgraph::graph::DiGraph;
use std::collections::HashMap;

pub use petgraph::graph::NodeIndex;

/// Information stored per Edge in a `RafftGraph`
#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub struct RafftEdgeInfo {
    /// number of base pairs being added by transitioning on this edge.
    pub basepairs: usize,
    /// lower index of the innermost base pair of the stack
    pub inner_i: usize,
    /// upper index of the innermost base pair of the stack
    pub inner_j: usize,
    /// base pairing score computed during window sliding
    pub score: usize,
    /// energy change gained by this transition
    pub energychange: i32,
}

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
    pub(crate) inner: DiGraph<RafftNodeInfo<'a>, RafftEdgeInfo>,
    node_table: HashMap<String, NodeIndex>,
    root: NodeIndex,
}

impl<'a> RafftGraph<'a> {
    /// Construct new graph containing only the root node
    /// with `n` being the sequence length.
    pub fn new(root: EncodedSequence<'a>) -> Self {
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
        e: RafftEdgeInfo,
        info: RafftNodeInfo<'a>,
    ) -> NodeIndex {
        let depth = self.inner[parent].depth + 1;

        let structure_string = info.structure.to_string();

        let node_index = if let Some(index) = self.node_table.get(&structure_string) {
            *index
        } else {
            let index = self.inner.add_node(info);
            self.node_table.insert(structure_string, index);

            index
        };

        self.inner.update_edge(parent, node_index, e);
        node_index
    }
}
