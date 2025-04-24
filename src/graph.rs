use std::{
    collections::{HashMap, HashSet, VecDeque},
    fs::read_to_string,
    process::id,
};

use ndarray::NewAxis;
use petgraph::graph;

/// 默认的节点权重
pub const DEFAULT_NODE_WEIGHT: f64 = 1.;
/// 默认的边权重
pub const DEFAULT_EDGE_WEIGHT: f64 = 1.;

#[derive(Debug)]
pub struct Graph {
    pub(crate) node_count: usize,
    pub(crate) edge_count: usize,
    pub(crate) nei_list: Vec<HashMap<usize, f64>>,
    pub(crate) node_weight: Vec<f64>,
    pub(crate) id_protein: Vec<String>,
}

impl Graph {
    pub fn new(n_count: usize) -> Self {
        Self {
            node_count: n_count,
            edge_count: Default::default(),
            nei_list: vec![Default::default(); n_count],
            node_weight: Default::default(),
            id_protein: Default::default(),
        }
    }

    pub fn new_from(edges: Vec<(usize, usize, f64)>) -> Self {
        let node_count = edges
            .iter()
            .map(|(a, b, _)| a.max(b))
            .max()
            .unwrap()
            .to_owned()
            + 1;

        let mut edge_count = 0;
        let mut nei_list = vec![<HashMap<usize, f64>>::new(); node_count];
        edges.clone().into_iter().for_each(|(a, b, w)| {
            nei_list[a].insert(b, w);
            nei_list[b].insert(a, w);
            edge_count += 1;
        });

        Self {
            node_count,
            edge_count,
            nei_list,
            node_weight: vec![DEFAULT_NODE_WEIGHT; node_count],
            id_protein: Default::default(),
        }
    }

    pub fn add_edge(&mut self, a: usize, b: usize, w: f64) {
        if a >= self.node_count || b >= self.node_count {
            panic!("out of bound index");
        }
        self.nei_list[a].insert(b, w);
        self.nei_list[b].insert(a, w);
        self.edge_count += 1;
    }

    pub fn has_edge(&self, a: usize, b: usize) -> bool {
        if a >= self.node_count || b >= self.node_count {
            panic!("out of bound index");
        }
        self.nei_list[a].contains_key(&b)
    }

    pub fn remove_edge(&mut self, a: usize, b: usize) {
        if self.has_edge(a, b) {
            self.nei_list[a].remove(&b);
            self.nei_list[b].remove(&a);
            self.edge_count -= 1;
        }
    }

    // 从节点解析出连边数据
    pub fn subgraph(&self, nodes: &HashSet<String>) -> Self {
        let id_protein: Vec<String> = nodes.clone().into_iter().collect();
        let protein_id = id_protein
            .iter()
            .enumerate()
            .map(|(id, protein)| (protein.to_string(), id))
            .collect::<HashMap<String, usize>>();
        let mut edges_id = Vec::<(usize, usize, f64)>::new();
        self.nei_list.iter().enumerate().for_each(|(a, nei)| {
            nei.iter().for_each(|(b, w)| {
                let a = protein_id.get(&self.id_protein[a]);
                let b = protein_id.get(&self.id_protein[*b]);
                if let (Some(a), Some(b)) = (a, b) {
                    edges_id.push((*a, *b, *w));
                }
            });
        });

        let mut graph: Self = Self::new_from(edges_id);
        graph.id_protein = id_protein;
        graph.edge_count /= 2;

        graph
    }

    pub fn new_from_file(file: &str, weighted: bool) -> Self {
        let contents = read_to_string(file).expect("Failed to open ppi file!");
        let mut pair_protein = Vec::<(usize, usize, f64)>::new();
        let mut protein_id = HashMap::<String, usize>::new();
        let mut id_protein = Vec::<String>::new();
        contents.lines().for_each(|line| {
            let pair = line.split_whitespace().collect::<Vec<_>>();
            let w: f64 = if weighted {
                match pair[2].parse() {
                    Err(_) => 1.,
                    Ok(w) => w,
                }
            } else {
                1.
            };
            // update id and protein
            let a = if let Some(id) = protein_id.get(pair[0]) {
                *id
            } else {
                let id = protein_id.len();
                protein_id.insert(pair[0].to_owned(), id);
                id_protein.push(pair[0].to_owned());
                id
            };

            let b = if let Some(id) = protein_id.get(pair[1]) {
                *id
            } else {
                let id = protein_id.len();
                protein_id.insert(pair[1].to_owned(), id);
                id_protein.push(pair[1].to_owned());
                id
            };
            pair_protein.push((a, b, w));
        });
        let mut graph = Self::new_from(pair_protein);
        graph.id_protein = id_protein;

        graph
    }

    pub fn jaccard(&self, a: usize, b: usize) -> f64 {
        let a_nei = self.nei_list[a].keys().collect::<HashSet<_>>();
        let b_nei = self.nei_list[b].keys().collect::<HashSet<_>>();
        let comsize = a_nei.intersection(&b_nei).count() as f64;
        let capsize = a_nei.union(&b_nei).count() as f64;
        comsize / capsize
    }

    pub fn jaccard_plus(&self, a: usize, b: usize) -> f64 {
        let mut a_nei = self.nei_list[a].keys().collect::<HashSet<_>>();
        let mut b_nei = self.nei_list[b].keys().collect::<HashSet<_>>();
        a_nei.insert(&a);
        b_nei.insert(&b);

        let comsize = a_nei.intersection(&b_nei).count() as f64;
        let capsize = a_nei.union(&b_nei).count() as f64;
        comsize / capsize
    }

    fn compute_core_number(&self) -> Vec<usize> {
        let mut degree: Vec<usize> = self.nei_list.iter().map(|neis| neis.len()).collect();
        let mut core: Vec<usize> = vec![0; self.node_count];
        let mut visited = vec![false; self.node_count];

        // 每次都遍历所有节点，找到当前最小度的未访问节点
        for _ in 0..self.node_count {
            let mut min_deg = usize::MAX;
            let mut min_node = None;
            for node in 0..self.node_count {
                if !visited[node] && degree[node] < min_deg {
                    min_deg = degree[node];
                    min_node = Some(node);
                }
            }

            if let Some(u) = min_node {
                visited[u] = true;
                core[u] = degree[u];

                for (&v, _) in &self.nei_list[u] {
                    if !visited[v] {
                        degree[v] -= 1;
                    }
                }
            }
        }

        core
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};

    use petgraph::graph;

    use super::{Graph, DEFAULT_EDGE_WEIGHT};

    #[test]
    fn test_graph_new_from() {
        let edges = [(1, 2, 0.1), (0, 1, 0.3), (1, 4, 0.5), (2, 4, 0.6)];
        let graph = Graph::new_from(edges.to_vec());
        println!("{:?}", graph);
    }

    #[test]
    fn test_graph_add_edge() {
        let mut graph = Graph::new(5);
        graph.add_edge(1, 2, DEFAULT_EDGE_WEIGHT);
        assert!(graph.has_edge(1, 2));
        graph.add_edge(0, 2, DEFAULT_EDGE_WEIGHT);
        assert!(graph.has_edge(0, 2));
        graph.add_edge(1, 4, DEFAULT_EDGE_WEIGHT);
        assert!(graph.has_edge(1, 4));
        graph.add_edge(2, 3, DEFAULT_EDGE_WEIGHT);
        assert!(graph.has_edge(2, 3));
        assert!(graph.node_count == 5);
        assert!(graph.edge_count == 4);
        println!("{:?}", graph);

        let res = graph.compute_core_number();
        println!("{:?}", res);
    }

    #[test]
    fn test_compute_core_number() {
        let mut graph = Graph::new_from(vec![(0, 1, 0.), (0, 2, 0.), (1, 2, 0.), (0, 3, 0.)]);
        let res = graph.compute_core_number();
        println!("{:?}", res);

        // jaccard and jaccard_plus
        println!("jaccrd: {}", graph.jaccard(0, 3));
        println!("jaccard_plus: {}", graph.jaccard_plus(0, 3));
    }

    #[test]
    fn test_graph_new_from_file() {
        let graph = Graph::new_from_file("./data/collins/collins.txt", true);
        println!("{:?}", graph.id_protein);
        println!("{}, {}", graph.node_count, graph.edge_count);

        let nodes = HashSet::from([
            "YAL001C".to_string(),
            "YBR123C".to_string(),
            "YDR362C".to_string(),
            "YGR047C".to_string(),
            "YOR110W".to_string(),
            "YPL007C".to_string(),
            "YAL002W".to_string(),
            "YLR148W".to_string(),
            "YLR396C".to_string(),
            "YMR231W".to_string(),
            "YPL045W".to_string(),
            "YAL003W".to_string(),
            "YEL034W".to_string(),
            "YGR285C".to_string(),
            "YHR064C".to_string(),
            "YKL081W".to_string(),
            "YLR249W".to_string(),
            "YPL048W".to_string(),
        ]);
        let subgraph = graph.subgraph(&nodes);
        println!("{:?}", subgraph);
    }
}
