use crate::graph::Graph;

pub trait Centrality {
    fn degree_centrality(&self) -> Vec<usize>;
    fn degree_centrality_extended(&self) -> Vec<usize>;
}

impl Centrality for Graph {
    fn degree_centrality(&self) -> Vec<usize> {
        let mut degree = self
            .nei_list
            .iter()
            .enumerate()
            .map(|(a, nei)| (a, nei.len()))
            .collect::<Vec<_>>();
        degree.sort_by_key(|a| a.1);
        degree.reverse();
        degree.into_iter().map(|(a, _)| a).collect()
    }

    fn degree_centrality_extended(&self) -> Vec<usize> {
        let mut res = self
            .nei_list
            .iter()
            .enumerate()
            .map(|(a, nei)| (a, nei.values().sum()))
            .collect::<Vec<(usize, f64)>>();
        res.sort_by(|a, b| a.partial_cmp(b).unwrap());
        res.reverse();
        res.into_iter().map(|(a, _)| a).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::graph::Graph;

    use super::Centrality;

    #[test]
    fn test_degree_centrality() {
        let edges = [(1, 2, 0.1), (0, 1, 0.3), (1, 4, 0.5), (2, 4, 0.6)];
        let graph = Graph::new_from(edges.to_vec());
        let res = graph.degree_centrality();
        assert_eq!(res, vec![1, 4, 2, 0, 3]);

        let res = graph.degree_centrality_extended();
        println!("{:?}", res);
    }
}
