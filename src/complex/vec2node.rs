use petgraph::{graph::NodeIndex, Graph};
use rand::Rng;

pub trait GraphEmbed {
    fn node2vec(&self);
    fn deepwalk(&self, walk_length: usize, num_walks: usize);
}

impl GraphEmbed for Graph<usize, f32> {
    fn deepwalk(&self, walk_length: usize, num_walks: usize) {
        let mut walks = Vec::new();

        for node in self.node_indices() {
            for _ in 0..num_walks {
                let walk = random_walk(&self, node, 20);
                walks.push(walk);
            }
        }

        println!("{:?}", walks);
    }

    fn node2vec(&self) {
        todo!()
    }
}

fn random_walk(graph: &Graph<usize, f32>, start: NodeIndex, walk_length: usize) -> Vec<NodeIndex> {
    let mut rng = rand::rng();
    let mut path = vec![start];

    for _ in 0..walk_length {
        let current = *path.last().unwrap();
        let neighbors: Vec<_> = graph.neighbors(current).collect();

        if neighbors.is_empty() {
            break; // 如果没有邻居，终止游走
        }

        // 随机选择一个邻居
        let next = neighbors[rng.random_range(0..neighbors.len())];
        path.push(next);
    }

    path
}

#[test]
fn main() {
    let mut graph: Graph<usize, f32> = Graph::new(); // The petgraph type must be defined like this, with the definition of edge direction.
    let node1 = graph.add_node(1);
    let node2 = graph.add_node(2);
    let node3 = graph.add_node(3);
    let node4 = graph.add_node(4);
    let node5 = graph.add_node(5);
    let node6 = graph.add_node(6);
    let node7 = graph.add_node(7);
    let node8 = graph.add_node(8);
    let node9 = graph.add_node(9);
    let node10 = graph.add_node(10);

    graph.add_edge(node1, node2, 1.0);
    graph.add_edge(node1, node3, 1.0);
    graph.add_edge(node1, node4, 1.0);
    graph.add_edge(node2, node5, 1.0);
    graph.add_edge(node2, node6, 1.0);
    graph.add_edge(node3, node7, 1.0);
    graph.add_edge(node3, node8, 1.0);
    graph.add_edge(node4, node9, 1.0);
    graph.add_edge(node4, node10, 1.0);
    graph.add_edge(node5, node6, 1.0);
    graph.add_edge(node7, node8, 1.0);

    graph.deepwalk(20, 5);
}
