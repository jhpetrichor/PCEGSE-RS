use ndarray::prelude::*;
use rand::prelude::*;
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::collections::HashMap;

// Function to generate a random walk on a graph
fn generate_random_walk(start_node: usize, graph: &HashMap<usize, Vec<usize>>, walk_length: usize, p: f64, q: f64) -> Vec<usize> {
    let mut walk = vec![start_node];
    let mut prev_node = start_node;
    let mut rng = thread_rng();

    for _ in 1..walk_length {
        let neighbors = graph.get(&prev_node).unwrap();
        let mut weights = neighbors.iter().map(|&node| {
            let weight = if node == walk[walk.len() - 2] {
                1.0 / p
            } else if graph.get(&node).unwrap().contains(&walk[walk.len() - 2]) {
                1.0
            } else {
                1.0 / q
            };
            weight
        }).collect::<Vec<f64>>();

        let sum_weights = weights.iter().sum::<f64>();
        weights = weights.iter().map(|&weight| weight / sum_weights).collect();

        let next_node = *neighbors.choose_weighted(&mut rng, &weights).unwrap();
        walk.push(next_node);
        prev_node = next_node;
    }

    walk
}

// Function to generate multiple random walks on a graph
fn generate_random_walks(graph: &HashMap<usize, Vec<usize>>, num_walks_per_node: usize, walk_length: usize, p: f64, q: f64) -> Array2<usize> {
    let num_nodes = graph.len();
    let mut walks = Array2::zeros((num_nodes, num_walks_per_node * walk_length));

    for node in 0..num_nodes {
        for i in 0..num_walks_per_node {
            let walk = generate_random_walk(node, graph, walk_length, p, q);
            for j in 0..walk_length {
                walks[[node, i * walk_length + j]] = walk[j];
            }
        }
    }

    walks
}

fn main() {
    // Define a simple graph with 4 nodes and 5 edges
    let mut graph = HashMap::new();
    graph.insert(0, vec![1, 2]);
    graph.insert(1, vec![0, 2, 3]);
    graph.insert(2, vec![0, 1, 3]);
    graph.insert(3, vec![1, 2]);

    // Generate 10 random walks of length 80 for each node using node2vec
    let num_walks_per_node = 10;
    let walk_length = 80;
    let p = 1.0;
    let q = 1.0;
    let walks = generate_random_walks(&graph, num_walks_per_node, walk_length, p, q);

    // Train a skip-gram model on the walks to obtain node embeddings
    // ...
}