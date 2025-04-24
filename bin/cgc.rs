use essential_protein::complex::{
    cgc::{get_complexes, split_graph},
    graph::Graph,
};
use std::collections::VecDeque;

fn main() {
    let collins = Graph::new("collins", false);
    let mut queue = VecDeque::from([collins]);
    let mut splitted_graph = Vec::new();
    while !queue.is_empty() {
        split_graph(&mut queue, &mut splitted_graph);
    }

    let mut result = Vec::new();

    splitted_graph.into_iter().for_each(|g| {
        get_complexes(&g, &mut result, 0.8);
    });

    println!("{}", result.len());
    result.into_iter().for_each(|res| {
        if res.proteins.len() >= 3 {
            res.proteins.into_iter().for_each(|p| {
                print!("{}\t", p);
            });
            println!();
        }
    });

    // splitted_graph.iter().for_each(|g| {
    //     if g.id_protein.len() >= 3 {
    //     for i in g.id_protein.iter() {
    //         print!("{}\t", i);
    //     }
    //     println!();
    // }
    // });
    // println!("{}", splitted_graph.len())
}
