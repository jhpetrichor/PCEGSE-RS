use core::fmt;
use std::{collections::HashMap, fs::read_to_string, rc::Rc};

#[derive(Debug)]
pub struct Graph {
    pub protein_id: HashMap<String, usize>,
    pub id_protein: Vec<String>,
    pub edges: Vec<(usize, usize, f32)>,
    pub ebd: Rc<HashMap<String, Vec<f32>>>,
}

impl fmt::Display for Graph {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "node_count: {}, edges_count: {}!",
            self.id_protein.len(),
            self.edges.len()
        )
    }
}

impl Graph {
    pub fn new(ppi: &str, weighted: bool) -> Self {
        let ppi_file = format!("data/{}/{}.txt", ppi, ppi);
        let ebd_file = format!("data/{}/{}_ebd.txt", ppi, ppi);

        let (protein_id, id_protein, edges) = read_edges(&ppi_file, weighted);
        let ebd = read_embeding(&ebd_file);

        Self {
            protein_id,
            id_protein,
            edges,
            ebd: Rc::new(ebd),
        }
    }

    // do not stable
    pub fn subgraph<I, S>(&self, nodes: I) -> Self
    where
        I: IntoIterator<Item = S> + Clone,
        S: ToString,
    {
        let id_protein = nodes.into_iter().map(|s| s.to_string()).collect::<Vec<_>>();
        let mut protein_id = HashMap::<String, usize>::new();
        let mut edges = Vec::<(usize, usize, f32)>::new();

        for (k, v) in id_protein.iter().enumerate() {
            protein_id.insert(v.to_string(), k);
        }

        self.edges.iter().for_each(|(a, b, w)| {
            match (
                protein_id.get(&self.id_protein[*a as usize]),
                protein_id.get(&self.id_protein[*b as usize]),
            ) {
                (Some(u), Some(v)) => edges.push((*u, *v, *w)),
                _ => {}
            }
        });

        Self {
            protein_id,
            id_protein,
            edges,
            ebd: self.ebd.clone(),
        }
    }
}

fn read_edges(
    ppi: &str,
    weighted: bool,
) -> (
    HashMap<String, usize>,
    Vec<String>,
    Vec<(usize, usize, f32)>,
) {
    let lines = read_to_string(ppi).expect("Failed to open ppi File!");

    let mut protein_id = HashMap::<String, usize>::new();
    let mut id_protein = Vec::<String>::new();
    let mut edges = Vec::<(usize, usize, f32)>::new();

    lines.lines().for_each(|line| {
        let edge = line.split_whitespace().collect::<Vec<_>>();

        let a = match protein_id.get(edge[0]) {
            Some(id) => *id,
            None => {
                protein_id.insert(edge[0].to_string(), protein_id.len());
                id_protein.push(edge[0].to_string());
                id_protein.len() - 1
            }
        };

        let b = match protein_id.get(edge[1]) {
            Some(id) => *id,
            None => {
                protein_id.insert(edge[1].to_string(), protein_id.len());
                id_protein.push(edge[1].to_string());
                id_protein.len() - 1
            }
        };

        let w = match weighted {
            false => 1.,
            true => edge[2].parse().expect("Faile to parse weight!"),
        };

        edges.push((a, b, w));
    });
    (protein_id, id_protein, edges)
}

fn read_embeding(ebd_file: &str) -> HashMap<String, Vec<f32>> {
    let lines = read_to_string(ebd_file).expect("Failed to open embeding File!");

    let mut ebd = HashMap::<String, Vec<f32>>::new();

    lines.lines().for_each(|l| {
        let line = l.split_whitespace().collect::<Vec<_>>();
        let ebd_values = line[1..]
            .into_iter()
            .map(|c| c.parse::<f32>().expect("Failed to parse embeding values!"))
            .collect::<Vec<f32>>();
        ebd.insert(line[0].to_string(), ebd_values);
    });

    ebd
}

#[cfg(test)]
mod tests {
    use super::Graph;

    #[test]
    fn graph_test() {
        let collins = Graph::new("collins", true);
        println!("{}", collins);

        let nodes = [
            "YAL001C", "YBR123C", "YDR362C", "YGR047C", "YOR110W", "YPL007C",
        ];

        let subgraph = collins.subgraph(nodes);
        println!("{}", subgraph);
    }
}
