/// 主要写基因表达相关逻辑
use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
};

use crate::graph::Graph;

pub struct GeneExpress {
    // 蛋白质的基因表达谱
    pub(crate) express: HashMap<String, Vec<f64>>,
    // 表达谱对应的均值和方差
    pub(crate) mean_variance: HashMap<String, (f64, f64)>,
}

impl GeneExpress {
    // 读取蛋白质的基因表达谱而非所有
    pub fn new(file: &str, proteins: &HashSet<String>) -> Self {
        let express = read_gene_expression(file, proteins);
        let mean_variance = express
            .iter()
            .map(|(a, b)| {
                let res = get_mean_variance(b);
                (a.to_string(), res)
            })
            .collect::<HashMap<_, _>>();

        Self {
            express,
            mean_variance,
        }
    }

    // EDPIN，为关键蛋白质和非关键蛋白质设置不同的基因活性阈值
    pub fn calculate_active_threshold(&self, eps: &HashSet<String>) -> HashMap<String, f64> {
        let mut threshold = HashMap::<String, f64>::new();

        self.mean_variance.iter().for_each(|(a, (mean, variance))| {
            let f = variance / (1.0 + variance);
            let k = match eps.contains(a) {
                false => 1.,
                true => -1.,
            };
            let active = mean + k * variance.sqrt() * f;
            threshold.insert(a.clone(), active);
        });

        threshold
    }

    fn count(&self) -> usize {
        let mut count = 0;
        for i in self.express.iter() {
            count = i.1.len();
            break;
        }
        count
    }
}

pub fn get_dpins(g: &Graph) -> Vec<Graph> {
    let proteins = g.id_protein.iter().map(|c| c.to_string()).collect();
    let gep = GeneExpress::new("./data/gene-expression.txt", &proteins);
    let eps = read_essential_protein();
    // 计算蛋白质的活跃值
    let active = gep.calculate_active_threshold(&eps);

    // 更新节点
    let mut nodes = vec![HashSet::<String>::new(); gep.count()];
    for p in proteins.iter() {
        // 不含有该蛋白质的基因表达
        if let Some(exp) = gep.express.get(p) {
            for (i, value) in exp.iter().enumerate() {
                if value.ge(active.get(p).unwrap()) {
                    nodes[i].insert(p.clone());
                }
            }
        } else {
            // 在所有子网络中保留
            nodes.iter_mut().for_each(|c| {
                c.insert(p.clone());
            });
        }
    }

    nodes.into_iter().map(|nodes| g.subgraph(&nodes)).collect()
}

pub fn read_essential_protein() -> HashSet<String> {
    let lines = read_to_string("./data/essential proteins.ref")
        .expect("Failed to read essential protein file!");
    lines.lines().map(|c| c.to_string()).collect()
}

fn read_gene_expression(file: &str, eps: &HashSet<String>) -> HashMap<String, Vec<f64>> {
    let lines = read_to_string(file)
        .expect("Failed to read gene expression file!")
        .lines()
        .map(|s| s.to_owned())
        .collect::<Vec<_>>();

    let mut protein_express = HashMap::<String, Vec<f64>>::new();
    for line in lines {
        let line = line.split_whitespace().collect::<Vec<_>>();
        let protein = line[0].to_owned();
        if !eps.contains(&protein) {
            continue;
        }
        let exp = line[3..]
            .into_iter()
            .map(|c| c.parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        protein_express.insert(protein, exp);
    }

    protein_express
}

// 计算均值和方差
fn get_mean_variance(data: &[f64]) -> (f64, f64) {
    assert!(data.len() != 0);

    let sum: f64 = data.iter().sum();
    let mean = sum / data.len() as f64;

    let variance: f64 = data
        .iter()
        .map(|value| {
            let diff = value - mean;
            diff * diff
        })
        .sum::<f64>()
        / data.len() as f64;

    (mean, variance)
}

#[cfg(test)]
mod tests {
    use crate::graph::Graph;

    use super::get_dpins;

    #[test]
    fn test_build_dpins() {
        let graph = Graph::new_from_file("./data/collins/collins.txt", true);

        let dpins = get_dpins(&graph);
        for dpin in dpins {
            print!(
                "node count: {}, edge cout: {}",
                dpin.node_count, dpin.edge_count
            );
            let edge = dpin.nei_list.iter().map(|c| c.len()).sum::<usize>();
            println!("\t{}", edge / 2);
        }
    }
}
