use std::collections::{HashMap, HashSet};

use ndarray_rand::rand::random;

use crate::{
    eva::{update_by_cohesion, Complex}, gene_expression::read_essential_protein, graph::Graph
};

impl Graph {
    pub fn all_cc(&self) -> Vec<f64> {
        (0..self.node_count)
            .map(|node| self.clustering_coefficient(node))
            .collect()
    }

    pub fn clustering_coefficient(&self, node: usize) -> f64 {
        // 获取节点的邻居
        let neighbors = &self.nei_list[node];
        let neighbor_count = neighbors.len();

        // 如果节点没有邻居，聚类系数为 0
        if neighbor_count < 2 {
            return 0.0;
        }

        let mut actual_edges = 0;

        // 遍历每对邻居，检查它们之间是否有边
        for (i, _) in neighbors {
            for (j, _) in neighbors {
                if i != j && self.nei_list[*i].contains_key(&j) {
                    actual_edges += 1;
                }
            }
        }

        // 每对邻居之间的最大可能边数
        let max_possible_edges = neighbor_count * (neighbor_count - 1) / 2;

        // 聚类系数公式
        actual_edges as f64 / (max_possible_edges as f64 * 2.0)
    }

    pub fn cc_sort(&self) -> Vec<usize> {
        let cc = self.all_cc();
        let mut seed = (0..self.node_count).collect::<Vec<_>>();
        seed.sort_by(|a, b| cc[*a].partial_cmp(&cc[*b]).unwrap());
        seed.reverse();
        seed
    }
}

pub fn pcegs(graph: &Graph, bate: f64) -> Vec<Complex<String>> {
    // 计算节点的权重
    let node_weight = graph
        .nei_list
        .iter()
        .map(|c| c.values().sum::<f64>())
        .collect::<Vec<_>>();
    // 节点对之间的权重
    let mut a = HashMap::<(usize, usize), f64>::new();
    graph.nei_list.iter().enumerate().for_each(|(n, nei)| {
        for (v, w) in nei.iter() {
            let value = node_weight[n] * node_weight[*v] / (1.0 + (1.0 + w).log10()).powi(2);
            a.insert((n, *v), value);
            a.insert((*v, n), value);
        }
    });

    // 节点的影响力
    let mut node_influ = HashMap::new();
    graph.nei_list.iter().enumerate().for_each(|(n, nei)| {
        let att = nei
            .iter()
            .map(|(m, _)| a.get(&(*m, n)).unwrap())
            .sum::<f64>();
        node_influ.insert(n, att);
    });

    // step 1. 获取种子节点
    let seeds = graph.cc_sort();
    let mut visited = HashSet::<usize>::new();

    let mut complexes = Vec::new();
    for seed in seeds {
        if visited.contains(&seed) {
            continue;
        }
        // 核心
        let mut core = graph.nei_list[seed]
            .keys()
            .map(|c| c.to_owned())
            .collect::<HashSet<usize>>();
        // 候选附属,即核心对应的邻居
        let mut all_neis = HashSet::<usize>::new();
        core.iter().for_each(|c| {
            if let Some(ns) = graph.nei_list.get(*c) {
                ns.keys().for_each(|m| {
                    if !core.contains(m) {
                        all_neis.insert(*m);
                    }
                });
            }
        });
        core.insert(seed);
        // nei.insert(seed);
        // let core = nei;
        // 更新访问过的节点
        core.iter().for_each(|n| {
            visited.insert(*n);
        });
        // 根据情况是否将其加入到
        let mut real_attach = HashSet::<usize>::new();
        // 计算社区对节点的吸引力
        for n in all_neis {
            let sum_a = core
                .iter()
                .map(|d| match a.get(&(*d, n)) {
                    Some(v) => *v,
                    None => 0.,
                })
                .sum::<f64>();

            match (sum_a / node_influ[&n]).ge(&bate) {
                true => {
                    real_attach.insert(n);
                }
                false => {}
            };
        }
        let proteins = core
            .union(&real_attach)
            .map(|a| a.clone())
            .collect::<HashSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        if proteins.len() >= 3 {
            let cohesion = graph.calculate_cohesion(&proteins);
            let complex = Complex { proteins, cohesion };
            complexes.push(complex);
        }
    }
    let complexes_id = update_by_cohesion(complexes);
    complexes_id
        .into_iter()
        .map(|c| {
            let proteins = c
                .proteins
                .into_iter()
                .map(|id| graph.id_protein[id].to_string())
                .collect::<Vec<_>>();
            Complex {
                proteins,
                cohesion: c.cohesion,
            }
        })
        .collect::<Vec<_>>()
}

pub fn pcegs_essential(graph: &Graph, bate: f64) -> Vec<Complex<String>> {
    // 计算节点的权重
    let node_weight = graph
        .nei_list
        .iter()
        .map(|c| c.values().sum::<f64>())
        .collect::<Vec<_>>();
    // 节点对之间的权重
    let mut a = HashMap::<(usize, usize), f64>::new();
    graph.nei_list.iter().enumerate().for_each(|(n, nei)| {
        for (v, w) in nei.iter() {
            let value = node_weight[n] * node_weight[*v] / (1.0 + (1.0 + w).log10()).powi(2);
            a.insert((n, *v), value);
            a.insert((*v, n), value);
        }
    });

    // 节点的影响力
    let mut node_influ = HashMap::new();
    graph.nei_list.iter().enumerate().for_each(|(n, nei)| {
        let att = nei
            .iter()
            .map(|(m, _)| a.get(&(*m, n)).unwrap())
            .sum::<f64>();
        node_influ.insert(n, att);
    });

    // step 1. 获取种子节点
    let seeds = graph.cc_sort();
    let mut visited = HashSet::<usize>::new();

    let eps = read_essential_protein();
    let eps = (0..graph.id_protein.len() - 1).into_iter().filter(|i| eps.contains(&graph.id_protein[*i])).collect::<HashSet<_>>();

    let mut complexes = Vec::new();
    for seed in seeds {
        if visited.contains(&seed) || !eps.contains(&seed) {
            continue;
        }
        // 核心
        let mut core = graph.nei_list[seed].keys().map(|a| a.clone()).collect::<HashSet<_>>().intersection(&eps).map(|a| a.clone()).collect::<HashSet<_>>();

        // 候选附属,即核心对应的邻居
        let mut all_neis = HashSet::<usize>::new();
        core.iter().for_each(|c| {
            if let Some(ns) = graph.nei_list.get(*c) {
                ns.keys().for_each(|m| {
                    if !core.contains(m) {
                        all_neis.insert(*m);
                    }
                });
            }
        });
        core.insert(seed);
        // nei.insert(seed);
        // let core = nei;
        // 更新访问过的节点
        core.iter().for_each(|n| {
            visited.insert(*n);
        });
        // 根据情况是否将其加入到
        let mut real_attach = HashSet::<usize>::new();
        // 计算社区对节点的吸引力
        for n in all_neis {
            let sum_a = core
                .iter()
                .map(|d| match a.get(&(*d, n)) {
                    Some(v) => *v,
                    None => 0.,
                })
                .sum::<f64>();

            match (sum_a / node_influ[&n]).ge(&bate) {
                true => {
                    real_attach.insert(n);
                }
                false => {}
            };
        }
        let proteins = core
            .union(&real_attach)
            .map(|a| a.clone())
            .collect::<HashSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        if proteins.len() >= 3 {
            let cohesion = graph.calculate_cohesion(&proteins);
            let complex = Complex { proteins, cohesion };
            complexes.push(complex);
        }
    }
    let complexes_id = update_by_cohesion(complexes);
    complexes_id
        .into_iter()
        .map(|c| {
            let proteins = c
                .proteins
                .into_iter()
                .map(|id| graph.id_protein[id].to_string())
                .collect::<Vec<_>>();
            Complex {
                proteins,
                cohesion: c.cohesion,
            }
        })
        .collect::<Vec<_>>()
}


#[cfg(test)]
mod tests {
    use crate::graph::Graph;

    use super::pcegs;

    #[test]
    fn test_cc() {
        let mut g = Graph::new_from(vec![
            (0, 1, 0.1),
            (0, 2, 0.1),
            (0, 3, 0.1),
            (1, 2, 0.1),
            (1, 3, 0.1),
            (2, 3, 0.1),
        ]);
        println!("{}", g.clustering_coefficient(0));
        println!("{:?}", g.all_cc());

        g.remove_edge(1, 2);
        g.remove_edge(1, 3);
        println!("{}", g.clustering_coefficient(0));
        println!("{:?}", g.all_cc());

        let seed = g.cc_sort();
        println!("{:?}", seed)
    }

    #[test]
    fn test_pcegs() {
        let graph = Graph::new_from_file("./data/collins/collins.txt", true);
        pcegs(&graph, 0.4);
    }
    
}
