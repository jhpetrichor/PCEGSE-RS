use crate::graph::Graph;
use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
    fs::read_to_string,
};

pub const IS_A_WEIGHT: f64 = 0.8;
pub const PART_OF_WEIGHT: f64 = 0.6;

#[derive(Debug)]
pub struct Dag {
    edges: Vec<HashMap<usize, f64>>,
    protein_go: HashMap<String, HashSet<usize>>,
    go_child: HashMap<usize, HashSet<usize>>,

    // sim_term只用作计算使用
    sim_term: HashMap<(usize, usize), f64>, // 计算时候更新
    sim_term_child: HashMap<(usize, usize), f64>, // 计算时更新
}

impl Dag {
    pub fn new() -> Self {
        let (edges, protein_go, go_child) = read_go_file();
        let mut dag = Self {
            edges: vec![Default::default(); edges.len()],
            protein_go,
            sim_term: Default::default(), // 计算是更新
            go_child,
            sim_term_child: Default::default(), // 计算时更新
        };

        edges.into_iter().for_each(|(a, b, w)| {
            dag.add_edge(a, b, w);
        });

        dag
    }

    // 根据父项和子项计算蛋白质的功能相似性
    pub fn get_sim_ancestor_child(&mut self, p1: &str, p2: &str) -> f64 {
        (self.get_function_sim(p1, p2) + self.get_function_sim_child(p1, p2)) / 2.0
    }

    // 获取蛋白质的功能相似性（仅仅根据父项术语）
    pub fn get_function_sim(&mut self, p1: &str, p2: &str) -> f64 {
        // 获取蛋白质对应的go术语, 如果没有
        let gos1 = self.protein_go.get(p1);
        let gos2 = self.protein_go.get(p2);
        // 没有共同术语则为0.
        if gos1.is_none() || gos2.is_none() {
            return 0.;
        }

        let gos1 = gos1.unwrap().clone();
        let gos2 = gos2.unwrap().clone();
        if gos1.is_empty() || gos2.is_empty() {
            return 0.;
        }

        // 有共同术语，则需要计算两两之间的最大值
        let mut sum1 = 0.;
        for go1 in gos1.iter() {
            let mut max_sim = 0.;
            for go2 in gos2.iter() {
                let sim = self.get_sim(*go1, *go2);
                max_sim = f64::max(sim, max_sim);
            }
            sum1 += max_sim;
        }

        let mut sum2 = 0.;
        for go1 in gos2.iter() {
            let mut max_sim = 0.;
            for go2 in gos1.iter() {
                let sim = self.get_sim(*go1, *go2);
                max_sim = f64::max(sim, max_sim);
            }
            sum2 += max_sim;
        }

        (sum1 + sum2) / (gos1.len() + gos2.len()) as f64
    }

    // 获取蛋白质的功能相似性（仅仅根据子项术语）
    pub fn get_function_sim_child(&mut self, p1: &str, p2: &str) -> f64 {
        // 获取蛋白质对应的go术语, 如果没有
        let gos1 = self.protein_go.get(p1).unwrap().clone();
        let gos2 = self.protein_go.get(p2).unwrap().clone();
        // 没有共同术语则为0.
        if gos1.intersection(&gos2).count() == 0 {
            return 0.;
        }

        // 有共同术语，则需要计算两两之间的最大值
        let mut sum1 = 0.;
        for go1 in gos1.iter() {
            let mut max_sim = 0.;
            for go2 in gos2.iter() {
                let sim = self.get_sim_child(*go1, *go2);
                max_sim = f64::max(sim, max_sim);
            }
            sum1 += max_sim;
        }

        let mut sum2 = 0.;
        for go1 in gos2.iter() {
            let mut max_sim = 0.;
            for go2 in gos1.iter() {
                let sim = self.get_sim_child(*go1, *go2);
                max_sim = f64::max(sim, max_sim);
            }
            sum2 += max_sim;
        }

        (sum1 + sum2) / (gos1.len() + gos2.len()) as f64
    }

    pub fn add_edge(&mut self, a: usize, b: usize, r: f64) {
        self.edges[a].insert(b, r);
    }

    /// 获取某个节点的语义值（包括祖先节点的语义贡献）
    fn calculate_semantic_value(&self, node: usize, memo: &mut HashMap<usize, f64>) {
        if let Some(parents) = self.edges.get(node) {
            for (&parent, w) in parents.iter() {
                // 是否已经计算过？如果计算过是否要更新
                let sv = memo[&node] * w;
                // todo 可以添加条件剪枝
                memo.entry(parent)
                    .and_modify(|c| {
                        *c = f64::max(sv, *c);
                    })
                    .or_insert(sv);
                self.calculate_semantic_value(parent, memo);
            }
        }
    }

    // 获取所有祖先节点对自身的语义贡献值
    fn get_semantic_value(&mut self, t: usize) -> HashMap<usize, f64> {
        let mut res = HashMap::<usize, f64>::from([(t, 1.)]);
        self.calculate_semantic_value(t, &mut res);
        res
    }

    // A new method to measure the semantic similarity of GO terms
    // https://doi.org/10.1093/bioinformatics/btm087
    #[allow(unused)]
    fn get_sim_wang(&mut self, a: usize, b: usize) -> f64 {
        // 先查找
        match self.sim_term.get(&(a, b)) {
            Some(sim) => return *sim,
            None => {}
        }

        let a_ancestor = self.get_all_ancestors(a);
        let b_ancestor = self.get_all_ancestors(b);
        let common_ances = b_ancestor.intersection(&a_ancestor).collect::<Vec<_>>();
        // 在一个DAG中
        let sim = if common_ances.is_empty() {
            0.
        } else {
            // let sim = 0.;
            let sv_a = self.get_semantic_value(b);
            let sv_b = self.get_semantic_value(a);

            common_ances
                .into_iter()
                .map(|c| sv_a[c] + sv_b[c])
                .sum::<f64>()
                / (sv_a.values().sum::<f64>() + sv_b.values().sum::<f64>())
        };

        self.sim_term.insert((a, b), sim);
        self.sim_term.insert((b, a), sim);

        sim
    }

    fn get_sim(&mut self, a: usize, b: usize) -> f64 {
        // 先查找，如果计算过，则查找返回
        if let Some(sim) = self.sim_term.get(&(a, b)) {
            return *sim;
        }

        // let a_ancestor = self.get_all_ancestors(a);
        // let b_ancestor = self.get_all_ancestors(b);
        // let res = self.find_lca_with_paths(a, b);
        let sim = match self.find_lca_with_paths(a, b) {
            None => 0.,
            Some((_, a_shortest, b_shortest)) => {
                match a_shortest.len().cmp(&b_shortest.len()) {
                    // a 是 b的祖先，计算b的语义贡献值
                    Ordering::Less => {
                        let sv = self.get_semantic_value(b);
                        sv[&a_shortest.last().unwrap()]
                            / b_shortest.into_iter().map(|i| sv[&i]).sum::<f64>()
                    }
                    // b 是 a的祖先
                    Ordering::Greater => {
                        let sv = self.get_semantic_value(a);
                        sv[&b_shortest.last().unwrap()]
                            / a_shortest.into_iter().map(|i| sv[&i]).sum::<f64>()
                    }
                    // a是b的堂兄弟
                    Ordering::Equal => {
                        let sv_a = self.get_semantic_value(a);
                        let sv_b = self.get_semantic_value(b);
                        f64::max(
                            sv_a[&a_shortest.last().unwrap()],
                            sv_b[&b_shortest.last().unwrap()],
                        ) / f64::max(
                            a_shortest.into_iter().map(|i| sv_a[&i]).sum::<f64>(),
                            b_shortest.into_iter().map(|i| sv_b[&i]).sum::<f64>(),
                        )
                    }
                }
            }
        };

        // 保存，避免多次计算
        self.sim_term.insert((a, b), sim);
        self.sim_term.insert((b, a), sim);

        sim
    }

    // ! Previous Version
    // fn get_sim(&mut self, a: usize, b: usize) -> f64 {
    //     // 先查找
    //     match self.sim_term.get(&(a, b)) {
    //         Some(sim) => return *sim,
    //         None => {}
    //     }

    //     let a_ancestor = self.get_all_ancestors(a);
    //     let b_ancestor = self.get_all_ancestors(b);
    //     let common_ances = b_ancestor.intersection(&a_ancestor).collect::<Vec<_>>();
    //     // 在一个DAG中
    //     let sim = if common_ances.is_empty() {
    //         0.
    //     } else if b_ancestor.contains(&a) {
    //         let sv = self.get_semantic_value(b);
    //         let sum = sv.values().sum::<f64>();
    //         let mut s = 0.;
    //         common_ances.into_iter().for_each(|c| s += sv[c]);
    //         // println!("if: {}", (s + sv[&b]) / sum );
    //         (s + sv[&b]) / sum * 2.
    //     } else if b_ancestor.contains(&b) {
    //         // b 是 a 的祖先
    //         let sv = self.get_semantic_value(a);
    //         let sum = sv.values().sum::<f64>();
    //         let mut s = 0.;
    //         common_ances.into_iter().for_each(|c| s += sv[c]);
    //         // println!("if else: {}", (s + sv[&a]) / sum);
    //         (s + sv[&a]) / sum* 2.
    //     } else {
    //         let sv_a = self.get_semantic_value(a);
    //         let sv_b = self.get_semantic_value(b);
    //         let sum = sv_a.values().sum::<f64>() + sv_b.values().sum::<f64>();
    //         let mut s = 0.;
    //         common_ances
    //             .into_iter()
    //             .for_each(|c| s += f64::max(sv_a[c], sv_b[&c]));
    //         // println!("else: {}", s / sum);
    //         s / sum* 2.
    //     };

    //     // 保存，避免多次计算
    //     self.sim_term.insert((a, b), sim);
    //     self.sim_term.insert((b, a), sim);

    //     sim
    // }

    fn get_sim_child(&mut self, a: usize, b: usize) -> f64 {
        if let Some(sim) = self.sim_term_child.get(&(a, b)) {
            return *sim;
        }

        // let a_child = self.go_child.get(&a);
        // let b_child = self.go_child.get(&b);
        let a_child = self.get_all_child(a);
        let b_child = self.get_all_child(b);
        if a_child.is_empty() || b_child.is_empty() {
            return 0.0;
        }

        let unionmsize = a_child.union(&b_child).count() as f64;
        let commonsize = a_child.intersection(&b_child).count() as f64;

        let sim = unionmsize / commonsize;
        self.sim_term_child.insert((a, b), sim);
        self.sim_term_child.insert((b, a), sim);

        sim
    }

    fn get_all_child(&self, t: usize) -> HashSet<usize> {
        let mut visited = HashSet::new();
        self.dfs_child(t, &mut visited);
        visited
    }

    fn dfs_child(&self, node: usize, visited: &mut HashSet<usize>) {
        if let Some(childs) = self.go_child.get(&node) {
            for child in childs.iter() {
                if visited.insert(*child) {
                    self.dfs(*child, visited);
                }
            }
        }
    }

    fn get_all_ancestors(&self, t: usize) -> HashSet<usize> {
        let mut visited = HashSet::new();
        self.dfs(t, &mut visited);
        visited
    }

    /// 深度优先搜索，向上查找祖先
    fn dfs(&self, node: usize, visited: &mut HashSet<usize>) {
        if let Some(parents) = self.edges.get(node) {
            for (&parent, _) in parents.iter() {
                if visited.insert(parent) {
                    self.dfs(parent, visited);
                }
            }
        }
    }

    pub fn get_ancestors(&self, start: usize) -> HashSet<usize> {
        let mut visited = HashSet::new();
        let mut stack = vec![start];

        while let Some(node) = stack.pop() {
            if visited.insert(node) {
                if let Some(parents) = self.edges.get(node) {
                    for &parent in parents.keys() {
                        stack.push(parent);
                    }
                }
            }
        }

        visited
    }

    pub fn get_ancestor_distances(&self, start: usize) -> HashMap<usize, usize> {
        let mut distances = HashMap::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back((start, 0));

        while let Some((node, dist)) = queue.pop_front() {
            if distances.insert(node, dist).is_none() {
                if let Some(parents) = self.edges.get(node) {
                    for &parent in parents.keys() {
                        queue.push_back((parent, dist + 1));
                    }
                }
            }
        }

        distances
    }

    pub fn find_lca(&self, node1: usize, node2: usize) -> Option<usize> {
        let dist1 = self.get_ancestor_distances(node1);
        let dist2 = self.get_ancestor_distances(node2);

        let common: HashSet<_> = dist1
            .keys()
            .copied()
            .collect::<HashSet<_>>()
            .intersection(&dist2.keys().copied().collect())
            .copied()
            .collect();

        common
            .into_iter()
            .min_by_key(|&ancestor| dist1[&ancestor] + dist2[&ancestor])
    }

    /// 获取祖先距离和前驱信息，用于构造路径
    fn get_ancestor_paths(&self, start: usize) -> (HashMap<usize, usize>, HashMap<usize, usize>) {
        let mut distances = HashMap::new();
        let mut predecessors = HashMap::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back((start, 0));

        while let Some((node, dist)) = queue.pop_front() {
            if distances.insert(node, dist).is_none() {
                if let Some(parents) = self.edges.get(node) {
                    for &parent in parents.keys() {
                        if !distances.contains_key(&parent) {
                            predecessors.insert(parent, node);
                            queue.push_back((parent, dist + 1));
                        }
                    }
                }
            }
        }

        (distances, predecessors)
    }

    /// 构造从 start 到 end 的路径（子 → 祖先）
    fn build_path(
        &self,
        predecessors: &HashMap<usize, usize>,
        start: usize,
        end: usize,
    ) -> Vec<usize> {
        let mut path = vec![end];
        let mut current = end;
        while current != start {
            if let Some(&prev) = predecessors.get(&current) {
                path.push(prev);
                current = prev;
            } else {
                break; // 没有路径连通
            }
        }
        path.reverse();
        path
    }

    /// 返回最近公共祖先和从两个节点到该祖先的路径
    pub fn find_lca_with_paths(
        &self,
        node1: usize,
        node2: usize,
    ) -> Option<(usize, Vec<usize>, Vec<usize>)> {
        let (dist1, pred1) = self.get_ancestor_paths(node1);
        let (dist2, pred2) = self.get_ancestor_paths(node2);

        let common: HashSet<_> = dist1
            .keys()
            .copied()
            .collect::<HashSet<_>>()
            .intersection(&dist2.keys().copied().collect())
            .copied()
            .collect();

        let lca = common
            .into_iter()
            .min_by_key(|&ancestor| dist1[&ancestor] + dist2[&ancestor])?;

        let path1 = self.build_path(&pred1, node1, lca);
        let path2 = self.build_path(&pred2, node2, lca);

        Some((lca, path1, path2))
    }
}

pub fn weight_by_dag(graph: &mut Graph) {
    let mut dag = Dag::new();

    graph.nei_list.iter_mut().enumerate().for_each(|(a, nei)| {
        nei.iter_mut().for_each(|(b, w)| {
            *w = dag.get_function_sim(graph.id_protein[a].as_str(), graph.id_protein[*b].as_str());
        });
    })
}

pub fn weight_by_dag_topo(graph: &mut Graph, alpha: f64) {
    let mut dag = Dag::new();

    // let mut func_sim = HashMap::<(usize, usize), f64>::new();
    // 暂时存储拓扑相似性，以避开借用检查机制
    let mut topo_sim = HashMap::<(usize, usize), f64>::new();
    graph.nei_list.iter().enumerate().for_each(|(a, nei)| {
        nei.iter().for_each(|(b, _)| {
            let sim = graph.jaccard(a, *b);
            topo_sim.insert((a, *b), sim);
            topo_sim.insert((*b, a), sim);
        });
    });

    let mut edge_remove = Vec::new();

    // 功能相似性
    graph.nei_list.iter_mut().enumerate().for_each(|(a, nei)| {
        nei.iter_mut().for_each(|(b, w)| {
            let sim =
                dag.get_function_sim(graph.id_protein[a].as_str(), graph.id_protein[*b].as_str());
            if sim.le(&0.1) {
                edge_remove.push((a, *b));
            }

            let topo_sim = topo_sim.get(&(a, *b)).unwrap();
            *w = alpha * sim + (1. - alpha) * topo_sim;
        });
    });

    // 删除连接强度太低的边
    edge_remove.into_iter().for_each(|(a, b)| {
        graph.remove_edge(a, b);
    });
}

fn read_go_file() -> (
    Vec<(usize, usize, f64)>,
    HashMap<String, HashSet<usize>>,
    HashMap<usize, HashSet<usize>>,
) {
    let mut go_terms = Vec::<String>::new();
    let mut go_term_id = HashMap::<String, usize>::new();
    let mut edges = Vec::<(usize, usize, f64)>::new();
    let mut go_child: HashMap<usize, HashSet<usize>> = HashMap::new();

    // 读取 IS_A 关系
    let contents = read_to_string("./data/is_a.txt").expect("Failed to read IS_A file!");
    for line in contents.lines() {
        let line = line.split_whitespace().collect::<Vec<_>>();
        let line: Vec<usize> = line
            .into_iter()
            .map(|s| {
                if let Some(id) = go_term_id.get(s) {
                    *id
                } else {
                    go_term_id.insert(s.to_string(), go_term_id.len());
                    go_terms.push(s.to_string());
                    go_terms.len() - 1
                }
            })
            .collect();
        // 更新边
        for i in 1..line.len() {
            edges.push((line[0], line[i], IS_A_WEIGHT));
            go_child
                .entry(line[i])
                .and_modify(|c| {
                    c.insert(line[0]);
                })
                .or_insert(HashSet::from([line[0]]));
        }
    }

    // 读取 PART_OF 关系
    let contents = read_to_string("./data/part_of.txt").expect("Failed to read IS_A file!");
    for line in contents.lines() {
        let line = line.split_whitespace().collect::<Vec<_>>();
        let line: Vec<usize> = line
            .into_iter()
            .map(|s| {
                if let Some(id) = go_term_id.get(s) {
                    *id
                } else {
                    go_term_id.insert(s.to_string(), go_term_id.len());
                    go_terms.push(s.to_string());
                    go_terms.len() - 1
                }
            })
            .collect();
        // 更新边
        for i in 1..line.len() {
            edges.push((line[0], line[i], PART_OF_WEIGHT));
            go_child
                .entry(line[i])
                .and_modify(|c| {
                    c.insert(line[0]);
                })
                .or_insert(HashSet::from([line[0]]));
        }
    }

    // 读取蛋白质到go term的映射
    let mut protein_go = HashMap::<String, HashSet<usize>>::new();
    let contents = read_to_string("./data/go_slim.txt").expect("Failed to read go slim file!");
    for line in contents.lines() {
        let line = line.split_whitespace().collect::<Vec<_>>();
        let mut gos = HashSet::<usize>::new();
        line[1..].into_iter().for_each(|g| {
            if let Some(id) = go_term_id.get(*g) {
                gos.insert(*id);
            }
        });
        protein_go.insert(line[0].to_string(), gos);
    }

    return (edges, protein_go, go_child);
}

#[cfg(test)]
mod tests {
    use super::Dag;
    use crate::graph::Graph;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn test_dag_from_file() {
        let mut dag = Dag::new();
        println!("{}", dag.protein_go.len());

        println!("{}", dag.get_function_sim("YAL011W", "YBR231C"));
    }

    #[test]
    fn test_get_all_ancestor() {
        let a = vec![
            HashMap::from([(1, 0.6), (2, 0.8), (3, 0.8)]),
            HashMap::from([(3, 0.8), (4, 0.8)]),
            HashMap::from([(3, 0.6)]),
            HashMap::from([(4, 0.6)]),
            HashMap::new(),
        ];

        let protien_go = HashMap::from([
            ("a".to_string(), HashSet::from([1, 2])),
            ("b".to_string(), HashSet::from([3, 1, 4])),
        ]);
        let dag = Dag {
            edges: a,
            protein_go: protien_go,
            sim_term: Default::default(),
            go_child: HashMap::from([
                (4, HashSet::from([1, 3])),
                (3, HashSet::from([0, 1, 2])),
                (2, HashSet::from([0])),
                (1, HashSet::from([0])),
            ]),
            sim_term_child: Default::default(),
        };

        let anc0 = dag.get_ancestors(0);
        println!("0 ancestor: {:?}", anc0);
        let anc1 = dag.get_ancestors(1);
        println!("0 ancestor: {:?}", anc1);
        println!("dis 0: {:?}", dag.get_ancestor_distances(0));
        println!("dis 1: {:?}", dag.get_ancestor_distances(1));
        println!("LCA 0 1: {:?}", dag.find_lca(0, 1));

        println!("shortest path and lca {:?}", dag.find_lca_with_paths(1, 2));
    }

    #[test]
    fn test_graph_weight() {
        let graph = Graph::new_from_file("./data/collins/collins.txt", false);
        let mut dag = Dag::new();
        for (a, neis) in graph.nei_list.iter().enumerate() {
            for (b, _) in neis {
                println!(
                    "{}\t{}\t{:.4?}\t{:.4?}\t{:.4?}\t{:.4?}",
                    graph.id_protein[a],
                    graph.id_protein[*b],
                    dag.get_function_sim(
                        graph.id_protein[a].as_str(),
                        graph.id_protein[*b].as_str()
                    ),
                    dag.get_function_sim_child(
                        graph.id_protein[a].as_str(),
                        graph.id_protein[*b].as_str()
                    ),
                    dag.get_sim_ancestor_child(
                        graph.id_protein[a].as_str(),
                        graph.id_protein[*b].as_str()
                    ),
                    graph.jaccard(a, *b)
                );
            }
        }
    }
}
