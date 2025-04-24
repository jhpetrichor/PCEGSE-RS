use crate::graph::Graph;
use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
};

pub const IS_A_WEIGHT: f64 = 0.8;
pub const PART_OF_WEIGHT: f64 = 0.6;

#[derive(Debug)]
pub struct Dag {
    pub edges: Vec<HashMap<usize, f64>>,
    pub protein_go: HashMap<String, HashSet<usize>>,
    pub go_child: HashMap<usize, HashSet<usize>>,

    pub(crate) sim_term: HashMap<(usize, usize), f64>,
    pub(crate) sim_term_child: HashMap<(usize, usize), f64>,
}

impl Dag {
    pub fn new() -> Self {
        let (edges, protein_go, go_child) = read_go_file();
        let mut dag = Self {
            edges: vec![Default::default(); edges.len()],
            protein_go,
            sim_term: Default::default(),
            go_child: Default::default(),
            sim_term_child: Default::default(),
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

    fn get_sim(&mut self, a: usize, b: usize) -> f64 {
        // 先查找
        match self.sim_term.get(&(a, b)) {
            Some(sim) => return *sim,
            None => {}
        }

        let a_ancestor = self.get_all_ancestors(a);
        let b_ancestor = self.get_all_ancestors(b);
        let common_ances = b_ancestor.intersection(&a_ancestor).collect::<Vec<_>>();
        // 在一个DAG中
        // a是b的祖先
        let sim = if b_ancestor.contains(&a) {
            let sv = self.get_semantic_value(b);
            let sum = sv.values().sum::<f64>();
            let mut s = 0.;
            common_ances.into_iter().for_each(|c| s += sv[c]);
            // println!("if: {}", (s + sv[&b]) / sum );
            (s + sv[&b]) / sum
        } else if b_ancestor.contains(&b) {
            // b 是 a 的祖先
            let sv = self.get_semantic_value(a);
            let sum = sv.values().sum::<f64>();
            let mut s = 0.;
            common_ances.into_iter().for_each(|c| s += sv[c]);
            // println!("if else: {}", (s + sv[&a]) / sum);
            (s + sv[&a]) / sum
        } else {
            let sv_a = self.get_semantic_value(a);
            let sv_b = self.get_semantic_value(b);
            let sum = sv_a.values().sum::<f64>() + sv_b.values().sum::<f64>();
            let mut s = 0.;
            common_ances
                .into_iter()
                .for_each(|c| s += f64::max(sv_a[&a], sv_b[&b]));
            // println!("else: {}", s / sum);
            s / sum
        };

        // 保存，避免多次计算
        self.sim_term.insert((a, b), sim);
        self.sim_term.insert((b, a), sim);

        sim
    }

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

    // delete edge
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
    use std::collections::{HashMap, HashSet};

    use crate::graph::Graph;

    use super::Dag;

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
        let mut dag = Dag {
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

        let res0 = dag.get_all_ancestors(0);
        assert_eq!(res0, HashSet::from([1, 2, 4, 3]));

        let res1 = dag.get_all_ancestors(1);
        assert_eq!(res1, HashSet::from([3, 4]));

        let res4 = dag.get_all_ancestors(4);
        assert!(res4.is_empty());

        let sv = dag.get_semantic_value(0);
        println!("{:?}", sv);
        let sv = dag.get_semantic_value(2);
        println!("{:?}", sv);

        let sim = dag.get_sim(1, 0);
        println!("{:?}", sim);

        let sim = dag.get_sim(1, 2);
        println!("ances: {:?}", sim);
        let sim = dag.get_sim_child(1, 2);
        println!("child: {:?}", sim);
        // let sim_ances_child = dag.get

        // function sim
        let s = dag.get_function_sim("a", "b");
        println!("sim of a and b is: {}", s);
        let s = dag.get_sim_ancestor_child("a", "b");
        println!("sim of a and b is: {}", s);
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
