use clustering::kmeans;
use core::panic;
use once_cell::sync::Lazy;
use std::{
    cmp,
    collections::{BTreeSet, HashMap, HashSet, VecDeque},
    vec,
};

use super::graph::Graph;
use crate::complex::union::UnionFind;

pub const MAX_SIZE: usize = 20;
pub const MAX_ITER: usize = 10000;

#[derive(Clone, PartialEq)]
pub struct Complex {
    pub proteins: BTreeSet<String>,
    pub cohesion: f32,
}

impl PartialOrd for Complex {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        self.cohesion.partial_cmp(&other.cohesion)
    }
}

impl Complex {
    fn similarity(&self, oth: &Self) -> f32 {
        let a = self.proteins.iter().collect::<HashSet<_>>();
        let b = oth.proteins.iter().collect::<HashSet<_>>();
        let count = a.intersection(&b).count();

        count as f32 / std::cmp::max(a.len(), b.len()) as f32
    }
}

fn calculate_complex_cohesion(graph: &Graph, proteins: &BTreeSet<String>) -> f32 {
    // 示例：简单返回内部边的密度（你可以替换为真实逻辑）
    let indices: HashSet<usize> = proteins
        .iter()
        .filter_map(|name| graph.protein_id.get(name).cloned())
        .collect();

    let mut sum = HashMap::<usize, f32>::new();
    let mut count = HashMap::<usize, usize>::new();

    graph.edges.iter().into_iter().for_each(|(a, b, w)| {
        if indices.contains(a) && indices.contains(b) {
            count.entry(*a).and_modify(|c| *c += 1).or_insert(1);
            count.entry(*b).and_modify(|c| *c += 1).or_insert(1);

            sum.entry(*a).and_modify(|c| *c += *w).or_insert(*w);
            sum.entry(*b).and_modify(|c| *c += w).or_insert(*w);
        }
    });

    let mut cohesion = 0.;
    indices.iter().for_each(|p| {
        cohesion += sum[p] * (count[p] as f32 + 1.) / indices.len() as f32;
    });

    cohesion / indices.len() as f32
}

fn update_result(new_complex: Complex, result: &mut Vec<Complex>, threshold: f32) {
    // 示例：只要 cohesion 超过阈值就保留
    for complex in result.iter() {
        match new_complex
            .similarity(&complex)
            .partial_cmp(&threshold)
            .unwrap()
        {
            cmp::Ordering::Less => continue,
            _ => return,
        };
    }

    result.push(new_complex);
}

pub fn split_graph(queue: &mut VecDeque<Graph>, graphs: &mut Vec<Graph>) {
    println!("split graph");
    let current = queue.pop_front().unwrap();

    let mut union = current.get_union_finder();
    // 图不联通，则获取联通的子图
    if !union.is_connected() {
        // println!("not connected");
        current.get_components().into_iter().for_each(|nodes| {
            let sg_nodes = nodes
                .into_iter()
                .map(|n| current.id_protein[n].clone())
                .collect::<Vec<_>>();
            queue.push_back(current.subgraph(sg_nodes));
        });
        return;
    }

    if current.id_protein.len() <= MAX_SIZE {
        graphs.push(current);
        return;
    }

    current.split().into_iter().for_each(|g| {
        queue.push_back(g);
    });
}

pub fn get_complexes(graph: &Graph, result: &mut Vec<Complex>, similarity_threshold: f32) {
    println!("get_complexes ... graph {}", graph.protein_id.len());

    let n = graph.id_protein.len();
    let mut connectivity = vec![vec![false; n]; n];
    let mut complex_record: HashSet<u32> = HashSet::new();
    let mut queue: VecDeque<u32> = VecDeque::new();

    // 构建邻接矩阵
    for &(a, b, _) in &graph.edges {
        connectivity[a][b] = true;
        connectivity[b][a] = true;
    }

    // 初始化单点状态
    for i in 0..n {
        let state = 1u32 << i;
        queue.push_back(state);
        complex_record.insert(state);
    }

    // BFS 枚举所有连通子集
    while let Some(state) = queue.pop_front() {
        for i in 0..n {
            if state & (1 << i) != 0 {
                for j in 0..n {
                    if !connectivity[i][j] || (state & (1 << j)) != 0 {
                        continue;
                    }
                    let new_state = state | (1 << j);
                    if complex_record.contains(&new_state) {
                        continue;
                    }
                    complex_record.insert(new_state);
                    queue.push_back(new_state);
                }
            }
        }
    }

    // 生成实际蛋白质复合体
    // let mut complex_result = BTreeSet::new();
    for &state in &complex_record {
        let mut proteins = BTreeSet::new();
        for i in 0..n {
            if state & (1 << i) != 0 {
                proteins.insert(graph.id_protein[i].clone());
            }
        }

        // if proteins.len() >= 3 {
        //     let cohesion = calculate_complex_cohesion(graph, &proteins);
        //     complex_result.insert(Complex { proteins, cohesion });
        // }
    }

    // 降序排序并筛选加入结果集
    // complex_result.sort_by(|a, b| b.cohesion.partial_cmp(&a.cohesion).unwrap());
    // for i in complex_result[0..complex_result.len() / 100].iter() {
    //     result.push(i.clone());
    // }
    // if complex_result.is_empty() {
    //     return;
    // }
    // let mut res = vec![0];
    // for i in 1..complex_result.len() {
    //     let mut bool= true;

    //     for j  in res.iter() {
    //         match complex_result[i].similarity(&complex_result[*j]).partial_cmp(&similarity_threshold).unwrap() {
    //             cmp::Ordering::Less => continue,
    //             _ => {
    //                 bool = false;
    //                 break;
    //             }
    //         }
    //     }

    //     if bool {
    //         res.push(i);
    //     }
    // }

    // for idx in res {
    //     result.push(complex_result[idx].clone());
    // }
}

impl Graph {
    // todo 尝试谱聚类，或者直接适用于图的k-means 聚类算法
    pub fn split(&self) -> Vec<Graph> {
        let mut vecs: Vec<Vec<f32>> = vec![];

        self.id_protein.iter().for_each(|n| {
            vecs.push(self.ebd.get(n).unwrap().clone());
        });

        // 默认二分类
        let clustering: clustering::Clustering<'_, Vec<f32>> = kmeans(2, &vecs, MAX_ITER);
        // protein_id -->> cluster_id
        let membership = clustering.membership;

        let (mut nodes1, mut nodes2) = (Vec::new(), Vec::new());
        membership
            .into_iter()
            .enumerate()
            .for_each(|(protein, m)| match m {
                0 => nodes1.push(self.id_protein[protein].clone()),
                1 => nodes2.push(self.id_protein[protein].clone()),
                _ => panic!("kmeans error!"),
            });

        vec![self.subgraph(nodes1), self.subgraph(nodes2)]
    }

    pub fn get_union_finder(&self) -> UnionFind {
        let edges = self
            .edges
            .iter()
            .map(|(a, b, _)| (*a, *b))
            .collect::<Vec<_>>();
        UnionFind::new_from(self.id_protein.len(), edges)
    }

    pub fn get_components(&self) -> Vec<Vec<usize>> {
        let edges = self
            .edges
            .iter()
            .map(|(a, b, _)| (*a, *b))
            .collect::<Vec<_>>();
        let mut union = UnionFind::new_from(self.id_protein.len(), edges);
        union.get_components()
    }
}
