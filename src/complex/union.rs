/// 并查集实现，以遍历Graph的连通分量
/// 帮助实现Graph的分裂操作
use std::collections::HashMap;

pub struct UnionFind {
    parent: Vec<usize>,
    count: usize,
}

impl UnionFind {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            count: n,
        }
    }

    pub fn new_from(n: usize, pairs: Vec<(usize, usize)>) -> Self {
        let mut union_find = Self::new(n);
        pairs.into_iter().for_each(|(a, b)| {
            union_find.union(a, b);
        });

        union_find
    }

    pub fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    pub fn union(&mut self, x: usize, y: usize) {
        let fx = self.find(x);
        let fy = self.find(y);
        if fx != fy {
            self.parent[fx] = fy;
            self.count -= 1;
        }
    }

    pub fn is_connected(&self) -> bool {
        self.count == 1
    }

    /// 获取所有连通分量
    pub fn get_components(&mut self) -> Vec<Vec<usize>> {
        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
        let n = self.parent.len();
        for i in 0..n {
            let root = self.find(i); // 会做路径压缩
            groups.entry(root).or_default().push(i);
        }
        groups.into_values().collect()
    }
}

#[test]
fn main() {
    let n = 6;
    let edges = vec![(0, 1), (1, 2), (3, 4)];

    let mut uf = UnionFind::new(n);
    for (u, v) in edges {
        uf.union(u, v);
    }

    let components = uf.get_components();
    println!("Connected components:");
    for (i, comp) in components.iter().enumerate() {
        println!("Component {}: {:?}", i + 1, comp);
    }
}
