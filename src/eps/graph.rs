use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    vec,
};

use ndarray::{Array1, Array2};

use crate::eps::{
    errors::Result,
    options::Options,
    utils::{get_all_go_terms, read_clique, read_label, read_ppi_file, read_protein_go, BiMap},
};

/// 使用随机游走，因此采用邻接表存储
pub struct Graph {
    pub p_count: usize,
    pub label: Vec<bool>, // 需要参与最终的验证，必需要拿到
    pub essential_proteins: HashSet<usize>,
    pub ppg: Array2<f64>, // 一阶图
    pub pcg: Array2<f64>, // 一阶高阶关联图
    pub cpg: Array2<f64>,
    pub ccg: Array2<f64>, // 高阶图
    pub p2g: Array2<f64>, // 蛋白质到go的关联图
    pub g2p: Array2<f64>, // go term到go term的关联图
}

impl Graph {
    pub fn new(pcount: usize, ccount: usize, gcount: usize) -> Self {
        Self {
            p_count: 0,
            label: Vec::new(),
            essential_proteins: HashSet::new(),
            ppg: Array2::<f64>::zeros((pcount, pcount)),
            pcg: Array2::<f64>::zeros((pcount, ccount)),
            cpg: Array2::<f64>::zeros((ccount, pcount)),
            ccg: Array2::<f64>::zeros((ccount, ccount)),
            p2g: Array2::<f64>::zeros((pcount, gcount)),
            g2p: Array2::<f64>::zeros((gcount, gcount)),
        }
    }

    // 从各文件中解析出多阶网络关联图
    pub fn new_with(options: Options) -> Result<Self> {
        // 读取PPI文件，拿到蛋白质对应的id， 一阶网络的邻居
        let (pid, pnei) = read_ppi_file(&options.ppi_file).unwrap();
        let c = get_matrix((pid.len(), pid.len()), &pnei);

        let (cid, pcn) = read_clique(&options.clique_file, &pid, options.min_size).unwrap();
        // 这个函数目前可以处理低阶节点到高阶节点的游走概率
        let u = get_matrix((pid.len(), cid.len()), &pcn);
        let d = get_matrix_d((cid.len(), pid.len()), &cid, &pnei);
        // 高阶游走概率
        let gid = &get_all_go_terms();
        // 蛋白质到go的映射
        let (ptg, gtp) = read_protein_go(&pid, gid)?;
        let (p2g, g2p) = get_go_protein_matrix(&ptg, &gtp, pid.len(), 26653);

        let label = read_label(&options.label_file, &pid).unwrap();
        let mut essential_proteins = HashSet::new();
        label.iter().enumerate().for_each(|(p, f)| {
            if *f {
                essential_proteins.insert(p);
            }
        });

        Ok(Self {
            p_count: pid.len(),
            label,
            essential_proteins,
            ppg: c,
            pcg: u,
            cpg: d,
            ccg: Array2::default((1, 1)),
            p2g,
            g2p,
        })
    }

    /// 返回重要性排序，越重要的越靠前
    pub fn essentiality_rank(&self, scores: &Array1<f64>) -> Vec<usize> {
        let mut result: Vec<usize> = (0..self.label.len()).into_iter().collect();
        result.sort_by(|a, b| scores[*a].partial_cmp(&scores[*b]).unwrap());
        result.reverse();
        result
    }

    /// 生成待验证的结果
    pub fn get_result(&self, scores: &Array1<f64>) -> Vec<bool> {
        // 全部验证
        let count = self.label.iter().filter(|c| c == &&true).count();
        let aa = self.essentiality_rank(scores);
        let mut res = vec![false; self.label.len()];
        for i in 0..count {
            res[aa[i]] = true;
        }
        res
    }
}

// go-terms数量为26653
#[allow(unused)]
fn get_go_protein_matrix(
    ptg: &BTreeMap<usize, BTreeSet<usize>>,
    gtp: &BTreeMap<usize, BTreeSet<usize>>,
    p_count: usize,
    g_count: usize,
) -> (Array2<f64>, Array2<f64>) {
    let mut ptg_m = Array2::zeros((p_count, 26653));
    let mut gtp_m = Array2::zeros((26653, p_count));

    for (p, gos) in ptg.into_iter() {
        let score = 1. / gos.len() as f64;
        for go in gos.into_iter() {
            ptg_m[[*p, *go]] = score;
        }
    }

    for (go, proteins) in gtp.into_iter() {
        let score = 1. / proteins.len() as f64;
        for protein in proteins.into_iter() {
            gtp_m[[*go, *protein]] = score;
        }
    }

    (ptg_m, gtp_m)
}

// walkdown matrix
// 向下游走矩阵，根据蛋白质在复合物中的重要性程度随机游走
fn get_matrix_d(
    dim: (usize, usize),
    cid: &BiMap<BTreeSet<usize>>,
    pnei: &BTreeMap<usize, BTreeSet<usize>>,
) -> Array2<f64> {
    let mut m = Array2::<f64>::zeros(dim);
    let cid_vec = cid.into_vec();
    for (id, c) in cid_vec.into_iter().enumerate() {
        // for (id, c) in cid.into
        let mut score = HashMap::<usize, f64>::new();
        let mut sum = 0.;
        for p in c.iter() {
            if let Some(nei) = pnei.get(&p) {
                let com_len = nei.intersection(&c).collect::<Vec<_>>().len();
                sum += com_len as f64;
                score.insert(*p, com_len as f64);
            } else {
                score.insert(*p, 0.);
            }
        }
        // 更新转移概率
        for p in c.iter() {
            m[[id, *p]] = score[p] / sum;
        }
    }

    m
}

fn get_matrix((a, b): (usize, usize), nei: &BTreeMap<usize, BTreeSet<usize>>) -> Array2<f64> {
    let mut m = Array2::<f64>::zeros((a, b));
    nei.into_iter().for_each(|(a, b)| {
        let p = 1.0 / b.len() as f64;
        b.into_iter().for_each(|j| {
            m[[*a, *j]] = p;
        });
    });
    m
}

// 得到转移矩阵
#[allow(unused)]
fn parse_martix(
    dim: (usize, usize),
    edge_list: &Vec<(usize, usize)>,
    directed: bool,
) -> Array2<f64> {
    let mut matrix = Array2::<f64>::zeros(dim);
    match directed {
        false => edge_list.into_iter().for_each(|(a, b)| {
            matrix[[*a, *b]] = 1.;
            matrix[[*b, *a]] = 1.
        }),
        true => edge_list
            .into_iter()
            .for_each(|(a, b)| matrix[[*a, *b]] = 1.),
    }

    matrix
}

#[cfg(test)]
mod tests {
    use std::{
        collections::{BTreeMap, BTreeSet},
        vec,
    };

    use ndarray::{arr1, Array1, Array2};

    use super::{get_matrix, parse_martix, Graph};

    #[test]
    fn parse_matrix() {
        let edge: Vec<(usize, usize)> = vec![(1, 2), (1, 0), (0, 2)];
        let matrix = parse_martix((3, 3), &edge, false);
        assert_eq!(
            matrix,
            Array2::<f64>::ones((3, 3)) - Array2::<f64>::from_diag(&arr1(&[1f64, 1., 1.]))
        );

        let matrix = parse_martix((3, 3), &edge, true);
        let mut res = Array2::<f64>::zeros((3, 3));
        res[[1, 2]] = 1.;
        res[[1, 0]] = 1.;
        res[[0, 2]] = 1.;
        assert_eq!(matrix, res);
    }

    #[test]
    fn get_matrix_test() {
        let pnei = BTreeMap::from([
            (0, BTreeSet::from([1usize, 2])),
            (1, BTreeSet::from([0])),
            (2, BTreeSet::from([0])),
        ]);
        let m = get_matrix((3, 3), &pnei);
        dbg!("{}", m);
    }

    #[test]
    fn rank_teset() {
        let mut g = Graph::new(4, 3, 5);
        g.label = vec![false; 4];
        let score = Array1::<f64>::from_vec(vec![0.3, 0.2, 0.4, 0.1]);
        // println!("{}", score);
        let result = g.essentiality_rank(&score);
        // println!("{:?}", result);
        assert_eq!(vec![2, 0, 1, 3], result);
    }
}
