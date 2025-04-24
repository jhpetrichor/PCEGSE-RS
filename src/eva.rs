use std::{
    collections::HashSet,
    fmt::{self, Display},
    fs::read_to_string,
    hash::Hash,
};

use crate::graph::Graph;

pub const THRESHOLD_OS: f64 = 0.2;
pub const COMPLEX_REF: &str = "./data/complex.txt";

pub const OVERLAP_SCORE: f64 = 0.6;

const MAX_SIZE: usize = 20;
const MIN_SIZE: usize = 3;

#[derive(Debug, Clone)]
pub struct Complex<T> {
    pub proteins: Vec<T>,
    pub cohesion: f64,
}

impl<T> Complex<T>
where
    T: Eq + Hash + Clone + Ord,
{
    pub fn new(proteins: Vec<T>, cohesion: f64) -> Self {
        Self { proteins, cohesion }
    }

    pub fn overlap_score(&self, other: &Self) -> (f64, f64) {
        let set_self: HashSet<_> = self.proteins.iter().collect();
        let set_other: HashSet<_> = other.proteins.iter().collect();
        let comsize = set_self.intersection(&set_other).count() as f64;

        let denom = (self.proteins.len() * other.proteins.len()) as f64;
        let os = if denom == 0.0 {
            0.0
        } else {
            comsize.powi(2) / denom
        };
        (comsize, os)
    }

    pub fn is_overlapped(&self, other: &Self) -> bool {
        let set_self: HashSet<_> = self.proteins.iter().collect();
        let set_other: HashSet<_> = other.proteins.iter().collect();
        let comsize = set_self.intersection(&set_other).count() as f64;

        let denom = usize::max(set_other.len(), set_self.len()) as f64;

        (comsize / denom).ge(&OVERLAP_SCORE)
    }
}

impl<T> Display for Complex<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {:.4}",
            self.proteins
                .iter()
                .map(|p| format!("{}", p))
                .collect::<Vec<String>>()
                .join("\t"),
            self.cohesion
        )
    }
}

impl Graph {
    pub fn calculate_cohesion(&self, cluster: &Vec<usize>) -> f64 {
        let mut cohesion = 0.;

        for i in cluster.iter() {
            let mut count = 0.;
            let mut weight = 0.;
            for j in cluster.iter() {
                if let Some(w) = self.nei_list[*i].get(j) {
                    count += 1.;
                    weight += *w;
                }
            }
            cohesion += weight * (count + 1.) / cluster.len() as f64;
        }

        cohesion / cluster.len() as f64
    }
}

pub fn update_by_cohesion<T>(mut complexes: Vec<Complex<T>>) -> Vec<Complex<T>>
where
    T: Eq + Hash + Clone + Ord,
{
    // 从大到小排序
    complexes.sort_by(|a, b| a.cohesion.partial_cmp(&b.cohesion).unwrap());
    complexes.reverse();

    let mut res = Vec::<Complex<T>>::new();
    complexes.into_iter().for_each(|c| {
        let mut flag = true;
        for m in res.iter() {
            if c.is_overlapped(&m) {
                flag = false;
                break;
            }
        }

        if flag {
            res.push(c);
        }
    });
    // res
    // 进保留前二分之一
    res[0..res.len() / 2]
        .into_iter()
        .map(|c| c.clone())
        .collect::<Vec<Complex<T>>>()
}

// !!!! not stable
pub fn confucion_matrix<T>(refc: Vec<Complex<T>>, prec: Vec<Complex<T>>)
where
    T: Eq + Hash + Clone + Ord,
{
    let mut matrix_os = vec![vec![f64::default(); prec.len()]; refc.len()];
    let mut matrix_num = vec![vec![usize::default(); prec.len()]; refc.len()];
    for i in 0..refc.len() {
        for j in 0..prec.len() {
            let (commisze, os) = refc[i].overlap_score(&prec[i]);
            matrix_os[i][j] = os;
            matrix_num[i][j] = commisze as usize;
        }
    }

    // precision
    let mut count = 0;
    for i in 0..prec.len() {
        for j in 0..refc.len() {
            if matrix_os[j][i].ge(&THRESHOLD_OS) {
                count += 1;
                break;
            }
        }
    }
    let precison = count as f64 / prec.len() as f64;
    println!("precision: {:.4?}", precison);

    let mut count = 0;
    for i in 0..refc.len() {
        for j in 0..prec.len() {
            if matrix_os[j][i].ge(&THRESHOLD_OS) {
                count += 1;
                break;
            }
        }
    }
    let reccall = count as f64 / refc.len() as f64;
    println!("Recall: {:.4?}", reccall);
    let f1 = 2.0 * precison * reccall / (reccall + precison);
    println!("fscore: {:.4?}", f1);

    // 分子
    // let sum1 = matrix_num.iter().map(|v| *v.iter().max().unwrap()).sum::<usize>();
    let sum1 = matrix_num
        .iter()
        .map(|v| *v.iter().max().unwrap())
        .sum::<usize>();
    let sum2 = refc.iter().map(|c| c.proteins.len()).sum::<usize>();
    let sn = sum1 as f64 / sum2 as f64;
    println!("Sn: {:.4?}", sn);

    let mut sum1 = 0;
    for i in 0..prec.len() {
        let mut max = 0;
        for j in 0..refc.len() {
            max = usize::max(max, matrix_num[j][i]);
        }
        sum1 += max;
    }
    let sum2 = prec.iter().map(|c| c.proteins.len()).sum::<usize>();
    let ppv = sum1 as f64 / sum2 as f64;
    println!("PPV: {:.4}", ppv);
    let acc = ppv.sqrt();
    println!("Acc: {:.4?}", acc);
}

// 泛型读取 complex 的函数
fn read_complex<T>(complex_file: &str, proteins: &HashSet<T>) -> Vec<Complex<T>>
where
    T: Clone + Eq + Hash + From<String> + Ord,
{
    let contents = read_to_string(complex_file).expect("Failed to read complex file!");
    let mut complexes = Vec::new();

    for line in contents.lines() {
        let filtered: Vec<T> = line
            .split_whitespace()
            .map(|s| s.to_string())
            .map(T::from)
            .filter(|item| proteins.contains(item))
            .collect();

        if filtered.len() >= MIN_SIZE && filtered.len() <= MAX_SIZE {
            complexes.push(Complex::new(filtered, 0.0));
        }
    }

    complexes
}

#[cfg(test)]
mod tests {
    use std::{collections::HashSet, fmt::Debug};

    use crate::{eva::COMPLEX_REF, graph::Graph};

    use super::{confucion_matrix, read_complex, update_by_cohesion, Complex};

    #[test]
    fn test_update_by_cohesion() {
        let complexes = vec![
            Complex::new(vec![1, 2, 3, 4, 5], 0.8),
            Complex::new(vec![1, 2, 3, 4], 0.9),
        ];
        let res = update_by_cohesion(complexes);
        println!("{:?}", res);
    }

    #[test]
    fn test_complex() {
        let graph = Graph::new_from_file("./data/collins/collins.txt", false);
        let proteins = graph.id_protein.clone().into_iter().collect::<HashSet<_>>();
        let complex = read_complex(COMPLEX_REF, &proteins);
        // println!("{:?}", complex);
        println!("{:?}", complex.len());
    }

    #[test]
    fn test_complex_os() {
        let c1 = Complex::new(vec![1, 2, 3, 4], Default::default());
        let c2 = Complex::new(vec![2, 3, 4], Default::default());
        println!("{:?}", c1.overlap_score(&c2));
    }
    #[test]
    fn test_confusion_matrix() {
        let refc = vec![
            Complex::new(vec!["A".to_string(), "B".to_string(), "C".to_string()], 0.8),
            Complex::new(vec!["B".to_string(), "C".to_string(), "D".to_string()], 0.9),
        ];

        let prec = vec![
            Complex::new(vec!["A".to_string(), "B".to_string(), "C".to_string()], 0.7),
            Complex::new(vec!["C".to_string(), "D".to_string(), "E".to_string()], 0.8),
        ];
        println!("{}", prec[0]);

        confucion_matrix(refc, prec);
    }

    #[test]
    fn test_confusion_matrix2() {
        let graph = Graph::new_from_file("./data/collins/collins.txt", false);
        let proteins = graph.id_protein.clone().into_iter().collect::<HashSet<_>>();
        let refc = read_complex("./data/collins/complex.txt", &proteins);
        let prec = read_complex("rs_pcegs.txt", &proteins);

        println!("refc: {}, prec: {}", refc.len(), prec.len());
        confucion_matrix(refc, prec);
    }

    #[test]
    fn test_cohesion() {
        let g = Graph::new_from(vec![(0, 1, 0.6), (1, 2, 0.8), (3, 4, 0.1)]);

        let cluster = vec![0, 1, 2, 4];
        let cohesion = g.calculate_cohesion(&cluster);
        println!("cohesion: {}", cohesion)
    }
}
