use log::debug;
use ndarray::{Array1, Array2};

use crate::eps::graph::Graph;

pub trait RandomWalk {
    const MAXITER: usize = usize::MAX;
    fn get_rw_matrix(&self, alpha: f64, beta: f64) -> Array2<f64>;
    fn random_walk(&self, init_state: Array1<f64>, alpha: f64, beta: f64) -> Array1<f64>;
}

impl RandomWalk for Graph {
    fn get_rw_matrix(&self, alpha: f64, beta: f64) -> Array2<f64> {
        debug!("{}", self.ppg);
        let c = self.ppg.clone() * 0.85
            + 0.15 / self.ppg.dim().0 as f64 * Array2::<f64>::ones(self.ppg.dim());
        let w = self.pcg.dot(&self.cpg);
        let b = self.p2g.dot(&self.g2p); // 蛋白质到go terms之间的游走矩阵

        c * alpha + w * beta + b * (1. - alpha - beta)
    }

    fn random_walk(&self, init_state: Array1<f64>, alpha: f64, beta: f64) -> Array1<f64> {
        println!("a: {}, b: {}", alpha, beta);
        let matrix = self.get_rw_matrix(alpha, beta);
        let mut cu_state = init_state;
        for _ in 0..Self::MAXITER {
            let temp_state = cu_state.dot(&matrix);
            if (&temp_state - &cu_state).mapv(|x| x.abs()).sum() < 1e-10 {
                break;
            }
            cu_state = temp_state;
            // if i % 100 == 0 {
            //     info!("round {i}: {}", cu_state);
            // }
        }

        cu_state
    }
}
