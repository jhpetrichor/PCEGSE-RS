use std::env;

use essential_protein::eps::{
    errors::Result,
    evaluation::{evalustion, top_n},
    graph::Graph,
    options::Options,
    rw::RandomWalk,
};

use log::{info, LevelFilter};
use ndarray::Array1;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    let ppi = args[1].as_str();

    // 设置目录级别
    env_logger::builder().filter_level(LevelFilter::Info).init();

    let mut options = Options::new(ppi);
    options.min_size = 3;
    info!("{:?}", options);
    let g = Graph::new_with(options)?;

    let mut a = Array1::<f64>::zeros(g.p_count);
    a[0] = 0.5;
    a[5] = 0.4;
    a[9] = 0.1;

    let iter = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

    for i in &iter {
        for j in &iter {
            if i + j > 1.0 {
                continue;
            }
            let scores = g.random_walk(a.clone(), *i, *j);
            let essen_score = g.essentiality_rank(&scores);

            evalustion(&g.essential_proteins, &essen_score);
            top_n(&g.essential_proteins, &essen_score);
        }
    }

    Ok(())
}
