use std::{collections::HashSet, env, fs::read_to_string};

use essential_protein::eps::evaluation::{evalustion, pr, roc, top_n};

fn read_essential_protein() -> HashSet<String> {
    let mut essentil_proteins = HashSet::<String>::new();

    let contents = read_to_string("./data/essential proteins.ref").expect("读取关键蛋白质失败");
    for line in contents.lines() {
        let ps = line.split_whitespace().collect::<Vec<_>>();
        essentil_proteins.insert(ps.first().unwrap().to_string());
    }

    essentil_proteins
}

// read prediction result
fn read_pre_file(pre_file: &str) -> Vec<String> {
    let mut res = Vec::<String>::new();

    for line in read_to_string(pre_file).expect("读取预测结果失败").lines() {
        let ps = line.split_whitespace().collect::<Vec<_>>();
        res.push(ps.first().unwrap().to_string());
    }
    res
}

fn main() {
    let arags = env::args().collect::<Vec<String>>();

    let refprotein = read_essential_protein();
    let preprotein = read_pre_file(arags[1].as_str());

    // 真是有的关键蛋白质
    let refprotein = refprotein
        .intersection(&preprotein.clone().into_iter().collect())
        .map(|s| s.clone())
        .collect::<HashSet<_>>();

    evalustion(&refprotein, &preprotein);
    top_n(&refprotein, &preprotein);
    let roc = roc(&refprotein, &preprotein, 50);
    println!("{:?} \n {:?}", roc[0], roc[1]);
    let pr = pr(&refprotein, &preprotein, 50);
    println!("{:?} \n {:?}", pr[0], pr[1]);
}
