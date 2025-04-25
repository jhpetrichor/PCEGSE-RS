use std::{fs::File, io::Write};

use essential_protein::{
    dag::weight_by_dag_topo, eva::Complex, gene_expression::get_dpins, graph::Graph, pcegs,
};

fn main() {
    let graph = Graph::new_from_file("./data/krogan_core/krogan_core.txt", false);
    // 动态网络
    let dpins = get_dpins(&graph);

    let mut complexes = Vec::new();
    for mut dp in dpins.into_iter() {
        // weight_by_dag(&mut dp);
        weight_by_dag_topo(&mut dp, 0.5);
        let res = pcegs::pcegs(&dp, 0.4);
        // res.into_iter().for_each(|c| println!("{}", c));
        complexes.extend(res);
    }
    let file = format!("result/krogan_core.txt");
    write_file(file, complexes);

    // let mut alpha = 0.;

    // while (1.04).ge(&alpha) {
    //     let mut beta = 0.;
    //     while (1.04).ge(&beta) {
    //         let graph = Graph::new_from_file("./data/krogan_core/krogan_core.txt", true);
    //         // 动态网络
    //         let dpins = get_dpins(&graph);

    //         let mut complexes = Vec::new();
    //         for mut dp in dpins.into_iter() {
    //             // weight_by_dag(&mut dp);
    //             weight_by_dag_topo(&mut dp, alpha);
    //             let res = pcegs::pcegs(&dp, beta);
    //             // res.into_iter().for_each(|c| println!("{}", c));
    //             complexes.extend(res);
    //         }
    //         let file = format!("result/krogan_core_{:.2}_{:.2}.txt", alpha, beta);
    //         write_file(file, complexes);
    //         beta += 0.05;
    //     }
    //     alpha += 0.05;
    // }
}

fn write_file(file: String, complexes: Vec<Complex<String>>) {
    let mut file = File::create(file).expect("Failed to crate faile!");
    // file.write_all(buf)
    for complex in complexes {
        if complex.cohesion.le(&0.2) {
            continue;
        }
        let mut str = format!(
            "{} {:.4}",
            complex
                .proteins
                .iter()
                .map(|p| format!("{}", p))
                .collect::<Vec<String>>()
                .join("\t"),
            complex.cohesion
        );
        str.push('\n');
        file.write_all(str.as_bytes())
            .expect("Failed to write to file!");
    }
}

#[test]
fn teste_eva() {
    let mut commands = String::new();

    let mut a = 0.;
    while (1.04).ge(&a) {
        let mut b = 0.;
        while (1.04).ge(&b) {
            commands = format!(
                "{} python evaluation\\match.py -n gavin result\\gavin_{:.2?}_{:.2}.txt \n",
                commands, a, b
            );
            b += 0.05;
        }
        a += 0.5;
    }
    let mut file = File::create("gavin.ps1").expect("msg");
    file.write_all(commands.as_bytes()).unwrap();
}
