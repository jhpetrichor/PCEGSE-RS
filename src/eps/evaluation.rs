use std::{collections::HashSet, hash::Hash};

pub fn evalustion<T>(refprotein: &HashSet<T>, preprotein: &[T])
where
    T: Clone + Hash + Eq,
{
    let (mut tn, mut tp, mut ff, mut fp) = (0, 0, 0, 0);

    preprotein[0..refprotein.len() + 1]
        .into_iter()
        .for_each(|c| {
            if refprotein.contains(c) {
                tp += 1;
            } else {
                fp += 1;
            }
        });

    preprotein[refprotein.len() + 1..]
        .into_iter()
        .for_each(|c| {
            if refprotein.contains(c) {
                ff += 1;
            } else {
                tn += 1;
            }
        });

    // println!("my {}, {}, {}, {}", tn, tp, ff, fp);

    let sn = tp as f64 / (tp + ff) as f64;
    let sp = tn as f64 / (tn + fp) as f64;
    let ppv = tp as f64 / (tp + fp) as f64;
    let npv = tn as f64 / (ff + tn) as f64;
    let f = 2. * sn * ppv / (sn + ppv);
    let acc = (tp + tn) as f64 / (tp + tn + fp + ff) as f64;

    println!(
        "{:.4} {:.4} {:.4} {:.4} {:.4} {:.4}",
        sn, sp, ppv, npv, f, acc
    );
}

pub fn top_n<T>(refprotein: &HashSet<T>, preprotein: &[T])
where
    T: Eq + Hash + Clone,
{
    let top_n = [
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
    ];
    let pre_cent = [1, 5, 10, 15, 20, 25];

    let mut top = Vec::new();
    let mut top_add = Vec::new();
    top_n.into_iter().for_each(|k| {
        let pre = preprotein[0..k]
            .into_iter()
            .map(|c| c.to_owned())
            .collect::<HashSet<_>>()
            .intersection(&refprotein)
            .count();
        top.push(pre);
        if top_add.is_empty() {
            top_add.push(pre);
        } else {
            top_add.push(top.last().unwrap() - top[top.len() - 2]);
        }
    });
    top_add.push(top.last().unwrap().clone());
    println!("top: {:?}", top);
    println!("top_add: {:?}", top_add);

    // println!("top%{}")
    top.clear();
    top_add.clear();
    pre_cent.into_iter().for_each(|p| {
        let pre = preprotein[0..preprotein.len() * p / 100]
            .into_iter()
            .map(|c| c.to_owned())
            .collect::<HashSet<_>>()
            .intersection(&refprotein)
            .count();
        top.push(pre);
        if top_add.is_empty() {
            top_add.push(pre);
        } else {
            top_add.push(top.last().unwrap() - top[top.len() - 2]);
        }
    });

    println!("top%: {:?}", top);
    println!("top%_add: {:?}", top_add);
}

// 生成n个点
pub fn roc<T>(refprotein: &HashSet<T>, preprotein: &[T], n: usize) -> [Vec<f64>; 2]
where
    T: Eq + Hash + Clone,
{
    let mut roc = [Vec::new(), Vec::new()];

    let k = (0..=100).step_by(100 / n).collect::<Vec<_>>();

    for i in k.iter() {
        let (mut tn, mut tp, mut ff, mut fp) = (0, 0, 0, 0);
        preprotein[0..preprotein.len() * i / 100]
            .into_iter()
            .for_each(|c| {
                if refprotein.contains(c) {
                    tp += 1;
                } else {
                    fp += 1;
                }
            });

        preprotein[preprotein.len() * i / 100..]
            .into_iter()
            .for_each(|c| {
                if refprotein.contains(c) {
                    ff += 1;
                } else {
                    tn += 1;
                }
            });
        roc[0].push(fp as f64 / (fp + tn) as f64);
        roc[1].push(tp as f64 / (tp + ff) as f64);
    }
    roc
}

// precision_recall曲线
pub fn pr<T>(refprotein: &HashSet<T>, preprotein: &[T], n: usize) -> [Vec<f64>; 2]
where
    T: Eq + Hash + Clone,
{
    let mut pr = [Vec::new(), Vec::new()];

    let k = (0..=100).step_by(100 / n).collect::<Vec<_>>();

    for i in k.iter() {
        let (mut tn, mut tp, mut ff, mut fp) = (0, 0, 0, 0);
        preprotein[0..preprotein.len() * i / 100]
            .into_iter()
            .for_each(|c| {
                if refprotein.contains(c) {
                    tp += 1;
                } else {
                    fp += 1;
                }
            });

        preprotein[preprotein.len() * i / 100..]
            .into_iter()
            .for_each(|c| {
                if refprotein.contains(c) {
                    ff += 1;
                } else {
                    tn += 1;
                }
            });

        if tp + fp == 0 {
            pr[0].push(0.);
        } else {
            pr[0].push(tp as f64 / (tp + fp) as f64);
        }
        if tp + ff == 0 {
            pr[1].push(0.);
        } else {
            pr[1].push(tp as f64 / (tp + ff) as f64);
        }
    }
    pr
}
