use log::error;

use crate::eps::errors::{Errors, Result};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs,
    hash::Hash,
    vec,
};

#[derive(Debug)]
pub struct BiMap<E> {
    entity_id: HashMap<E, usize>,
    id_entity: HashMap<usize, E>,
}

impl<E> BiMap<E>
where
    E: Eq + Hash + Clone + Ord + Default,
{
    pub fn new() -> Self {
        Self {
            entity_id: HashMap::new(),
            id_entity: HashMap::new(),
        }
    }

    pub fn len(&self) -> usize {
        assert_eq!(self.id_entity.len(), self.entity_id.len());
        self.id_entity.len()
    }

    pub fn contains_entity(&self, entity: &E) -> bool {
        self.entity_id.contains_key(entity)
    }

    pub fn get_id(&self, entity: &E) -> Option<&usize> {
        self.entity_id.get(entity)
    }

    pub fn get_entity(&self, id: &usize) -> Option<&E> {
        self.id_entity.get(id)
    }

    /// 插入，如果已经存在，则返回id值
    pub fn insert(&mut self, entity: E) -> usize {
        match self.entity_id.get(&entity) {
            None => {
                let id = self.entity_id.len();
                self.entity_id.insert(entity.clone(), id);
                self.id_entity.insert(id, entity);

                return id;
            }
            Some(id) => return id.to_owned(),
        }
    }

    pub fn remove_by_entity(&mut self, entity: &E) {
        if let Some(id) = self.entity_id.remove(entity) {
            self.id_entity.remove(&id);
        }
    }

    pub fn remove_by_id(&mut self, id: &usize) {
        if let Some(entity) = self.id_entity.remove(id) {
            self.entity_id.remove(&entity);
        }
    }

    pub fn ids(&self) -> BTreeSet<&usize> {
        self.id_entity.keys().collect()
    }

    pub fn entitys(&self) -> BTreeSet<&E> {
        self.entity_id.keys().collect()
    }

    pub fn into_vec(&self) -> Vec<E> {
        let mut v = vec![E::default(); self.entity_id.len()];
        for (id, entity) in self.id_entity.iter() {
            v[*id] = entity.clone();
        }
        v
    }
}

pub fn read_ppi_file(file: &str) -> Result<(BiMap<String>, BTreeMap<usize, BTreeSet<usize>>)> {
    let reader = match fs::read_to_string(file) {
        Ok(reader) => reader,
        Err(e) => {
            return {
                error!("{e}");
                Err(crate::eps::errors::Errors::FailedToReadPPIFile)
            }
        }
    };

    let mut bimap = BiMap::<String>::new();
    let mut neighbor = BTreeMap::<usize, BTreeSet<usize>>::new();
    reader.lines().into_iter().for_each(|e| {
        let pair = e.split_whitespace().collect::<Vec<_>>();
        let a = bimap.insert(pair[0].to_string());
        let b = bimap.insert(pair[1].to_string());
        // edge.push((a, b));
        neighbor
            .entry(a)
            .and_modify(|c| {
                c.insert(b);
            })
            .or_insert_with(|| BTreeSet::from([b]));
        neighbor
            .entry(b)
            .and_modify(|c| {
                c.insert(a);
            })
            .or_insert_with(|| BTreeSet::from([a]));
    });

    Ok((bimap, neighbor))
}

// 读取所有go term的索引
pub fn get_all_go_terms() -> BiMap<String> {
    let reader1 =
        fs::read_to_string("./data/is_a_child_ancestor.txt").expect("Failed to open go terms");
    let reader2 =
        fs::read_to_string("./data/part_of_child_ancestor.txt").expect("Failed to open go terms");
    let reader = format!("{}\n{}", reader1, reader2);
    let mut bimap = BiMap::new();
    reader.lines().into_iter().for_each(|line| {
        let line = line.split_whitespace().collect::<Vec<_>>();
        line.into_iter().for_each(|c| {
            bimap.insert(c.to_string());
        });
    });

    bimap
}

// 读取蛋白质及其相关联的go term
pub fn read_protein_go(
    pid: &BiMap<String>,
    gid: &BiMap<String>,
) -> Result<(
    BTreeMap<usize, BTreeSet<usize>>,
    BTreeMap<usize, BTreeSet<usize>>,
)> {
    let reader = match fs::read_to_string("./data/protein-go.txt") {
        Ok(reader) => reader,
        Err(e) => {
            error!("{}", e);
            return Err(Errors::FailedToReadPPIFile);
        }
    };
    let mut protein_go = BTreeMap::new();
    let mut go_protein = BTreeMap::<usize, BTreeSet<usize>>::new();
    for line in reader.lines() {
        let line = line.split_whitespace().collect::<Vec<_>>();
        let protein = line[0].to_string();
        if !pid.contains_entity(&protein) {
            continue;
        }
        // 蛋白质
        let protein = *pid.get_id(&protein).unwrap();
        let mut gos = BTreeSet::new();
        for g in line[1..line.len()].into_iter() {
            if !gid.contains_entity(&g.to_string()) {
                continue;
            }
            // go term
            let term = *gid.get_id(&g.to_string()).unwrap();
            gos.insert(term);
            go_protein
                .entry(term)
                .and_modify(|s| {
                    s.insert(protein);
                })
                .or_insert_with(|| BTreeSet::from([protein]));
        }
        protein_go.insert(protein, gos);
    }

    Ok((protein_go, go_protein))
}

// 读取子项到父项的映射, 这一步需要构建go term到id的映射
// 返回go terms的映射和go terms之间的连接关系
pub fn read_go_edge(file: &str, gid: &BiMap<String>) -> Result<BTreeMap<usize, BTreeSet<usize>>> {
    let reader = match fs::read_to_string(file) {
        Ok(reader) => reader,
        Err(e) => {
            error!("{}", e);
            return Err(Errors::FailedToReadPPIFile);
        }
    };
    // todo 后继可以考虑同时在父项和子项节点中游走
    let mut gancestor = BTreeMap::<usize, BTreeSet<usize>>::new();

    for line in reader.lines() {
        // println!("{line}");
        let gos = line.split_whitespace().collect::<Vec<_>>();
        if !gid.contains_entity(&gos[0].to_string()) {
            println!("not contains, {}", gos[0]);
            continue;
        }
        let child = gid.get_id(&gos[0].to_string()).unwrap();
        for i in 1..gos.len() {
            if !gid.contains_entity(&gos[i].to_string()) {
                continue;
            }
            let ancestor = gid.get_id(&gos[i].to_string()).unwrap();
            // gedge.push((*child, *ancestor));
            gancestor
                .entry(*child)
                .and_modify(|c| {
                    c.insert(*ancestor);
                })
                .or_insert_with(|| BTreeSet::from([*ancestor]));
            // gancestor
            //     .entry(*ancestor)
            //     .and_modify(|c| {
            //         c.insert(*child);
            //     })
            //     .or_insert_with(|| BTreeSet::from([*child]));
        }
    }

    Ok(gancestor)
    // Ok((bimap, edge))
}

// 高阶节点到id的映射
// 添加蛋白质到高阶节点的映射
pub fn read_clique(
    file: &str,
    pid: &BiMap<String>,
    min_size: usize,
) -> Result<(BiMap<BTreeSet<usize>>, BTreeMap<usize, BTreeSet<usize>>)> {
    let reader = match fs::read_to_string(file) {
        Ok(reader) => reader,
        Err(e) => {
            error!("{}", e);
            return Err(Errors::FailedToReadPPIFile);
        }
    };

    let mut p_c = BTreeMap::<usize, BTreeSet<usize>>::new();
    let mut hnode = BiMap::new();
    let proteins = pid.entitys();
    // reader.lines().into_iter().for_each(|line| {
    for line in reader.lines() {
        let c = line
            .split_whitespace()
            .filter(|a| proteins.contains(&a.to_string()))
            .collect::<Vec<_>>();
        let clique: BTreeSet<usize> = c.iter().map(|p| pid.entity_id[&p.to_string()]).collect();
        if clique.len() <= min_size {
            continue;
        }
        let c_id = hnode.insert(clique.clone());
        clique.into_iter().for_each(|p| {
            p_c.entry(p)
                .and_modify(|m| {
                    m.insert(c_id);
                })
                .or_insert_with(|| BTreeSet::from([c_id]));
        });
    }
    //蛋白质到高阶节点的映射， 高阶邻居

    Ok((hnode, p_c))
}

/// 从文件中读取标签，返回蛋白质对应的id --> label
/// label(bool) 表示蛋白质是否为关键蛋白质
pub fn read_label(file: &str, pid: &BiMap<String>) -> Result<Vec<bool>> {
    let reader = match fs::read_to_string(file) {
        Ok(reader) => reader,
        Err(e) => {
            error!("{}", e);
            return Err(Errors::FailedToReadLableFile);
        }
    };
    let mut result_label = vec![false; pid.len()];
    for line in reader.lines() {
        let line = line.split_whitespace().collect::<Vec<_>>();
        let (protein, label) = (line[0].parse::<String>(), line[1].parse::<usize>());
        match (protein, label) {
            (Ok(protein), Ok(label)) => {
                let id = pid.get_id(&protein).unwrap();
                let l = if label == 1 { true } else { false };
                result_label[*id] = l;
            }
            _ => return Err(Errors::FailedToParseLabel),
        }
    }

    Ok(result_label)
}

#[cfg(test)]
mod tests {
    use crate::eps::utils::{read_clique, read_go_edge, read_label, read_protein_go};

    use super::{get_all_go_terms, read_ppi_file, BiMap};

    #[test]
    fn bimap_into_vec_test() {
        let mut bimap = BiMap::<String>::new();
        bimap.insert(String::from("a"));
        bimap.insert(String::from("b"));
        bimap.insert(String::from("c"));
        bimap.insert(String::from("d"));

        let v = bimap.into_vec();
        assert_eq!(v.len(), 4);
        assert_eq!(
            v,
            vec![
                String::from("a"),
                String::from("b"),
                String::from("c"),
                String::from("d")
            ]
        );
        println!("{:?}", v);
    }

    #[test]
    fn bimap_test() {
        let mut bmap = BiMap::<String>::new();
        let id = bmap.insert("Hello".to_string());
        assert_eq!(id, 0);
        let id = bmap.get_id(&"Hello".to_string());
        assert_eq!(id, Some(&0));
        let entity = bmap.get_entity(&0);
        assert_eq!(entity, Some(&"Hello".to_string()));

        let id = bmap.insert("World".to_string());
        assert_eq!(id, 1);
        let id = bmap.get_id(&"World".to_string());
        assert_eq!(id, Some(&1));
        let entity = bmap.get_entity(&1);
        assert_eq!(entity, Some(&"World".to_string()));

        bmap.remove_by_entity(&"World".to_string());
        assert_eq!(None, bmap.get_entity(&1));
        assert_eq!(None, bmap.get_id(&"World".to_string()));

        bmap.remove_by_id(&0);
        assert_eq!(None, bmap.get_entity(&0));
        assert_eq!(None, bmap.get_id(&"Hello".to_string()));
    }

    #[test]
    fn read_ppi_test() {
        let (pid, pedge) =
            read_ppi_file("./data/collins/collins.txt").expect("Failed to read ppi file");
        assert_eq!(pid.len(), 1622);
        assert_eq!(pedge.len(), 1622);
    }

    #[test]
    fn read_label_test() {
        let (pid, _) =
            read_ppi_file("./data/collins/collins.txt").expect("Failed to read ppi_file");
        let label = read_label("./data/collins/collins_label.txt", &pid);
        assert_eq!(label.unwrap().len(), pid.len());
        // println!("{:?}", label);
    }

    #[test]
    fn read_clique_test() {
        let (pid, _) =
            read_ppi_file("./data/collins/collins.txt").expect("Failed to read ppi file!");

        let clique_res = read_clique("./data/collins/collins_clique.txt", &pid, 3);
        assert!(clique_res.is_ok());
        let clique = clique_res.unwrap().0;
        assert_eq!(clique.len(), 352);
    }

    #[test]
    fn read_go_file_test() {
        let go_id = get_all_go_terms();
        assert_eq!(26653, go_id.len());
        let g_edge = read_go_edge("./data/is_a_child_ancestor.txt", &go_id)
            .expect("Failed to open DAG go file!");
        assert_eq!(g_edge.len(), 26478);
    }

    #[test]
    fn get_all_go_terms_test() {
        let go_terms = get_all_go_terms();
        assert_eq!(26653, go_terms.len());
    }

    #[test]
    fn read_protein_go_test() {
        // 读取蛋白质pid
        let (pid, _p_edge) =
            read_ppi_file("./data/collins/collins.txt").expect("Failed to open ppi file!");
        // 读取go term id
        let gid = get_all_go_terms();
        // 读取蛋白质对应的go_id
        let p_go = read_protein_go(&pid, &gid).expect("Failed to open ");
        assert_eq!(p_go.0.len(), 1622);
        assert_eq!(p_go.1.len(), 2781);
        // println!("{:?}", p_go);
    }
}
