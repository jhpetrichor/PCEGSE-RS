const DATA_PREFIX: &str = "./data/";

#[derive(Debug)]
pub struct Options {
    pub ppi_file: String,
    pub clique_file: String,
    pub label_file: String,
    pub min_size: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            ppi_file: String::from("./data/collins/collins.txt"),
            clique_file: String::from("./data/collins/collins_clique.txt"),
            label_file: String::from("./data/collins/collins_labels"),
            min_size: 3,
        }
    }
}

impl Options {
    pub fn new(data_name: &str) -> Self {
        let mut options = Self::default();
        options.ppi_file = format!("{}{}/{}.txt", DATA_PREFIX, data_name, data_name);
        options.clique_file = format!("{}{}/{}_clique.txt", DATA_PREFIX, data_name, data_name);
        options.label_file = format!("{}{}/{}_label.txt", DATA_PREFIX, data_name, data_name);
        options
    }
}

pub fn get_essential_protein_count() -> usize {
    todo!()
}

#[cfg(test)]
mod test {
    use super::Options;

    #[test]
    fn options_new_test() {
        let options = Options::new("collins");
        assert_eq!(options.ppi_file, String::from("./data/collins/collins.txt"));
        assert_eq!(
            options.clique_file,
            String::from("./data/collins/collins_clique.txt")
        );
        assert_eq!(
            options.label_file,
            String::from("./data/collins/collins_label.txt")
        );
    }
}
