use std::fs::read_to_string;

#[allow(unused)]
pub fn read_write_compelx(data: &str) {
    //
    let (collins_complex, gavin_complex, krogan_complex) = (
        read_to_string(data).unwrap(),
        read_to_string(data).unwrap(),
        read_to_string(data).unwrap(),
    );
}
