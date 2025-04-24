use core::panic;
use std::{
    ops::{Add, AddAssign, Index},
    usize,
};

#[derive(Debug, PartialEq, Clone)]
pub struct Matrix<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
}

impl<T> Matrix<T>
where
    T: Add + AddAssign + Eq + Copy,
{
    pub fn new(data: impl Into<Vec<T>>, rows: usize, cols: usize) -> Self {
        Self {
            data: data.into(),
            rows,
            cols,
        }
    }
}

impl<T> Index<[usize; 2]> for Matrix<T> {
    type Output = T;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        let index = index[0] * self.cols + index[1];
        if index >= self.data.len() {
            panic!("Index out of bounds: {}, {}", index, self.data.len());
        }
        return &self.data[index];
    }
}

impl<T> Add for Matrix<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.cols != rhs.cols || self.rows != rhs.rows {
            panic!("Shape error");
        }
        let data = self
            .data
            .iter()
            .zip(rhs.data.iter())
            .map(|(a, b)| *a + *b)
            .collect::<Vec<_>>();

        Matrix {
            data,
            rows: self.rows,
            cols: self.cols,
        }
    }
}

impl<T> AddAssign for Matrix<T> {
    fn add_assign(&mut self, rhs: Self) {}
}

#[cfg(test)]
mod tests {
    use super::Matrix;

    #[test]
    fn matrix_new() {
        let mut a = Matrix::new([1, 2, 3, 4, 5, 6], 2, 3);
        assert_eq!(a[[0, 1]], 2);
        assert_eq!(a[[1, 1]], 5);

        let b = Matrix::new([1, 2, 3, 4, 5, 6], 2, 3);
        assert_eq!(
            Matrix::new([2, 4, 6, 8, 10, 12], 2, 3),
            a.clone() + b.clone()
        );

        a += b;
        assert_eq!(Matrix::new([2, 4, 6, 8, 10, 12], 2, 3), a);
    }

    #[test]
    #[should_panic]
    fn matrix_new1() {
        let a = Matrix::new([1, 2, 3, 4, 5, 6], 2, 3);
        assert_eq!(a[[0, 1]], 2);
        assert_eq!(a[[4, 4]], 5);
    }
}
