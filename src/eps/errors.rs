use std::result;

use thiserror::Error;

pub type Result<T> = result::Result<T, Errors>;

#[derive(Error, PartialEq, Eq, Debug)]
pub enum Errors {
    #[error("Failed to open PPI file")]
    FailedToReadPPIFile,

    #[error("Failed to parse string")]
    FailedToParseLabel,

    #[error("Failed to read label file")]
    FailedToReadLableFile,
}
