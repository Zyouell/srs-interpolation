//! This module contains the error types used in this library.
use ark_std::{
    error::Error,
    fmt::{Display, Formatter, Result},
};
/// Enum for errors that can occur during interpolation.
#[derive(Debug, Clone)]
pub enum InterpolationError {
    /// Invalid parameters were provided to the interpolation function.
    InvalidParameters(String),
    /// Error performing finite field operation.
    FieldError(String),
    /// The size of the provided SRS was not a power of two.
    SizeError,
}

impl Display for InterpolationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            InterpolationError::InvalidParameters(s) => {
                write!(f, "Invalid parameters: {}", s)
            }
            InterpolationError::FieldError(s) => {
                write!(f, "Field error: {}", s)
            }
            InterpolationError::SizeError => {
                write!(
                    f,
                    "Size error: the provided SRS size was not a power of two"
                )
            }
        }
    }
}

impl Error for InterpolationError {}
