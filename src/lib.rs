use std::error::Error;
pub mod config;
pub mod locator;
pub mod reference;

pub type BoxError = Box<dyn Error + Send + Sync>;
