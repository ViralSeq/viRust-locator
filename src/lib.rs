//! This crate provides a set of utilities for generating sequence locators for HIV/SIV sequences
//! resembling the LANL HIV-locator tool.

use std::error::Error;
pub mod config;
pub mod locator;
pub mod reference;

///! This crate provides a set of utilities for working with the [OpenTelemetry](https://opentelemetry.io/) ecosystem.
pub type BoxError = Box<dyn Error + Send + Sync>;
