//! This module defines the `Args` struct and related functionality for parsing and validating
//! command-line arguments for the `viRust-locator` CLI tool. The tool is a simple implementation
//! of LANL's HIV locator in Rust.
//!
//! # Structs
//!
//! - `Args`: Represents the command-line arguments for the tool. It includes fields for the query
//!   sequence, reference genome, query type, and algorithm choice.
//!
//! # Functions
//!
//! - `get_styles`: Configures and returns custom styles for the CLI output, including styles for
//!   usage, headers, literals, errors, and placeholders.
//!
//! - `Args::validate`: Validates the parsed arguments to ensure they meet the expected criteria,
//!   such as valid query types (`nt` or `aa`), valid reference genomes (`HXB2` or `SIVmm239`),
//!   and valid nucleotide or amino acid sequences.
//!
//! # Command-Line Arguments
//!
//! - `--query` (`-q`): Specifies the query sequence. This can be a nucleotide or amino acid
//!   sequence, depending on the `--type-query` argument.
//!
//! - `--reference` (`-r`): Specifies the reference genome. The default value is `HXB2`. Valid
//!   options are `HXB2` or `SIVmm239`.
//!
//! - `--type-query` (`-t`): Specifies the type of the query sequence. The default value is `nt`
//!   (nucleotide). Valid options are `nt` or `aa` (amino acid).
//!
//! - `--algorithm` (`-a`): Specifies the algorithm to use for the locator. The default value is `1`.
//!   Valid options are `1` (accurate but slower) or `2` (fast but less accurate, suitable for smaller
//!   query sequences).
//!
//! # Validation Rules
//!
//! - The `type_query` must be either `nt` or `aa`.
//! - The `algorithm` must be either `1` or `2`.
//! - The `reference` must be either `HXB2` or `SIVmm239`.
//! - For nucleotide sequences (`nt`):
//!   - The sequence must conform to the IUPAC nucleotide alphabet.
//!   - The sequence length must be greater than 3.
//! - For amino acid sequences (`aa`):
//!   - The sequence must conform to the IUPAC protein alphabet.
//!   - The sequence length must be greater than 3.
//!
//! # Errors
//!
//! The `Args::validate` function returns an error message if any of the validation rules are
//! violated, such as invalid query types, invalid sequences, or unsupported reference genomes.
use bio::alphabets;
use clap::builder::styling::{AnsiColor, Color};
use clap::builder::styling::{Style, Styles};
use clap::{ColorChoice, Parser};

#[derive(Parser, Debug)]
#[command(
    name = "viRust-locator",
    version = "0.1.0",
    about = "\x1b[1;91mSimple LANL's HIV locator tool implementation in Rust CLI\x1b[0m",
    color = ColorChoice::Always,
    styles = get_styles(),
)]
pub struct Args {
    /// Query sequence
    #[arg(short, long)]
    pub query: String,

    /// Reference genome, either HXB2 or SIVmm239
    #[arg(short, long, default_value = "HXB2")]
    pub reference: String,

    /// Type of query, either nt or aa
    #[arg(short, long, default_value = "nt")]
    pub type_query: String,

    /// algorithm for locator, 1 is accurate but slower, 2 is fast but less accurate, suitable for smaller query sequences
    #[arg(short, long, default_value_t = 1)]
    pub algorithm: u8,
}

pub fn get_styles() -> Styles {
    Styles::styled()
        .usage(
            Style::new()
                .bold()
                .underline()
                .fg_color(Some(Color::Ansi(AnsiColor::Yellow))),
        )
        .header(
            Style::new()
                .bold()
                .underline()
                .fg_color(Some(Color::Ansi(AnsiColor::Yellow))),
        )
        .literal(Style::new().fg_color(Some(Color::Ansi(AnsiColor::Green))))
        .invalid(
            Style::new()
                .bold()
                .fg_color(Some(Color::Ansi(AnsiColor::Red))),
        )
        .error(
            Style::new()
                .bold()
                .fg_color(Some(Color::Ansi(AnsiColor::Red))),
        )
        .valid(
            Style::new()
                .bold()
                .underline()
                .fg_color(Some(Color::Ansi(AnsiColor::Green))),
        )
        .placeholder(Style::new().fg_color(Some(Color::Ansi(AnsiColor::White))))
}

impl Args {
    pub fn validate(self) -> Result<Args, String> {
        if self.type_query != "nt" && self.type_query != "aa" {
            return Err("Type of query must be either 'nt' or 'aa'".to_string());
        }
        if self.algorithm != 1 && self.algorithm != 2 {
            return Err("Algorithm must be either 1 or 2".to_string());
        }
        if self.reference != "HXB2" && self.reference != "SIVmm239" {
            return Err("Reference genome must be either 'HXB2' or 'SIVmm239'".to_string());
        }
        if self.type_query == "nt" && !self.query.is_empty() {
            let alphabet = alphabets::dna::iupac_alphabet();
            if alphabet.is_word(self.query.as_bytes()) {
                if self.query.len() <= 3 {
                    return Err("Nucleotide sequence length too short".to_string());
                } else {
                    return Ok(self);
                }
            } else {
                return Err("Invalid nucleotide sequence: ".to_string() + &self.query);
            }
        } else if self.type_query == "aa" && !self.query.is_empty() {
            let alphabet = alphabets::protein::iupac_alphabet();
            if alphabet.is_word(self.query.as_bytes()) {
                if self.query.len() <= 3 {
                    return Err("Amino acid sequence length too short".to_string());
                } else {
                    return Ok(self);
                }
            } else {
                return Err("Invalid amino acid sequence: ".to_string() + &self.query);
            }
        }
        Ok(self)
    }
}
