//! End-to-end integration tests for viRust-locator binary
//!
//! These tests execute the actual binary and verify its behavior across different scenarios:
//! - Valid nucleotide and amino acid queries
//! - Different reference genomes (HXB2, SIVmm239)
//! - Different algorithms (1, 2)
//! - Multiple queries in a single run
//! - Error cases and edge conditions

use std::path::Path;
use std::process::Command;

/// Helper function to get the path to the binary executable
fn get_binary_path() -> &'static str {
    if cfg!(debug_assertions) {
        "./target/debug/virust-locator"
    } else {
        "./target/release/virust-locator"
    }
}

/// Helper function to run the viRust-locator binary with given arguments
fn run_virust_locator(args: &[&str]) -> (String, String, i32) {
    let output = Command::new(get_binary_path())
        .args(args)
        .output()
        .expect("Failed to execute binary");

    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    let exit_code = output.status.code().unwrap_or(-1);

    (stdout, stderr, exit_code)
}

/// Helper function to parse locator output
/// Expected format: start_pos end_pos similarity reverse_complement query_seq reference_match
fn parse_locator_output(output: &str) -> Option<(i32, i32, i32, bool, String, String)> {
    let line = output.trim();
    if line.is_empty() {
        return None;
    }

    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() != 6 {
        return None;
    }

    Some((
        parts[0].parse().ok()?,
        parts[1].parse().ok()?,
        parts[2].parse().ok()?,
        parts[3].parse().ok()?,
        parts[4].to_string(),
        parts[5].to_string(),
    ))
}

/// Test successful nucleotide query with default parameters (HXB2 reference, algorithm 1)
#[test]
fn test_nucleotide_query_success() {
    let (stdout, stderr, exit_code) = run_virust_locator(&["--query", "ATGCATGCATGC"]);

    assert_eq!(
        exit_code, 0,
        "Binary should exit with code 0 for valid input"
    );
    assert!(
        stderr.is_empty(),
        "No error messages should be printed for valid input"
    );
    assert!(!stdout.is_empty(), "Should produce output for valid query");

    let parsed = parse_locator_output(&stdout);
    assert!(parsed.is_some(), "Should be able to parse output format");

    let (start_pos, end_pos, similarity, _reverse_comp, query_seq, _ref_match) = parsed.unwrap();
    assert!(start_pos > 0, "Start position should be positive");
    assert!(
        end_pos > start_pos,
        "End position should be greater than start position"
    );
    assert!(
        similarity >= 0 && similarity <= 100,
        "Similarity should be between 0 and 100"
    );
    assert_eq!(
        query_seq, "ATGCATGCATGC",
        "Query sequence should match input"
    );
}

/// Test successful amino acid query
#[test]
fn test_amino_acid_query_success() {
    let (stdout, stderr, exit_code) =
        run_virust_locator(&["--query", "MHAC", "--type-query", "aa"]);

    assert_eq!(
        exit_code, 0,
        "Binary should exit with code 0 for valid amino acid input"
    );
    assert!(
        stderr.is_empty(),
        "No error messages should be printed for valid input"
    );
    assert!(
        !stdout.is_empty(),
        "Should produce output for valid amino acid query"
    );

    let parsed = parse_locator_output(&stdout);
    assert!(
        parsed.is_some(),
        "Should be able to parse amino acid output format"
    );

    let (start_pos, end_pos, similarity, _reverse_comp, query_seq, _ref_match) = parsed.unwrap();
    assert!(start_pos > 0, "Start position should be positive");
    assert!(
        end_pos > start_pos,
        "End position should be greater than start position"
    );
    assert!(
        similarity >= 0 && similarity <= 100,
        "Similarity should be between 0 and 100"
    );
    assert_eq!(query_seq, "MHAC", "Query sequence should match input");
}

/// Test with SIVmm239 reference genome
#[test]
fn test_sivmm239_reference() {
    let (stdout, stderr, exit_code) =
        run_virust_locator(&["--query", "ATGCATGCATGC", "--reference", "SIVmm239"]);

    assert_eq!(
        exit_code, 0,
        "Binary should exit with code 0 for SIVmm239 reference"
    );
    assert!(stderr.is_empty(), "No error messages should be printed");
    assert!(
        !stdout.is_empty(),
        "Should produce output for SIVmm239 reference"
    );

    let parsed = parse_locator_output(&stdout);
    assert!(
        parsed.is_some(),
        "Should be able to parse SIVmm239 output format"
    );
}

/// Test with algorithm 2 (fast mode)
#[test]
fn test_algorithm_2() {
    let (stdout, stderr, exit_code) =
        run_virust_locator(&["--query", "ATGCATGCATGC", "--algorithm", "2"]);

    assert_eq!(
        exit_code, 0,
        "Binary should exit with code 0 for algorithm 2"
    );
    assert!(stderr.is_empty(), "No error messages should be printed");
    assert!(!stdout.is_empty(), "Should produce output for algorithm 2");

    let parsed = parse_locator_output(&stdout);
    assert!(
        parsed.is_some(),
        "Should be able to parse algorithm 2 output format"
    );
}

/// Test with multiple queries
#[test]
fn test_multiple_queries() {
    let (stdout, stderr, exit_code) =
        run_virust_locator(&["--query", "ATGCATGCATGC", "GCATGCATGCAT"]);

    assert_eq!(
        exit_code, 0,
        "Binary should exit with code 0 for multiple queries"
    );
    assert!(stderr.is_empty(), "No error messages should be printed");
    assert!(
        !stdout.is_empty(),
        "Should produce output for multiple queries"
    );

    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(
        lines.len(),
        2,
        "Should produce one line of output per query"
    );

    for line in lines {
        let parsed = parse_locator_output(line);
        assert!(parsed.is_some(), "Each line should be parseable");
    }
}

/// Test error case: empty query
#[test]
fn test_error_empty_query() {
    let (stdout, stderr, exit_code) = run_virust_locator(&[]);

    assert_eq!(
        exit_code, 1,
        "Binary should exit with code 1 for empty query"
    );
    assert!(
        stdout.is_empty(),
        "No output should be produced for empty query"
    );
    assert!(
        stderr.contains("Query sequence cannot be empty"),
        "Should show appropriate error message"
    );
}

/// Test error case: sequence too short
#[test]
fn test_error_sequence_too_short() {
    let (stdout, stderr, exit_code) = run_virust_locator(&["--query", "AT"]);

    assert_eq!(
        exit_code, 1,
        "Binary should exit with code 1 for short sequence"
    );
    assert!(
        stdout.is_empty(),
        "No output should be produced for invalid input"
    );
    assert!(
        stderr.contains("Nucleotide sequence length too short"),
        "Should show appropriate error message"
    );
}

/// Test error case: invalid nucleotide sequence
#[test]
fn test_error_invalid_nucleotide() {
    let (stdout, stderr, exit_code) = run_virust_locator(&["--query", "ATGCXYZ"]);

    assert_eq!(
        exit_code, 1,
        "Binary should exit with code 1 for invalid sequence"
    );
    assert!(
        stdout.is_empty(),
        "No output should be produced for invalid input"
    );
    assert!(
        stderr.contains("Invalid nucleotide sequence"),
        "Should show appropriate error message"
    );
}

/// Test error case: invalid amino acid sequence  
#[test]
fn test_error_invalid_amino_acid() {
    let (stdout, stderr, exit_code) =
        run_virust_locator(&["--query", "MHACX123", "--type-query", "aa"]);

    assert_eq!(
        exit_code, 1,
        "Binary should exit with code 1 for invalid amino acid sequence"
    );
    assert!(
        stdout.is_empty(),
        "No output should be produced for invalid input"
    );
    assert!(
        stderr.contains("Invalid amino acid sequence"),
        "Should show appropriate error message"
    );
}

/// Test error case: invalid reference genome
#[test]
fn test_error_invalid_reference() {
    let (stdout, stderr, exit_code) =
        run_virust_locator(&["--query", "ATGCATGCATGC", "--reference", "INVALID"]);

    assert_eq!(
        exit_code, 1,
        "Binary should exit with code 1 for invalid reference"
    );
    assert!(
        stdout.is_empty(),
        "No output should be produced for invalid input"
    );
    assert!(
        stderr.contains("Reference genome must be either"),
        "Should show appropriate error message"
    );
}

/// Test error case: invalid query type
#[test]
fn test_error_invalid_query_type() {
    let (stdout, stderr, exit_code) =
        run_virust_locator(&["--query", "ATGCATGCATGC", "--type-query", "invalid"]);

    assert_eq!(
        exit_code, 1,
        "Binary should exit with code 1 for invalid query type"
    );
    assert!(
        stdout.is_empty(),
        "No output should be produced for invalid input"
    );
    assert!(
        stderr.contains("Type of query must be either"),
        "Should show appropriate error message"
    );
}

/// Test error case: invalid algorithm
#[test]
fn test_error_invalid_algorithm() {
    let (stdout, stderr, exit_code) =
        run_virust_locator(&["--query", "ATGCATGCATGC", "--algorithm", "3"]);

    assert_eq!(
        exit_code, 1,
        "Binary should exit with code 1 for invalid algorithm"
    );
    assert!(
        stdout.is_empty(),
        "No output should be produced for invalid input"
    );
    assert!(
        stderr.contains("Algorithm must be either 1 or 2"),
        "Should show appropriate error message"
    );
}

/// Test help flag
#[test]
fn test_help_flag() {
    let (stdout, stderr, exit_code) = run_virust_locator(&["--help"]);

    assert_eq!(exit_code, 0, "Binary should exit with code 0 for help flag");
    assert!(stderr.is_empty(), "No error messages for help flag");
    assert!(
        stdout.contains("Simple LANL's HIV locator tool"),
        "Should display help text"
    );
    assert!(
        stdout.contains("Usage:"),
        "Should display usage information"
    );
}

/// Test version flag
#[test]
fn test_version_flag() {
    let (stdout, stderr, exit_code) = run_virust_locator(&["--version"]);

    assert_eq!(
        exit_code, 0,
        "Binary should exit with code 0 for version flag"
    );
    assert!(stderr.is_empty(), "No error messages for version flag");
    assert!(!stdout.is_empty(), "Should display version information");
}

/// Test comprehensive scenario with all valid parameters
#[test]
fn test_comprehensive_scenario() {
    let (stdout, stderr, exit_code) = run_virust_locator(&[
        "--query",
        "ATGCATGCATGC",
        "GCATGCATGCAT",
        "--reference",
        "HXB2",
        "--type-query",
        "nt",
        "--algorithm",
        "1",
    ]);

    assert_eq!(
        exit_code, 0,
        "Binary should exit with code 0 for comprehensive test"
    );
    assert!(stderr.is_empty(), "No error messages should be printed");
    assert!(!stdout.is_empty(), "Should produce output");

    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(lines.len(), 2, "Should produce output for both queries");
}

/// Test that binary exists and is executable
#[test]
fn test_binary_exists() {
    let binary_path = get_binary_path();
    assert!(
        Path::new(binary_path).exists(),
        "Binary should exist at {}",
        binary_path
    );

    // Try to run the binary with no args to ensure it's executable
    let output = Command::new(binary_path).output();
    assert!(output.is_ok(), "Binary should be executable");
}
