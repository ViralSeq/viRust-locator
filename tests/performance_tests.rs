//! Performance and edge case tests for viRust-locator binary
//!
//! These tests focus on:
//! - Performance characteristics
//! - Edge cases with boundary conditions
//! - Stress testing with various sequence lengths
//! - Memory usage validation

use std::process::Command;
use std::time::Instant;

/// Helper function to get the path to the binary executable
fn get_binary_path() -> &'static str {
    if cfg!(debug_assertions) {
        "./target/debug/virust-locator"
    } else {
        "./target/release/virust-locator"
    }
}

/// Helper function to run the viRust-locator binary with given arguments and measure time
fn run_virust_locator_timed(args: &[&str]) -> (String, String, i32, std::time::Duration) {
    let start = Instant::now();
    let output = Command::new(get_binary_path())
        .args(args)
        .output()
        .expect("Failed to execute binary");
    let duration = start.elapsed();

    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    let exit_code = output.status.code().unwrap_or(-1);

    (stdout, stderr, exit_code, duration)
}

/// Test with minimum length sequence (4 nucleotides)
#[test]
fn test_minimum_length_sequence() {
    let (stdout, stderr, exit_code, _) = run_virust_locator_timed(&["--query", "ATGC"]);

    assert_eq!(exit_code, 0, "Should accept minimum length sequence");
    assert!(
        stderr.is_empty(),
        "No error messages for minimum length sequence"
    );
    assert!(
        !stdout.is_empty(),
        "Should produce output for minimum length sequence"
    );
}

/// Test with longer sequences to check performance
#[test]
fn test_long_sequence_performance() {
    // Create a longer sequence (120 nucleotides)
    let long_sequence = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC";

    let (stdout, stderr, exit_code, duration) =
        run_virust_locator_timed(&["--query", long_sequence]);

    assert_eq!(exit_code, 0, "Should handle long sequences");
    assert!(stderr.is_empty(), "No error messages for long sequences");
    assert!(
        !stdout.is_empty(),
        "Should produce output for long sequences"
    );

    // Performance should be reasonable (less than 5 seconds for a sequence this size)
    assert!(
        duration.as_secs() < 5,
        "Should complete in reasonable time: {:?}",
        duration
    );
}

/// Test with sequences containing IUPAC ambiguous nucleotides
#[test]
fn test_ambiguous_nucleotides() {
    // Test with IUPAC ambiguous nucleotides (R, Y, K, M, S, W, B, D, H, V, N)
    let ambiguous_sequence = "ATGCRYWKMSBDHVN";

    let (stdout, stderr, exit_code, _) = run_virust_locator_timed(&["--query", ambiguous_sequence]);

    assert_eq!(exit_code, 0, "Should handle IUPAC ambiguous nucleotides");
    assert!(
        stderr.is_empty(),
        "No error messages for ambiguous nucleotides"
    );
    assert!(
        !stdout.is_empty(),
        "Should produce output for ambiguous nucleotides"
    );
}

/// Test with sequences containing ambiguous amino acids
#[test]
fn test_ambiguous_amino_acids() {
    // Test with valid ambiguous amino acids (B, Z, X are valid in IUPAC)
    let ambiguous_aa = "MHACBZX";

    let (stdout, stderr, exit_code, _) =
        run_virust_locator_timed(&["--query", ambiguous_aa, "--type-query", "aa"]);

    assert_eq!(exit_code, 0, "Should handle ambiguous amino acids");
    assert!(
        stderr.is_empty(),
        "No error messages for ambiguous amino acids"
    );
    assert!(
        !stdout.is_empty(),
        "Should produce output for ambiguous amino acids"
    );
}

/// Test performance difference between algorithms
#[test]
fn test_algorithm_performance_comparison() {
    let test_sequence = "ATGCATGCATGCATGCATGCATGC";

    // Test algorithm 1 (accurate but slower)
    let (_, _, exit_code1, duration1) =
        run_virust_locator_timed(&["--query", test_sequence, "--algorithm", "1"]);

    // Test algorithm 2 (fast but less accurate)
    let (_, _, exit_code2, duration2) =
        run_virust_locator_timed(&["--query", test_sequence, "--algorithm", "2"]);

    assert_eq!(exit_code1, 0, "Algorithm 1 should succeed");
    assert_eq!(exit_code2, 0, "Algorithm 2 should succeed");

    // Note: Algorithm 2 should generally be faster, but for very short sequences
    // the difference might not be significant. We just ensure both complete in reasonable time.
    assert!(
        duration1.as_secs() < 10,
        "Algorithm 1 should complete in reasonable time"
    );
    assert!(
        duration2.as_secs() < 10,
        "Algorithm 2 should complete in reasonable time"
    );
}

/// Test multiple long sequences
#[test]
fn test_multiple_long_sequences() {
    let seq1 = "ATGCATGCATGCATGCATGCATGCATGC";
    let seq2 = "GCATGCATGCATGCATGCATGCATGCAT";
    let seq3 = "CATGCATGCATGCATGCATGCATGCATG";

    let (stdout, stderr, exit_code, duration) =
        run_virust_locator_timed(&["--query", seq1, seq2, seq3]);

    assert_eq!(exit_code, 0, "Should handle multiple long sequences");
    assert!(
        stderr.is_empty(),
        "No error messages for multiple long sequences"
    );
    assert!(
        !stdout.is_empty(),
        "Should produce output for multiple long sequences"
    );

    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(lines.len(), 3, "Should produce three lines of output");

    // Should complete in reasonable time even with multiple sequences
    assert!(
        duration.as_secs() < 10,
        "Should complete multiple sequences in reasonable time: {:?}",
        duration
    );
}

/// Test case sensitivity (sequences should be case-insensitive)
#[test]
fn test_case_insensitive_sequences() {
    let lowercase_seq = "atgcatgcatgc";
    let uppercase_seq = "ATGCATGCATGC";
    let mixed_case_seq = "AtGcAtGcAtGc";

    let (stdout1, stderr1, exit_code1, _) = run_virust_locator_timed(&["--query", lowercase_seq]);
    let (stdout2, stderr2, exit_code2, _) = run_virust_locator_timed(&["--query", uppercase_seq]);
    let (stdout3, stderr3, exit_code3, _) = run_virust_locator_timed(&["--query", mixed_case_seq]);

    assert_eq!(exit_code1, 0, "Should handle lowercase sequences");
    assert_eq!(exit_code2, 0, "Should handle uppercase sequences");
    assert_eq!(exit_code3, 0, "Should handle mixed case sequences");

    assert!(
        stderr1.is_empty() && stderr2.is_empty() && stderr3.is_empty(),
        "No error messages for case variations"
    );
    assert!(
        !stdout1.is_empty() && !stdout2.is_empty() && !stdout3.is_empty(),
        "Should produce output for all case variations"
    );
}

/// Test edge case: sequence at maximum typical length
#[test]
fn test_maximum_typical_length() {
    // HIV genome is approximately 9000-10000 base pairs, test with a significant portion
    let long_seq = "A".repeat(1000) + &"T".repeat(500) + &"G".repeat(500) + &"C".repeat(500);

    let (stdout, stderr, exit_code, duration) = run_virust_locator_timed(&["--query", &long_seq]);

    assert_eq!(exit_code, 0, "Should handle very long sequences");
    assert!(
        stderr.is_empty(),
        "No error messages for very long sequences"
    );
    assert!(
        !stdout.is_empty(),
        "Should produce output for very long sequences"
    );

    // Even very long sequences should complete in reasonable time
    assert!(
        duration.as_secs() < 30,
        "Should complete very long sequence in reasonable time: {:?}",
        duration
    );
}

/// Test both reference genomes with same sequence to compare results
#[test]
fn test_reference_genome_comparison() {
    let test_sequence = "ATGCATGCATGCATGC";

    let (stdout_hxb2, stderr_hxb2, exit_code_hxb2, _) =
        run_virust_locator_timed(&["--query", test_sequence, "--reference", "HXB2"]);

    let (stdout_siv, stderr_siv, exit_code_siv, _) =
        run_virust_locator_timed(&["--query", test_sequence, "--reference", "SIVmm239"]);

    assert_eq!(exit_code_hxb2, 0, "HXB2 reference should work");
    assert_eq!(exit_code_siv, 0, "SIVmm239 reference should work");

    assert!(
        stderr_hxb2.is_empty() && stderr_siv.is_empty(),
        "No errors for either reference"
    );
    assert!(
        !stdout_hxb2.is_empty() && !stdout_siv.is_empty(),
        "Both references should produce output"
    );

    // Results should be different for different references (different start/end positions)
    assert_ne!(
        stdout_hxb2, stdout_siv,
        "Different references should produce different results"
    );
}

/// Test stress scenario with many small queries
#[test]
fn test_stress_many_small_queries() {
    let mut args = vec!["--query"];
    let queries = vec!["ATGC", "GCTA", "TACG", "CGAT", "TGCA"];

    for query in &queries {
        args.push(query);
    }

    let (stdout, stderr, exit_code, duration) = run_virust_locator_timed(&args);

    assert_eq!(exit_code, 0, "Should handle many small queries");
    assert!(
        stderr.is_empty(),
        "No error messages for many small queries"
    );
    assert!(
        !stdout.is_empty(),
        "Should produce output for many small queries"
    );

    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(
        lines.len(),
        queries.len(),
        "Should produce output for each query"
    );

    // Should complete multiple small queries quickly
    assert!(
        duration.as_secs() < 5,
        "Should complete many small queries quickly: {:?}",
        duration
    );
}
