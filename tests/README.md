# End-to-End Testing for viRust-locator

This directory contains comprehensive end-to-end tests for the viRust-locator binary. The tests verify that the binary works correctly across various scenarios and edge cases.

## Test Files

### `integration_tests.rs`

Main integration tests that cover:

- Basic functionality with nucleotide and amino acid queries
- Different reference genomes (HXB2, SIVmm239)
- Different algorithms (1, 2)
- Multiple queries in single run
- Error handling and validation
- Command-line interface (help, version flags)

### `performance_tests.rs`

Performance and edge case tests that cover:

- Minimum and maximum sequence lengths
- Performance characteristics of different algorithms
- IUPAC ambiguous nucleotides and amino acids
- Case sensitivity handling
- Stress testing with multiple queries
- Very long sequence handling

## Running the Tests

### Run all integration tests:

```bash
cargo test --test integration_tests
```

### Run all performance tests:

```bash
cargo test --test performance_tests
```

### Run all end-to-end tests:

```bash
cargo test --test integration_tests --test performance_tests
```

### Run all tests (including unit tests):

```bash
cargo test
```

### Run tests with output visible:

```bash
cargo test -- --nocapture
```

### Run specific test:

```bash
cargo test test_nucleotide_query_success
```

## Test Structure

Each test follows this pattern:

1. **Setup**: Prepare test input arguments
2. **Execute**: Run the binary with test arguments
3. **Verify**: Check exit code, stdout, and stderr output
4. **Assert**: Validate expected behavior

## Test Helpers

The tests include helper functions:

- `get_binary_path()`: Gets the correct binary path for debug/release builds
- `run_virust_locator()`: Executes the binary and captures output
- `parse_locator_output()`: Parses the standard output format
- `run_virust_locator_timed()`: Measures execution time for performance tests

## Expected Output Format

The binary outputs locator results in this format:

```
start_pos end_pos similarity reverse_complement query_seq reference_match
```

Example:

```
1373    1384    75    false    ATGCATGCATGC    AAGCAGCCATGC
```

## Error Testing

Tests verify proper error handling for:

- Empty queries
- Sequences that are too short (≤3 characters)
- Invalid nucleotide sequences
- Invalid amino acid sequences
- Invalid reference genomes
- Invalid query types
- Invalid algorithm numbers

## Performance Expectations

Tests include performance assertions:

- Individual queries should complete in <5 seconds
- Multiple queries should complete in <10 seconds
- Very long sequences should complete in <30 seconds

## Test Coverage

The end-to-end tests cover:

- ✅ Valid nucleotide queries
- ✅ Valid amino acid queries
- ✅ HXB2 reference genome
- ✅ SIVmm239 reference genome
- ✅ Algorithm 1 (accurate)
- ✅ Algorithm 2 (fast)
- ✅ Multiple queries
- ✅ Error conditions
- ✅ Help and version flags
- ✅ IUPAC ambiguous characters
- ✅ Case sensitivity
- ✅ Performance characteristics
- ✅ Edge cases and boundary conditions

## Continuous Integration

These tests are designed to run in CI environments and will:

- Build the binary before testing
- Test both debug and release builds
- Validate all command-line interface features
- Ensure robust error handling
- Verify performance characteristics

## Adding New Tests

When adding new tests:

1. Place in appropriate test file (`integration_tests.rs` for functionality, `performance_tests.rs` for performance)
2. Follow naming convention: `test_<feature>_<scenario>`
3. Include proper assertions for exit code, stdout, and stderr
4. Add documentation comments explaining the test purpose
5. Consider both success and failure cases
