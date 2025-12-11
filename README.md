# viRust-Locator

A simplified version of LANL HIV Database HIV Locator tool

## Usage

Usage: `cargo run -- [OPTIONS] --query`

Options:

-q, --query Query sequence

-r, --reference Reference genome, either HXB2 or SIVmm239 [default: HXB2]

-t, --type-query <TYPE_QUERY> Type of query, either nt or aa [default: nt]

-a, --algorithm algorithm for locator, 1 is accurate but slower, 2 is fast but less accurate, suitable for smaller query sequences [default: 1]

-h, --help Print help

-V, --version Print version

## Examples

```bash
# Basic nucleotide query
cargo run -- --query "ATGCATGCATGC"

# Amino acid query with SIVmm239 reference
cargo run -- --query "MHAC" --type-query "aa" --reference "SIVmm239"

# Multiple queries with fast algorithm
cargo run -- --query "ATGCATGC" "GCATGCAT" --algorithm "2"
```

## Testing

This project includes comprehensive end-to-end tests for the binary:

### Run all tests:

```bash
cargo test
```

### Run integration tests only:

```bash
cargo test --test integration_tests
```

### Run performance tests only:

```bash
cargo test --test performance_tests
```

### Run the complete test suite:

```bash
./run_tests.sh
```

The test suite includes:

- ✅ Basic functionality (nucleotide/amino acid queries)
- ✅ Reference genome support (HXB2, SIVmm239)
- ✅ Algorithm variants (1: accurate, 2: fast)
- ✅ Multiple query support
- ✅ Error handling and validation
- ✅ Command-line interface (help, version)
- ✅ Performance characteristics
- ✅ Edge cases and boundary conditions
- ✅ IUPAC ambiguous character support
- ✅ Case sensitivity handling

See the [tests/README.md](tests/README.md) for detailed testing documentation.

## License

The package is available as open source under the terms of the [MIT License](https://opensource.org/licenses/MIT).

## Code of Conduct

Everyone interacting in the VirustLocator project's codebases, issue trackers, chat rooms and mailing lists is expected to follow the [code of conduct](https://github.com/[USERNAME]/virust-locator-ruby/blob/main/CODE_OF_CONDUCT.md).
