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

## License

The package is available as open source under the terms of the [MIT License](https://opensource.org/licenses/MIT).

## Code of Conduct

Everyone interacting in the VirustLocator project's codebases, issue trackers, chat rooms and mailing lists is expected to follow the [code of conduct](https://github.com/[USERNAME]/virust-locator-ruby/blob/main/CODE_OF_CONDUCT.md).
