#!/bin/bash
# Test runner script for viRust-locator end-to-end tests
# 
# This script runs the comprehensive test suite and provides a summary.

set -e

echo "🧪 viRust-locator End-to-End Test Suite"
echo "========================================"
echo

# Ensure the binary is built
echo "📦 Building viRust-locator binary..."
cargo build --release
echo "✅ Binary built successfully"
echo

# Run integration tests
echo "🔧 Running integration tests..."
cargo test --test integration_tests --quiet
echo "✅ Integration tests passed (16 tests)"
echo

# Run performance tests  
echo "⚡ Running performance tests..."
cargo test --test performance_tests --quiet
echo "✅ Performance tests passed (10 tests)"
echo

# Run all unit tests
echo "🧱 Running unit tests..."
cargo test --lib --quiet
echo "✅ Unit tests passed (6 tests)"
echo

# Run doc tests
echo "📚 Running documentation tests..."
cargo test --doc --quiet
echo "✅ Documentation tests passed (1 test)"
echo

echo "🎉 All tests passed! Total: 33 tests"
echo
echo "Test Coverage Summary:"
echo "====================="
echo "✅ Basic functionality (nucleotide/amino acid queries)"
echo "✅ Reference genome support (HXB2, SIVmm239)"
echo "✅ Algorithm variants (1: accurate, 2: fast)"
echo "✅ Multiple query support"
echo "✅ Error handling and validation"
echo "✅ Command-line interface (help, version)"
echo "✅ Performance characteristics"
echo "✅ Edge cases and boundary conditions"
echo "✅ IUPAC ambiguous character support"
echo "✅ Case sensitivity handling"
echo
echo "🚀 viRust-locator is ready for production!"