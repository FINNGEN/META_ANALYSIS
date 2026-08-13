#!/bin/bash

# Simple bash test runner for R scripts
# This provides a quick way to verify tests work without running full test suite

echo "============================================"
echo "R Scripts Test Verification"
echo "============================================"
echo ""

cd "$(dirname "$0")"

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R."
    exit 1
fi

echo "✓ Rscript found"
echo ""

# Check that test files exist
test_files=("test_miami.R" "test_qqplot.R" "test_qc.R")
for test_file in "${test_files[@]}"; do
    if [ -f "$test_file" ]; then
        echo "✓ Found: $test_file"
    else
        echo "✗ Missing: $test_file"
        exit 1
    fi
done

echo ""
echo "All test files found."
echo ""
echo "To run the full test suite:"
echo "  cd scripts/tests"
echo "  Rscript run_all_tests.R"
echo ""
echo "To run individual tests:"
echo "  cd scripts/tests"
echo "  Rscript test_miami.R"
echo "  Rscript test_qqplot.R"
echo "  Rscript test_qc.R"
echo ""
