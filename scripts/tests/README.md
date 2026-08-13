# Meta-Analysis Pipeline Test Suite

This directory contains comprehensive tests for both Python and R scripts used in the meta-analysis pipeline.

## Overview

The test suite includes:

### Python Scripts
- **harmonize.py** - Harmonize GWAS summary statistics with numeric precision validation
- **munge.py** - Process and format GWAS files with rounding precision tests  
- **meta_analysis.py** - Meta-analysis functions and statistical methods
- **test_precision_integration.py** - End-to-end integration tests for numeric precision

### R Scripts
- **miami.R** - Miami plot generation
- **qqplot.R** - QQ plot and Manhattan plot generation
- **qc.R** - Quality control plots and metrics for meta-analysis

## Testing Philosophy

The tests focus on **functional verification** - ensuring the scripts run successfully with reasonable inputs and produce expected output files. They:

✓ Create realistic test data (GWAS summary statistics)
✓ Run each script with various parameter combinations
✓ Verify output files are created
✓ Verify output files are non-empty
✓ Validate numeric precision preservation
✓ Automatically clean up test files

The tests **do NOT** validate:
- Plot visual quality
- Statistical correctness of complex calculations (beyond basic validation)
- Exact output values (except for precision tests)

## Test Coverage

| Script | Test Cases | Coverage |
|--------|------------|----------|
| **Python Scripts** | | |
| harmonize.py | 15+ | Precision preservation, beta/SE rounding, strand flipping |
| munge.py | 10+ | Rounding precision, column handling, data validation |
| meta_analysis.py | 50+ | All meta-analysis methods, heterogeneity, formatting |
| precision_integration.py | 5+ | End-to-end precision validation across pipeline |
| **R Scripts** | | |
| miami.R | 5 | All major options and input types |
| qqplot.R | 5 | All major options and multiple p-value columns |
| qc.R | 5 | Meta-analysis QC, LOO, weighted regression |

## Test Data

Test data files are stored in the `data/` subdirectory:
- `test_miami_basic.tsv` - Sample GWAS data for Miami plot tests

Each test creates additional temporary test data files that are automatically cleaned up after the test completes. Test data includes:
- Realistic GWAS summary statistics with multiple variants
- Significant hits (p < 5e-8) to test plot generation
- Multiple studies for meta-analysis tests
- Config JSON files for qc.R tests

---

## Python Tests

### Test Files

- **test_harmonize.py** - Tests for harmonize.py
  - Decimal place counting
  - Beta/SE precision preservation when flipping alleles
  - Frequency rounding after strand flipping
  - EAF validation
  - Numeric precision throughout harmonization

- **test_munge.py** - Tests for munge.py
  - Default rounding precision (6 decimals)
  - Beta value rounding
  - Standard error rounding
  - Frequency rounding
  - Column name handling

- **test_meta_analysis.py** - Tests for meta_analysis.py
  - Heterogeneity tests (no heterogeneity, with heterogeneity)
  - Sample size-based meta-analysis
  - Inverse variance meta-analysis
  - Variance weight meta-analysis
  - Multiple studies handling
  - Edge cases (single study, missing data, zero weights)
  - Numeric formatting

- **test_precision_integration.py** - Integration tests
  - End-to-end precision validation from munge → harmonize → meta-analysis
  - Validates that 6-decimal precision is maintained throughout pipeline
  - Tests realistic workflow scenarios

### Requirements

Python 3.6+ and the following packages:
- numpy
- scipy
- pandas (for some tests)
- Standard library: unittest, tempfile, gzip

### Running Python Tests

**Run all Python tests:**
```bash
cd scripts/tests
python3 -m unittest discover -s . -p "test_*.py"
```

**Run individual test files:**
```bash
cd scripts/tests
python3 test_harmonize.py
python3 test_munge.py
python3 test_meta_analysis.py
python3 test_precision_integration.py
```

**Run specific test class or method:**
```bash
cd scripts/tests
python3 -m unittest test_meta_analysis.TestMetaAnalysis.test_het_test_no_heterogeneity
```

### Python Test Framework

The Python tests use the standard `unittest` framework. Key features:

- **Import mechanism**: Tests add the parent directory to `sys.path` to import the scripts:
  ```python
  sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
  import harmonize
  ```

- **Test structure**: Each test file contains one or more test classes inheriting from `unittest.TestCase`

- **Assertions**: Standard unittest assertions (`assertEqual`, `assertAlmostEqual`, `assertTrue`, etc.)

- **Test isolation**: Tests are independent and don't share state

### Python Test Output

**Successful run:**
```
......................................................
----------------------------------------------------------------------
Ran 54 tests in 0.123s

OK
```

**Failed test:**
```
F.....
======================================================================
FAIL: test_beta_precision_6_decimals (test_munge.TestMungePrecision)
----------------------------------------------------------------------
AssertionError: Beta should be rounded to 6 decimal places

----------------------------------------------------------------------
Ran 6 tests in 0.045s

FAILED (failures=1)
```

---

## R Tests

### Test Files

- **test_miami.R** - Tests for miami.R
  - Basic functionality with p-values
  - With -log10 p-values
  - Auto-detect p-value type
  - Highlight option
  - Custom column names

- **test_qqplot.R** - Tests for qqplot.R
  - Basic functionality
  - Multiple p-value columns
  - Minrep column option
  - Custom loglog parameters
  - Lambda calculation validation

- **test_qc.R** - Tests for qc.R
  - Basic functionality
  - Leave-one-out analysis
  - Weighted regression
  - Custom p-value thresholds
  - Keep HLA region variants

### Master Test Runner

- **run_all_tests.R** - Runs all R test suites and provides a summary
- **verify_tests.sh** - Quick verification script to check test environment

### Requirements

The following R packages are required:
- data.table
- qqman
- optparse
- R.utils
- ggplot2
- rjson
- stringi
- ggpubr

These packages will be automatically installed by the scripts if not already present.

### Running R Tests

**Run all R test suites:**

```bash
cd scripts/tests
Rscript run_all_tests.R
```

**Run individual test suites:**

```bash
cd scripts/tests
Rscript test_miami.R
Rscript test_qqplot.R
Rscript test_qc.R
```

**Quick verification:**

```bash
cd scripts/tests
./verify_tests.sh
```

**Run a single test function:**

Each test file can also be sourced in an R session to run individual tests:

```r
# From scripts/tests directory
source("test_miami.R")

# Run individual test
test_miami_basic_pvalues()
```

### R Test Output

When tests pass, you'll see output like:

```
========================================
Running miami.R test suite
========================================

=== Test 1: Basic functionality with p-values ===
Creating test data...
Running miami.R...
✓ Test passed: Output file created with size 12345 bytes

...

========================================
Test Results:
Passed: 5
Failed: 0
========================================
```

### Failed Test Run

If a test fails, you'll see error messages indicating what went wrong:

```
✗ Test failed: Expected output file not created: test_output.png
```

---

## Exit Codes

Both Python and R test suites return standard exit codes for CI/CD integration:

- **0**: All tests passed
- **1**: One or more tests failed

Example usage in CI/CD:
```bash
cd scripts/tests

# Run R tests
Rscript run_all_tests.R
echo $?  # Returns 0 if passed, 1 if failed

# Run Python tests
python3 -m unittest discover -s . -p "test_*.py"
echo $?  # Returns 0 if passed, 1 if failed
```

---

## Debugging Failed Tests

If a test fails:

1. Check that all required R packages are installed
2. Look at the error message to identify which output file is missing
3. Try running the R script manually with the test parameters
4. Check that you're running from the correct directory
5. Ensure the R scripts have execute permissions

## Adding New Tests

### Adding R Tests

To add new R tests:

1. Choose the appropriate test file (test_miami.R, test_qqplot.R, or test_qc.R)
2. Add a new test function following the pattern:
   ```r
   test_<script>_<feature> <- function() {
     cat("\n=== Test N: Description ===\n")
     
     # Create test data
     # Run script
     # Check outputs
     # Cleanup
     
     return(TRUE)
   }
   ```
3. Add the test function to the `tests` list in the `main()` function
4. Run the test suite to verify it works

### Adding Python Tests

To add new Python tests:

1. Choose the appropriate test file (test_harmonize.py, test_munge.py, etc.) or create a new one
2. Add a new test method to the test class:
   ```python
   def test_new_feature(self):
       """Test description"""
       # Setup test data
       # Run function
       # Assert expected results
       self.assertEqual(expected, actual)
   ```
3. Run the test file to verify it works

## Continuous Integration

These tests can be integrated into CI/CD pipelines. Both Python and R tests return appropriate exit codes:

```bash
cd scripts/tests

# Run R tests
Rscript run_all_tests.R
if [ $? -ne 0 ]; then echo "R tests failed"; exit 1; fi

# Run Python tests
python3 -m unittest discover -s . -p "test_*.py"
if [ $? -ne 0 ]; then echo "Python tests failed"; exit 1; fi
```

## Notes

- Tests create temporary files in the scripts/tests directory
- Test data files are stored in scripts/tests/data/ subdirectory
- All temporary files are cleaned up after each test
- Tests run independently and don't affect each other
- R tests: Plot generation is tested but plot quality is not validated
- R tests: Tests use fixed random seeds for reproducibility
- Python tests: Precision validation ensures 6-decimal accuracy throughout pipeline

## Troubleshooting

### R Package Installation Issues

If R packages fail to install automatically, install them manually:

```r
install.packages(c("data.table", "qqman", "optparse", "R.utils", 
                   "ggplot2", "rjson", "stringi", "ggpubr"))
```

### Python Package Installation Issues

If Python packages are missing:

```bash
pip install numpy scipy pandas
```

### Permission Errors

Make sure the test scripts are executable:

```bash
chmod +x scripts/tests/run_all_tests.R
chmod +x scripts/tests/verify_tests.sh
chmod +x scripts/tests/*.R
```

### Memory Issues

If R tests fail due to memory, you can modify the test data size by editing the `n_variants` parameter in the test data creation functions.

### Python Import Errors

If Python tests cannot import the scripts, verify you're running from the correct directory:

```bash
cd scripts/tests
python3 test_harmonize.py
```

The tests use `sys.path.insert()` to import from the parent directory.

## Contact

For issues or questions about the tests, please contact the development team.
