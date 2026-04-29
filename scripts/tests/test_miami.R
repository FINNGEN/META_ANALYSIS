#!/usr/bin/env Rscript

# Comprehensive tests for miami.R script
# Tests that the script runs successfully with various reasonable inputs
# Note: We don't validate plot quality, just that the script completes and outputs exist

library(data.table)

# Helper function to create test data
create_test_sumstats <- function(filename, n_variants = 1000, pvalue_type = "p") {
  set.seed(42)
  
  # Create chromosomes 1-22 and X
  chrs <- c(rep(1:22, each = floor(n_variants / 23)), rep("X", n_variants %% 23))
  
  data <- data.table(
    `#CHR` = chrs,
    POS = sample(1:250000000, n_variants, replace = TRUE),
    pval1 = runif(n_variants, 0, 1),
    pval2 = runif(n_variants, 0, 1)
  )
  
  # Add some significant hits
  sig_indices <- sample(1:n_variants, min(50, n_variants / 10))
  data[sig_indices, pval1 := runif(.N, 0, 5e-8)]
  data[sig_indices, pval2 := runif(.N, 0, 5e-8)]
  
  if (pvalue_type == "mlogp") {
    # Convert to -log10 p-values
    data[, pval1 := -log10(pval1)]
    data[, pval2 := -log10(pval2)]
  }
  
  fwrite(data, filename, sep = "\t", quote = FALSE)
  return(data)
}

# Helper function to clean up test files
cleanup_files <- function(pattern) {
  files_to_remove <- list.files(pattern = pattern, full.names = TRUE)
  if (length(files_to_remove) > 0) {
    file.remove(files_to_remove)
  }
}

# Test 1: Basic functionality with p-values
test_miami_basic_pvalues <- function() {
  cat("\n=== Test 1: Basic functionality with p-values ===\n")
  
  test_file <- "data/test_miami_basic.tsv"
  output_prefix <- "test_miami_basic_output"
  
  # Cleanup any existing test files
  cleanup_files("test_miami_basic")
  
  # Create test data
  cat("Creating test data...\n")
  create_test_sumstats(test_file, n_variants = 1000, pvalue_type = "p")
  
  # Run miami.R
  cat("Running miami.R...\n")
  cmd <- sprintf(
    "Rscript ../miami.R --file %s --out %s --pval_cols pval1,pval2 --pvalue_type p",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("miami.R failed with exit code ", result)
  }
  
  # Check that output file exists
  expected_output <- paste0(output_prefix, "_pval1_miami.png")
  if (!file.exists(expected_output)) {
    stop("Expected output file not created: ", expected_output)
  }
  
  # Check file size > 0
  file_size <- file.info(expected_output)$size
  if (file_size == 0) {
    stop("Output file is empty: ", expected_output)
  }
  
  cat("✓ Test passed: Output file created with size", file_size, "bytes\n")
  
  # Cleanup
  cleanup_files("test_miami_basic")
  
  return(TRUE)
}

# Test 2: With -log10 p-values
test_miami_mlogp <- function() {
  cat("\n=== Test 2: With -log10 p-values ===\n")
  
  test_file <- "test_miami_mlogp.tsv"
  output_prefix <- "test_miami_mlogp_output"
  
  # Cleanup any existing test files
  cleanup_files("test_miami_mlogp")
  
  # Create test data with -log10 p-values
  cat("Creating test data with -log10 p-values...\n")
  create_test_sumstats(test_file, n_variants = 1000, pvalue_type = "mlogp")
  
  # Run miami.R
  cat("Running miami.R...\n")
  cmd <- sprintf(
    "Rscript ../miami.R --file %s --out %s --pval_cols pval1,pval2 --pvalue_type mlogp",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("miami.R failed with exit code ", result)
  }
  
  # Check that output file exists
  expected_output <- paste0(output_prefix, "_pval1_miami.png")
  if (!file.exists(expected_output)) {
    stop("Expected output file not created: ", expected_output)
  }
  
  # Check file size > 0
  file_size <- file.info(expected_output)$size
  if (file_size == 0) {
    stop("Output file is empty: ", expected_output)
  }
  
  cat("✓ Test passed: Output file created with size", file_size, "bytes\n")
  
  # Cleanup
  cleanup_files("test_miami_mlogp")
  
  return(TRUE)
}

# Test 3: Auto-detect p-value type
test_miami_auto_detect <- function() {
  cat("\n=== Test 3: Auto-detect p-value type ===\n")
  
  test_file <- "test_miami_auto.tsv"
  output_prefix <- "test_miami_auto_output"
  
  # Cleanup any existing test files
  cleanup_files("test_miami_auto")
  
  # Create test data with p-values (will be auto-detected)
  cat("Creating test data...\n")
  create_test_sumstats(test_file, n_variants = 1000, pvalue_type = "p")
  
  # Run miami.R with detect mode
  cat("Running miami.R with auto-detect...\n")
  cmd <- sprintf(
    "Rscript miami.R --file %s --out %s --pval_cols pval1,pval2 --pvalue_type detect",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("miami.R failed with exit code ", result)
  }
  
  # Check that output file exists
  expected_output <- paste0(output_prefix, "_pval1_miami.png")
  if (!file.exists(expected_output)) {
    stop("Expected output file not created: ", expected_output)
  }
  
  cat("✓ Test passed: Auto-detection worked\n")
  
  # Cleanup
  cleanup_files("test_miami_auto")
  
  return(TRUE)
}

# Test 4: With highlight option
test_miami_highlight <- function() {
  cat("\n=== Test 4: With highlight option ===\n")
  
  test_file <- "test_miami_highlight.tsv"
  output_prefix <- "test_miami_highlight_output"
  
  # Cleanup any existing test files
  cleanup_files("test_miami_highlight")
  
  # Create test data
  cat("Creating test data...\n")
  create_test_sumstats(test_file, n_variants = 1000, pvalue_type = "p")
  
  # Run miami.R with highlight
  cat("Running miami.R with highlight...\n")
  cmd <- sprintf(
    "Rscript ../miami.R --file %s --out %s --pval_cols pval1,pval2 --pvalue_type p --highlight",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("miami.R failed with exit code ", result)
  }
  
  # Check that output file exists
  expected_output <- paste0(output_prefix, "_pval1_miami.png")
  if (!file.exists(expected_output)) {
    stop("Expected output file not created: ", expected_output)
  }
  
  cat("✓ Test passed: Highlight mode worked\n")
  
  # Cleanup
  cleanup_files("test_miami_highlight")
  
  return(TRUE)
}

# Test 5: Custom column names
test_miami_custom_columns <- function() {
  cat("\n=== Test 5: Custom column names ===\n")
  
  test_file <- "test_miami_custom.tsv"
  output_prefix <- "test_miami_custom_output"
  
  # Cleanup any existing test files
  cleanup_files("test_miami_custom")
  
  # Create test data with custom column names
  cat("Creating test data with custom columns...\n")
  set.seed(42)
  n_variants <- 1000
  chrs <- c(rep(1:22, each = floor(n_variants / 23)), rep("X", n_variants %% 23))
  
  data <- data.table(
    chromosome = chrs,
    position = sample(1:250000000, n_variants, replace = TRUE),
    p_meta = runif(n_variants, 0, 1),
    p_study = runif(n_variants, 0, 1)
  )
  
  # Add some significant hits
  sig_indices <- sample(1:n_variants, min(50, n_variants / 10))
  data[sig_indices, p_meta := runif(.N, 0, 5e-8)]
  data[sig_indices, p_study := runif(.N, 0, 5e-8)]
  
  fwrite(data, test_file, sep = "\t", quote = FALSE)
  
  # Run miami.R with custom column names
  cat("Running miami.R with custom columns...\n")
  cmd <- sprintf(
    "Rscript ../miami.R --file %s --out %s --chr_col chromosome --pos_col position --pval_cols p_meta,p_study --pvalue_type p",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("miami.R failed with exit code ", result)
  }
  
  # Check that output file exists
  expected_output <- paste0(output_prefix, "_p_meta_miami.png")
  if (!file.exists(expected_output)) {
    stop("Expected output file not created: ", expected_output)
  }
  
  cat("✓ Test passed: Custom columns worked\n")
  
  # Cleanup
  cleanup_files("test_miami_custom")
  
  return(TRUE)
}

# Main test runner
main <- function() {
  cat("========================================\n")
  cat("Running miami.R test suite\n")
  cat("========================================\n")
  
  # Change to scripts directory if not already there
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("--file=", args, value = TRUE)
  if (length(file_flag) > 0) {
    script_dir <- dirname(normalizePath(sub("--file=", "", file_flag)))
    setwd(file.path(script_dir, ".."))
  }
  
  tests <- list(
    test_miami_basic_pvalues,
    test_miami_mlogp,
    test_miami_auto_detect,
    test_miami_highlight,
    test_miami_custom_columns
  )
  
  passed <- 0
  failed <- 0
  
  for (test in tests) {
    tryCatch({
      test()
      passed <- passed + 1
    }, error = function(e) {
      cat("✗ Test failed:", conditionMessage(e), "\n")
      failed <- failed + 1
    })
  }
  
  cat("\n========================================\n")
  cat("Test Results:\n")
  cat(sprintf("Passed: %d\n", passed))
  cat(sprintf("Failed: %d\n", failed))
  cat("========================================\n")
  
  if (failed > 0) {
    quit(status = 1)
  }
}

# Run tests if script is executed directly
if (!interactive()) {
  main()
}
