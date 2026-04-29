#!/usr/bin/env Rscript

# Comprehensive tests for qqplot.R script
# Tests that the script runs successfully with various reasonable inputs
# Note: We don't validate plot quality, just that the script completes and outputs exist

library(data.table)

# Helper function to create test data
create_test_sumstats <- function(filename, n_variants = 10000) {
  set.seed(42)
  
  # Create chromosomes 1-22 and X
  chrs <- c(rep(1:22, each = floor(n_variants / 23)), rep("X", n_variants %% 23))
  
  data <- data.table(
    CHR = chrs,
    BP = sample(1:250000000, n_variants, replace = TRUE),
    P = runif(n_variants, 0, 1)
  )
  
  # Add some significant hits (realistic GWAS distribution)
  # Most p-values are null, some are significant
  sig_indices <- sample(1:n_variants, min(100, n_variants / 100))
  data[sig_indices, P := runif(.N, 0, 5e-8)]
  
  # Add some moderately significant
  mod_indices <- sample(1:n_variants, min(500, n_variants / 20))
  data[mod_indices, P := runif(.N, 0, 1e-4)]
  
  fwrite(data, filename, sep = "\t", quote = FALSE)
  return(data)
}

# Helper function to create test data with multiple p-value columns
create_test_sumstats_multiple_pvals <- function(filename, n_variants = 10000) {
  set.seed(42)
  
  chrs <- c(rep(1:22, each = floor(n_variants / 23)), rep("X", n_variants %% 23))
  
  data <- data.table(
    CHR = chrs,
    BP = sample(1:250000000, n_variants, replace = TRUE),
    P1 = runif(n_variants, 0, 1),
    P2 = runif(n_variants, 0, 1),
    P3 = runif(n_variants, 0, 1)
  )
  
  # Add some significant hits to each
  for (pcol in c("P1", "P2", "P3")) {
    sig_indices <- sample(1:n_variants, min(50, n_variants / 200))
    data[sig_indices, (pcol) := runif(.N, 0, 5e-8)]
  }
  
  fwrite(data, filename, sep = "\t", quote = FALSE)
  return(data)
}

# Helper function to create test data with minrep format
create_test_sumstats_minrep <- function(filename, n_variants = 10000) {
  set.seed(42)
  
  chrs <- c(rep(1:22, each = floor(n_variants / 23)), rep("X", n_variants %% 23))
  positions <- sample(1:250000000, n_variants, replace = TRUE)
  
  data <- data.table(
    minrep_id = paste(chrs, positions, "A", "T", sep = ":"),
    P = runif(n_variants, 0, 1)
  )
  
  # Add some significant hits
  sig_indices <- sample(1:n_variants, min(100, n_variants / 100))
  data[sig_indices, P := runif(.N, 0, 5e-8)]
  
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

# Test 1: Basic functionality with default options
test_qqplot_basic <- function() {
  cat("\n=== Test 1: Basic functionality with default columns ===\n")
  
  test_file <- "test_qqplot_basic.tsv"
  output_prefix <- "test_qqplot_basic_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qqplot_basic")
  
  # Create test data
  cat("Creating test data...\n")
  create_test_sumstats(test_file, n_variants = 10000)
  
  # Run qqplot.R
  cat("Running qqplot.R...\n")
  cmd <- sprintf(
    "Rscript ../qqplot.R --file %s --out %s --chrcol CHR --bp_col BP --pval_col P",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qqplot.R failed with exit code ", result)
  }
  
  # Check that output files exist
  expected_outputs <- c(
    paste0(output_prefix, "_P_qqplot.png"),
    paste0(output_prefix, "_P_manhattan.png"),
    paste0(output_prefix, "_P_manhattan_loglog.png"),
    paste0(output_prefix, "_P_qquantiles.txt")
  )
  
  for (output in expected_outputs) {
    if (!file.exists(output)) {
      stop("Expected output file not created: ", output)
    }
    
    # Check file size > 0
    file_size <- file.info(output)$size
    if (file_size == 0) {
      stop("Output file is empty: ", output)
    }
    cat("✓ Created:", output, "(", file_size, "bytes)\n")
  }
  
  cat("✓ Test passed: All output files created\n")
  
  # Cleanup
  cleanup_files("test_qqplot_basic")
  
  return(TRUE)
}

# Test 2: Multiple p-value columns
test_qqplot_multiple_pvals <- function() {
  cat("\n=== Test 2: Multiple p-value columns ===\n")
  
  test_file <- "test_qqplot_multi.tsv"
  output_prefix <- "test_qqplot_multi_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qqplot_multi")
  
  # Create test data with multiple p-value columns
  cat("Creating test data with multiple p-values...\n")
  create_test_sumstats_multiple_pvals(test_file, n_variants = 10000)
  
  # Run qqplot.R with comma-separated p-value columns
  cat("Running qqplot.R with multiple p-value columns...\n")
  cmd <- sprintf(
    "Rscript ../qqplot.R --file %s --out %s --chrcol CHR --bp_col BP --pval_col P1,P2,P3",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qqplot.R failed with exit code ", result)
  }
  
  # Check that output files exist for all three p-value columns
  for (pcol in c("P1", "P2", "P3")) {
    expected_outputs <- c(
      paste0(output_prefix, "_", pcol, "_qqplot.png"),
      paste0(output_prefix, "_", pcol, "_manhattan.png"),
      paste0(output_prefix, "_", pcol, "_manhattan_loglog.png"),
      paste0(output_prefix, "_", pcol, "_qquantiles.txt")
    )
    
    for (output in expected_outputs) {
      if (!file.exists(output)) {
        stop("Expected output file not created: ", output)
      }
      cat("✓ Created:", output, "\n")
    }
  }
  
  cat("✓ Test passed: All output files created for multiple p-value columns\n")
  
  # Cleanup
  cleanup_files("test_qqplot_multi")
  
  return(TRUE)
}

# Test 3: With minrep_col option
test_qqplot_minrep <- function() {
  cat("\n=== Test 3: With minrep_col option ===\n")
  
  test_file <- "test_qqplot_minrep.tsv"
  output_prefix <- "test_qqplot_minrep_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qqplot_minrep")
  
  # Create test data with minrep format
  cat("Creating test data with minrep format...\n")
  create_test_sumstats_minrep(test_file, n_variants = 10000)
  
  # Add dummy CHR and BP columns (will be overwritten by minrep parsing)
  data <- fread(test_file)
  data[, CHR := "1"]
  data[, BP := 0]
  fwrite(data, test_file, sep = "\t", quote = FALSE)
  
  # Run qqplot.R with minrep_col
  cat("Running qqplot.R with minrep_col...\n")
  cmd <- sprintf(
    "Rscript ../qqplot.R --file %s --out %s --chrcol CHR --bp_col BP --pval_col P --minrep_col minrep_id",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qqplot.R failed with exit code ", result)
  }
  
  # Check that output files exist
  expected_outputs <- c(
    paste0(output_prefix, "_P_qqplot.png"),
    paste0(output_prefix, "_P_manhattan.png"),
    paste0(output_prefix, "_P_manhattan_loglog.png"),
    paste0(output_prefix, "_P_qquantiles.txt")
  )
  
  for (output in expected_outputs) {
    if (!file.exists(output)) {
      stop("Expected output file not created: ", output)
    }
  }
  
  cat("✓ Test passed: minrep_col option worked\n")
  
  # Cleanup
  cleanup_files("test_qqplot_minrep")
  
  return(TRUE)
}

# Test 4: Custom loglog parameters
test_qqplot_custom_loglog <- function() {
  cat("\n=== Test 4: Custom loglog parameters ===\n")
  
  test_file <- "test_qqplot_loglog.tsv"
  output_prefix <- "test_qqplot_loglog_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qqplot_loglog")
  
  # Create test data
  cat("Creating test data...\n")
  create_test_sumstats(test_file, n_variants = 10000)
  
  # Run qqplot.R with custom loglog parameters
  cat("Running qqplot.R with custom loglog parameters...\n")
  cmd <- sprintf(
    "Rscript ../qqplot.R --file %s --out %s --chrcol CHR --bp_col BP --pval_col P --loglog_pval 15 --loglog_ylim 200",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qqplot.R failed with exit code ", result)
  }
  
  # Check that output files exist
  expected_outputs <- c(
    paste0(output_prefix, "_P_qqplot.png"),
    paste0(output_prefix, "_P_manhattan.png"),
    paste0(output_prefix, "_P_manhattan_loglog.png"),
    paste0(output_prefix, "_P_qquantiles.txt")
  )
  
  for (output in expected_outputs) {
    if (!file.exists(output)) {
      stop("Expected output file not created: ", output)
    }
  }
  
  cat("✓ Test passed: Custom loglog parameters worked\n")
  
  # Cleanup
  cleanup_files("test_qqplot_loglog")
  
  return(TRUE)
}

# Test 5: Lambda calculation validation
test_qqplot_lambda <- function() {
  cat("\n=== Test 5: Lambda calculation validation ===\n")
  
  test_file <- "test_qqplot_lambda.tsv"
  output_prefix <- "test_qqplot_lambda_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qqplot_lambda")
  
  # Create test data
  cat("Creating test data...\n")
  create_test_sumstats(test_file, n_variants = 10000)
  
  # Run qqplot.R
  cat("Running qqplot.R...\n")
  cmd <- sprintf(
    "Rscript ../qqplot.R --file %s --out %s --chrcol CHR --bp_col BP --pval_col P",
    test_file, output_prefix
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qqplot.R failed with exit code ", result)
  }
  
  # Read and validate lambda values
  lambda_file <- paste0(output_prefix, "_P_qquantiles.txt")
  if (!file.exists(lambda_file)) {
    stop("Lambda quantiles file not created: ", lambda_file)
  }
  
  lambda_content <- readLines(lambda_file)
  cat("Lambda values:", lambda_content, "\n")
  
  # Basic validation: should have lambda values
  if (length(lambda_content) == 0) {
    stop("Lambda file is empty")
  }
  
  cat("✓ Test passed: Lambda calculation completed\n")
  
  # Cleanup
  cleanup_files("test_qqplot_lambda")
  
  return(TRUE)
}

# Main test runner
main <- function() {
  cat("========================================\n")
  cat("Running qqplot.R test suite\n")
  cat("========================================\n")
  
  # Change to scripts directory if not already there
  script_dir <- dirname(sys.frame(1)$ofile)
  if (script_dir != "") {
    setwd(file.path(script_dir, ".."))
  }
  
  tests <- list(
    test_qqplot_basic,
    test_qqplot_multiple_pvals,
    test_qqplot_minrep,
    test_qqplot_custom_loglog,
    test_qqplot_lambda
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
