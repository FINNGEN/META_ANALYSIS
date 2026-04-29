#!/usr/bin/env Rscript

# Comprehensive tests for qc.R script
# Tests that the script runs successfully with various reasonable inputs
# Note: We don't validate plot quality, just that the script completes and outputs exist

library(data.table)
library(rjson)

# Helper function to create test config JSON
create_test_config <- function(filename, studies = c("study1", "study2", "study3")) {
  config <- list(
    meta = lapply(studies, function(s) {
      list(
        name = s,
        file = paste0(s, ".tsv"),
        n_cases = 1000,
        n_controls = 10000,
        chr = "#CHR",
        pos = "POS",
        ref = "REF",
        alt = "ALT",
        effect = "beta",
        se = "sebeta",
        pval = "pval",
        af_alt = "af_alt"
      )
    })
  )
  
  write(toJSON(config), file = filename)
  return(config)
}

# Helper function to create test meta-analysis summary statistics
create_test_meta_sumstats <- function(filename, studies = c("study1", "study2", "study3"), 
                                      n_variants = 1000, method = "inv_var", 
                                      include_loo = FALSE) {
  set.seed(42)
  
  # Create chromosomes 1-22 and X
  chrs <- c(rep(1:22, each = floor(n_variants / 23)), rep("X", n_variants %% 23))
  
  data <- data.table(
    `#CHR` = chrs,
    POS = sample(1:250000000, n_variants, replace = TRUE),
    REF = sample(c("A", "C", "G", "T"), n_variants, replace = TRUE),
    ALT = sample(c("A", "C", "G", "T"), n_variants, replace = TRUE)
  )
  
  # Add columns for each study
  for (study in studies) {
    data[, paste0(study, "_pval") := runif(n_variants, 0, 1)]
    data[, paste0(study, "_beta") := rnorm(n_variants, 0, 0.1)]
    data[, paste0(study, "_sebeta") := runif(n_variants, 0.01, 0.05)]
    data[, paste0(study, "_af_alt") := runif(n_variants, 0.01, 0.5)]
  }
  
  # Add meta-analysis results
  meta_pval_col <- paste("all", method, "meta_p", sep = "_")
  meta_beta_col <- paste("all", method, "meta_beta", sep = "_")
  meta_sebeta_col <- paste("all", method, "meta_sebeta", sep = "_")
  meta_het_p_col <- paste("all", method, "het_p", sep = "_")
  
  data[, (meta_pval_col) := runif(n_variants, 0, 1)]
  data[, (meta_beta_col) := rnorm(n_variants, 0, 0.08)]
  data[, (meta_sebeta_col) := runif(n_variants, 0.01, 0.04)]
  data[, (meta_het_p_col) := runif(n_variants, 0, 1)]
  
  # Add some significant hits across studies
  sig_indices <- sample(1:n_variants, min(50, n_variants / 20))
  for (study in studies) {
    data[sig_indices, paste0(study, "_pval") := runif(.N, 0, 5e-8)]
  }
  data[sig_indices, (meta_pval_col) := runif(.N, 0, 5e-8)]
  
  # Add leave-one-out results if requested
  if (include_loo && length(studies) >= 3) {
    for (study in studies) {
      loo_pval_col <- paste("leave", study, method, "meta_p", sep = "_")
      loo_beta_col <- paste("leave", study, method, "meta_beta", sep = "_")
      loo_sebeta_col <- paste("leave", study, method, "meta_sebeta", sep = "_")
      
      data[, (loo_pval_col) := runif(n_variants, 0, 1)]
      data[, (loo_beta_col) := rnorm(n_variants, 0, 0.08)]
      data[, (loo_sebeta_col) := runif(n_variants, 0.01, 0.04)]
      
      # Add significant hits
      data[sig_indices, (loo_pval_col) := runif(.N, 0, 5e-8)]
    }
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

# Test 1: Basic functionality with minimal options
test_qc_basic <- function() {
  cat("\n=== Test 1: Basic functionality ===\n")
  
  test_file <- "test_qc_basic.tsv"
  config_file <- "test_qc_basic_conf.json"
  output_prefix <- "test_qc_basic_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qc_basic")
  
  # Create test config and data
  cat("Creating test config and data...\n")
  create_test_config(config_file, studies = c("study1", "study2", "study3"))
  create_test_meta_sumstats(test_file, studies = c("study1", "study2", "study3"), 
                           n_variants = 1000, method = "inv_var")
  
  # Run qc.R
  cat("Running qc.R...\n")
  cmd <- sprintf(
    "Rscript ../qc.R --file %s --out %s --conf %s --method inv_var --pheno test_pheno",
    test_file, output_prefix, config_file
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qc.R failed with exit code ", result)
  }
  
  # Check that output files exist
  expected_outputs <- c(
    paste0(output_prefix, ".5e-08.qc.tsv"),
    paste0(output_prefix, ".5e-06.qc.tsv"),
    paste0(output_prefix, ".5e-08.forest_plots.pdf"),
    paste0(output_prefix, ".5e-06.forest_plots.pdf"),
    paste0(output_prefix, ".5e-08.qc.pdf"),
    paste0(output_prefix, ".5e-06.qc.pdf")
  )
  
  for (output in expected_outputs) {
    if (!file.exists(output)) {
      # Some outputs may not be created if there are too few hits
      cat("⚠ Expected output file not created (may be OK if no significant hits):", output, "\n")
      next
    }
    
    # Check file size > 0
    file_size <- file.info(output)$size
    if (file_size == 0) {
      stop("Output file is empty: ", output)
    }
    cat("✓ Created:", output, "(", file_size, "bytes)\n")
  }
  
  # At minimum, QC TSV files should exist
  qc_tsv <- paste0(output_prefix, ".5e-08.qc.tsv")
  if (!file.exists(qc_tsv)) {
    stop("QC TSV file not created: ", qc_tsv)
  }
  
  # Read and validate QC file has expected structure
  qc_data <- fread(qc_tsv)
  if (!"pheno" %in% colnames(qc_data)) {
    stop("QC file missing 'pheno' column")
  }
  
  cat("✓ Test passed: Basic QC completed\n")
  
  # Cleanup
  cleanup_files("test_qc_basic")
  
  return(TRUE)
}

# Test 2: With leave-one-out option
test_qc_loo <- function() {
  cat("\n=== Test 2: With leave-one-out ===\n")
  
  test_file <- "test_qc_loo.tsv"
  config_file <- "test_qc_loo_conf.json"
  output_prefix <- "test_qc_loo_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qc_loo")
  
  # Create test config and data with LOO
  cat("Creating test config and data with leave-one-out...\n")
  create_test_config(config_file, studies = c("study1", "study2", "study3"))
  create_test_meta_sumstats(test_file, studies = c("study1", "study2", "study3"), 
                           n_variants = 1000, method = "inv_var", include_loo = TRUE)
  
  # Run qc.R with --loo flag
  cat("Running qc.R with leave-one-out...\n")
  cmd <- sprintf(
    "Rscript ../qc.R --file %s --out %s --conf %s --method inv_var --pheno test_pheno --loo",
    test_file, output_prefix, config_file
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qc.R failed with exit code ", result)
  }
  
  # Check that QC file exists
  qc_tsv <- paste0(output_prefix, ".5e-08.qc.tsv")
  if (!file.exists(qc_tsv)) {
    stop("QC TSV file not created: ", qc_tsv)
  }
  
  cat("✓ Test passed: Leave-one-out QC completed\n")
  
  # Cleanup
  cleanup_files("test_qc_loo")
  
  return(TRUE)
}

# Test 3: With weighted option
test_qc_weighted <- function() {
  cat("\n=== Test 3: With weighted regression ===\n")
  
  test_file <- "test_qc_weighted.tsv"
  config_file <- "test_qc_weighted_conf.json"
  output_prefix <- "test_qc_weighted_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qc_weighted")
  
  # Create test config and data
  cat("Creating test config and data...\n")
  create_test_config(config_file, studies = c("study1", "study2"))
  create_test_meta_sumstats(test_file, studies = c("study1", "study2"), 
                           n_variants = 1000, method = "inv_var")
  
  # Run qc.R with --weighted flag
  cat("Running qc.R with weighted regression...\n")
  cmd <- sprintf(
    "Rscript ../qc.R --file %s --out %s --conf %s --method inv_var --pheno test_pheno --weighted",
    test_file, output_prefix, config_file
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qc.R failed with exit code ", result)
  }
  
  # Check that QC file exists
  qc_tsv <- paste0(output_prefix, ".5e-08.qc.tsv")
  if (!file.exists(qc_tsv)) {
    stop("QC TSV file not created: ", qc_tsv)
  }
  
  cat("✓ Test passed: Weighted regression QC completed\n")
  
  # Cleanup
  cleanup_files("test_qc_weighted")
  
  return(TRUE)
}

# Test 4: Custom p-value thresholds
test_qc_custom_pval_thresh <- function() {
  cat("\n=== Test 4: Custom p-value thresholds ===\n")
  
  test_file <- "test_qc_pval.tsv"
  config_file <- "test_qc_pval_conf.json"
  output_prefix <- "test_qc_pval_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qc_pval")
  
  # Create test config and data
  cat("Creating test config and data...\n")
  create_test_config(config_file, studies = c("study1", "study2"))
  create_test_meta_sumstats(test_file, studies = c("study1", "study2"), 
                           n_variants = 1000, method = "inv_var")
  
  # Run qc.R with custom p-value thresholds
  cat("Running qc.R with custom p-value thresholds...\n")
  cmd <- sprintf(
    "Rscript ../qc.R --file %s --out %s --conf %s --method inv_var --pheno test_pheno --pval_thresh '1e-5, 1e-4'",
    test_file, output_prefix, config_file
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qc.R failed with exit code ", result)
  }
  
  # Check that QC files exist for custom thresholds
  expected_qc_files <- c(
    paste0(output_prefix, ".1e-05.qc.tsv"),
    paste0(output_prefix, ".0.0001.qc.tsv")
  )
  
  for (qc_file in expected_qc_files) {
    if (!file.exists(qc_file)) {
      stop("QC file not created for custom threshold: ", qc_file)
    }
    cat("✓ Created:", qc_file, "\n")
  }
  
  cat("✓ Test passed: Custom p-value thresholds worked\n")
  
  # Cleanup
  cleanup_files("test_qc_pval")
  
  return(TRUE)
}

# Test 5: Keep HLA region
test_qc_keep_hla <- function() {
  cat("\n=== Test 5: Keep HLA region variants ===\n")
  
  test_file <- "test_qc_hla.tsv"
  config_file <- "test_qc_hla_conf.json"
  output_prefix <- "test_qc_hla_output"
  
  # Cleanup any existing test files
  cleanup_files("test_qc_hla")
  
  # Create test config and data
  cat("Creating test config and data with HLA variants...\n")
  create_test_config(config_file, studies = c("study1", "study2"))
  
  # Create data with some HLA region variants (chr 6, 20-40Mb)
  set.seed(42)
  n_variants <- 1000
  studies <- c("study1", "study2")
  method <- "inv_var"
  
  data <- data.table(
    `#CHR` = c(rep(6, n_variants / 2), rep(1:22, each = floor(n_variants / 2 / 22))),
    POS = c(
      sample(20e6:40e6, n_variants / 2, replace = TRUE),  # HLA region
      sample(1:250000000, n_variants / 2, replace = TRUE)  # Other positions
    ),
    REF = sample(c("A", "C", "G", "T"), n_variants, replace = TRUE),
    ALT = sample(c("A", "C", "G", "T"), n_variants, replace = TRUE)
  )
  
  for (study in studies) {
    data[, paste0(study, "_pval") := runif(n_variants, 0, 1)]
    data[, paste0(study, "_beta") := rnorm(n_variants, 0, 0.1)]
    data[, paste0(study, "_sebeta") := runif(n_variants, 0.01, 0.05)]
    data[, paste0(study, "_af_alt") := runif(n_variants, 0.01, 0.5)]
  }
  
  meta_pval_col <- paste("all", method, "meta_p", sep = "_")
  meta_beta_col <- paste("all", method, "meta_beta", sep = "_")
  meta_sebeta_col <- paste("all", method, "meta_sebeta", sep = "_")
  meta_het_p_col <- paste("all", method, "het_p", sep = "_")
  
  data[, (meta_pval_col) := runif(n_variants, 0, 1)]
  data[, (meta_beta_col) := rnorm(n_variants, 0, 0.08)]
  data[, (meta_sebeta_col) := runif(n_variants, 0.01, 0.04)]
  data[, (meta_het_p_col) := runif(n_variants, 0, 1)]
  
  sig_indices <- sample(1:n_variants, min(50, n_variants / 20))
  for (study in studies) {
    data[sig_indices, paste0(study, "_pval") := runif(.N, 0, 5e-8)]
  }
  data[sig_indices, (meta_pval_col) := runif(.N, 0, 5e-8)]
  
  fwrite(data, test_file, sep = "\t", quote = FALSE)
  
  # Run qc.R with --keep_hla flag
  cat("Running qc.R with --keep_hla...\n")
  cmd <- sprintf(
    "Rscript ../qc.R --file %s --out %s --conf %s --method inv_var --pheno test_pheno --keep_hla",
    test_file, output_prefix, config_file
  )
  
  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
  
  # Check that command succeeded
  if (result != 0) {
    stop("qc.R failed with exit code ", result)
  }
  
  # Check that QC file exists
  qc_tsv <- paste0(output_prefix, ".5e-08.qc.tsv")
  if (!file.exists(qc_tsv)) {
    stop("QC TSV file not created: ", qc_tsv)
  }
  
  cat("✓ Test passed: Keep HLA option worked\n")
  
  # Cleanup
  cleanup_files("test_qc_hla")
  
  return(TRUE)
}

# Main test runner
main <- function() {
  cat("========================================\n")
  cat("Running qc.R test suite\n")
  cat("========================================\n")
  
  # Change to scripts directory if not already there
  script_dir <- dirname(sys.frame(1)$ofile)
  if (script_dir != "") {
    setwd(file.path(script_dir, ".."))
  }
  
  tests <- list(
    test_qc_basic,
    test_qc_loo,
    test_qc_weighted,
    test_qc_custom_pval_thresh,
    test_qc_keep_hla
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
