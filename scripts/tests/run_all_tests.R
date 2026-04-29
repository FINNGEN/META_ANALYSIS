#!/usr/bin/env Rscript

# Master test runner for all R scripts
# This script runs all test suites and provides a summary

cat("================================================================================\n")
cat("MASTER TEST RUNNER FOR R SCRIPTS\n")
cat("================================================================================\n\n")

# Change to scripts directory
script_dir <- dirname(sys.frame(1)$ofile)
if (script_dir != "") {
  setwd(script_dir)
}

# Define test scripts
test_scripts <- c(
  "test_miami.R",
  "test_qqplot.R", 
  "test_qc.R"
)

# Track results
results <- data.frame(
  script = character(),
  status = character(),
  stringsAsFactors = FALSE
)

# Run each test suite
for (test_script in test_scripts) {
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("Running:", test_script, "\n")
  cat(strrep("=", 80), "\n")
  
  test_path <- test_script
  
  if (!file.exists(test_path)) {
    cat("✗ Test file not found:", test_path, "\n")
    results <- rbind(results, data.frame(script = test_script, status = "NOT FOUND"))
    next
  }
  
  # Run the test script
  exit_code <- system(paste("Rscript", test_path), intern = FALSE, ignore.stderr = FALSE)
  
  if (exit_code == 0) {
    cat("\n✓ Test suite passed:", test_script, "\n")
    results <- rbind(results, data.frame(script = test_script, status = "PASSED"))
  } else {
    cat("\n✗ Test suite failed:", test_script, "(exit code:", exit_code, ")\n")
    results <- rbind(results, data.frame(script = test_script, status = "FAILED"))
  }
}

# Print summary
cat("\n")
cat(strrep("=", 80), "\n")
cat("TEST SUMMARY\n")
cat(strrep("=", 80), "\n\n")

print(results, row.names = FALSE)

passed <- sum(results$status == "PASSED")
failed <- sum(results$status == "FAILED")
not_found <- sum(results$status == "NOT FOUND")

cat("\n")
cat(sprintf("Total test suites: %d\n", nrow(results)))
cat(sprintf("Passed: %d\n", passed))
cat(sprintf("Failed: %d\n", failed))
cat(sprintf("Not found: %d\n", not_found))
cat("\n")

if (failed > 0 || not_found > 0) {
  cat("⚠ Some tests failed or were not found!\n")
  quit(status = 1)
} else {
  cat("✓ All test suites passed!\n")
  quit(status = 0)
}
