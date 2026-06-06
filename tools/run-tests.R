# Run the testthat suite against the installed (already-compiled) package.
# Using env = asNamespace("bqmm") so bare exported AND internal names resolve
# without re-running load_all() (which would recompile the Stan model).
# NOT_CRAN=true so the end-to-end Stan-fit tests actually run.
Sys.setenv(NOT_CRAN = "true")
library(bqmm)            # default library (where R CMD INSTALL now installs)
library(testthat)

pkg_dir <- "C:/Users/kyle/Documents/Projects/Personal/bayesian quantile mixed model/bqmm"
res <- test_dir(file.path(pkg_dir, "tests", "testthat"),
                env = asNamespace("bqmm"),
                reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat("\n=== TEST TOTALS ===\n")
cat("failed:  ", sum(df$failed), "\n")
cat("warnings:", sum(df$warning), "\n")
cat("skipped: ", sum(df$skipped > 0), "\n")
cat("passed:  ", sum(df$passed), "\n")
