Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(bqmm, lib.loc = file.path(Sys.getenv("TEMP"), "bqmm_lib")))
library(testthat)
pkg <- "C:/Users/kyle/Documents/Projects/Personal/bayesian quantile mixed model/bqmm"
res <- test_file(file.path(pkg, "tests/testthat/test-corr.R"),
                 env = asNamespace("bqmm"), reporter = "summary")
df <- as.data.frame(res)
cat("\nfailed:", sum(df$failed), " passed:", sum(df$passed),
    " skipped:", sum(df$skipped > 0), "\n")
