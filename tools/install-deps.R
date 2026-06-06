options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("rstan", "rstantools", "StanHeaders", "posterior", "RcppParallel",
          "BH", "RcppEigen", "bayesplot", "lqmm", "nlme", "testthat",
          "knitr", "rmarkdown", "devtools", "roxygen2")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
cat("Installing:", paste(need, collapse = ", "), "\n")
if (length(need)) install.packages(need, type = "binary")
cat("DONE pkg install\n")
for (p in pkgs)
  cat(sprintf("  %-13s %s\n", p,
              ifelse(requireNamespace(p, quietly = TRUE),
                     as.character(packageVersion(p)), "MISSING")))
