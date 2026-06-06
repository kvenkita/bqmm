# Full CRAN check, running \donttest examples (as CRAN does).
res <- rcmdcheck::rcmdcheck(
  path = ".",
  args = c("--as-cran", "--no-manual", "--run-donttest"),
  build_args = "--no-manual",
  error_on = "never",
  check_dir = file.path(Sys.getenv("TEMP"), "bqmm_check_cran")
)
cat("\n================ CHECK SUMMARY ================\n")
cat("errors:  ", length(res$errors), "\n")
cat("warnings:", length(res$warnings), "\n")
cat("notes:   ", length(res$notes), "\n\n")
if (length(res$errors))   { cat("---- ERRORS ----\n");   cat(res$errors,   sep = "\n\n") }
if (length(res$warnings)) { cat("---- WARNINGS ----\n"); cat(res$warnings, sep = "\n\n") }
if (length(res$notes))    { cat("---- NOTES ----\n");    cat(res$notes,    sep = "\n\n") }
