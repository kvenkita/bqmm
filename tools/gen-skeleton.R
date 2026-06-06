# Generate a throwaway rstantools skeleton at a FIXED path so we can copy its
# glue files. Wrapped in try() because rstan_create_package tries to open a
# Read-and-delete-me file in RStudio at the end (harmless under Rscript).
options(usethis.quiet = TRUE)
root <- file.path(Sys.getenv("TEMP"), "bqmmskel")
unlink(root, recursive = TRUE)

try(
  rstantools::rstan_create_package(
    path    = root,
    roxygen = FALSE,
    travis  = FALSE,
    license = FALSE
  ),
  silent = TRUE
)

cat("ROOT:", root, "\n")
cat("EXISTS:", dir.exists(root), "\n")
files <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
cat(paste0("  ", files), sep = "\n")
cat("\n")
