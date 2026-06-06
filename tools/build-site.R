# Build the pkgdown site into docs/ without reinstalling (uses the installed
# package). Run: Rscript tools/build-site.R
pkgdown::build_site(pkg = ".", install = FALSE, new_process = FALSE,
                    preview = FALSE)
cat("SITE BUILD DONE\n")
