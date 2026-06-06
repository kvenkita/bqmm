# Generate man/ pages and NAMESPACE from roxygen tags, loading the already-
# installed package from the default library (no recompile). Then print the
# generated NAMESPACE.
roxygen2::roxygenise(package.dir = ".")
cat("\n===== NAMESPACE =====\n")
cat(readLines("NAMESPACE"), sep = "\n")
cat("\n===== man/ =====\n")
cat(list.files("man"), sep = "\n")
