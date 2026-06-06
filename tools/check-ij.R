# Standalone check of compute_ij (pure numeric, no package build).
source("R/ij.R")
ok <- function(c, m) if (isTRUE(c)) cat("PASS:", m, "\n") else stop("FAIL: ", m)

set.seed(1)
S <- 2000; K <- 2; n <- 240
beta <- matrix(rnorm(S * K), S, K)
# make per-obs loglik correlated with beta so influences are non-trivial
a <- matrix(rnorm(n * K, sd = 0.3), K, n)
loglik <- beta %*% a + matrix(rnorm(S * n, sd = 0.5), S, n)
groups <- rep(1:24, each = 10)

Vind <- compute_ij(beta, loglik)
Vcl  <- compute_ij(beta, loglik, groups = groups)

ok(all(dim(Vind) == c(K, K)), "indep IJ is K x K")
ok(isSymmetric(Vind) && isSymmetric(Vcl), "IJ matrices symmetric")
ok(all(eigen(Vind, only.values = TRUE)$values > -1e-9), "indep IJ PSD")
ok(all(eigen(Vcl, only.values = TRUE)$values > -1e-9), "cluster IJ PSD")
ok(!isTRUE(all.equal(Vind, Vcl)), "cluster IJ differs from independence IJ")
cat("indep IJ diag:", round(diag(Vind), 5), "\n")
cat("cluster IJ diag:", round(diag(Vcl), 5), "\n")
cat("\ncompute_ij OK\n")
