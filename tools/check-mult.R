source("R/family.R"); source("R/adjust.R")
set.seed(1)
X <- cbind(1, rnorm(120)); r <- rald(120, 0, 1, 0.5); Sp <- diag(c(.02, .02))
g <- rep(1:12, each = 10)
ad <- compute_ywh_multiplicative(Sp, X, r, sigma = 1, tau = 0.5, groups = g)
stopifnot(all(dim(ad$vcov) == c(2, 2)),
          isSymmetric(ad$vcov),
          all(eigen(ad$vcov, only.values = TRUE)$values > -1e-9),
          isTRUE(all.equal(unname(ad$vcov), unname(Sp %*% ad$G %*% Sp))))
# sigma scaling
v1 <- compute_ywh_multiplicative(Sp, X, r, 1, 0.5)$vcov
v2 <- compute_ywh_multiplicative(Sp, X, r, 2, 0.5)$vcov
stopifnot(isTRUE(all.equal(as.numeric(v1 / v2), rep(4, 4))))
cat("compute_ywh_multiplicative OK\n")
