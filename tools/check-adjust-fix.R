# Quick standalone check of the revised YWH sandwich (no package build).
src <- function(f) source(file.path("R", f), local = .GlobalEnv)
invisible(lapply(c("family.R", "adjust.R"), src))
ok <- function(c, m) if (isTRUE(c)) cat("PASS:", m, "\n") else stop("FAIL: ", m)

set.seed(11)
n <- 400
X <- cbind(1, rnorm(n))
resid <- rald(n, 0, 1, 0.5)

# tiny bandwidth would have triggered the old broken fallback
sw_tiny <- compute_ywh_sandwich(X, resid, tau = 0.5, bandwidth = 1e-8)
sw_ok   <- compute_ywh_sandwich(X, resid, tau = 0.5)
ok(all(is.finite(sw_tiny$vcov)), "tiny-bandwidth vcov finite")
ok(qr(sw_tiny$vcov)$rank == 2, "tiny-bandwidth vcov full rank")
ratio <- diag(sw_tiny$vcov) / diag(sw_ok$vcov)
cat("  tiny/ok diag ratio:", round(ratio, 3), "\n")
ok(all(ratio > 0.2 & ratio < 5), "tiny-bandwidth same order of magnitude (no collapse)")
cat("  bandwidth grew from 1e-8 to", signif(sw_tiny$bandwidth, 3), "\n")

# 1/n scaling, intercept-only (closed form V = tau(1-tau)/(n f0^2), stable)
v <- function(n) as.numeric(compute_ywh_sandwich(matrix(1, n, 1),
                                                 rald(n, 0, 1, 0.5), 0.5)$vcov)
v1 <- mean(sapply(1:60, function(i) v(500)))
v8 <- mean(sapply(1:60, function(i) v(4000)))
cat("  v(500)/v(4000) ratio:", round(v1/v8, 2), "(expect ~8 for 1/n; 64 for 1/n^2)\n")
ok(v1/v8 > 5.5 & v1/v8 < 11, "variance scales like 1/n")

# analytic one-coef closed form still matches
N <- 4000; X1 <- matrix(1, N, 1); r <- rald(N, 0, 1, 0.5)
sw <- compute_ywh_sandwich(X1, r, 0.5)
h <- hall_sheather_bandwidth(N, 0.5); f0 <- mean(abs(r) < h)/(2*h)
ok(abs(as.numeric(sw$vcov) - 0.25/(N*f0^2)) < 1e-6, "one-coef closed form matches analytic")

# vs quantreg on identical residuals (sandwich agreement, averaged)
if (requireNamespace("quantreg", quietly = TRUE)) {
  set.seed(5); reps <- 30; rr <- matrix(NA, reps, 2); rq_se <- matrix(NA, reps, 2)
  for (i in 1:reps) {
    nn <- 800; xx <- rnorm(nn); y <- 1 + 2*xx + rald(nn, 0, 1, 0.5)
    f <- quantreg::rq(y ~ xx, tau = 0.5)
    rr[i,] <- diag(compute_ywh_sandwich(cbind(1, xx), residuals(f), 0.5)$vcov)
    s <- summary(f, se = "nid", cov = TRUE)
    rq_se[i,] <- diag(s$cov)
  }
  cat("  mean bqmm sandwich diag:", round(colMeans(rr), 5), "\n")
  cat("  mean quantreg nid  diag:", round(colMeans(rq_se), 5), "\n")
  ok(all(colMeans(rr)/colMeans(rq_se) > 0.7 & colMeans(rr)/colMeans(rq_se) < 1.4),
     "sandwich agrees with quantreg nid (averaged)")
}
cat("\nadjust.R fix checks passed.\n")
