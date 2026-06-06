source("R/family.R"); source("R/adjust.R")
set.seed(1)
cat("Powell uniform-kernel density bias as n grows (true f0 = dnorm(0) = 0.3989):\n\n")
for (n in c(2000, 20000, 200000, 1000000)) {
  resid <- rnorm(n, 0, 1); tau <- 0.5
  h <- hall_sheather_bandwidth(n, tau)
  f0_pow <- mean(abs(resid) < h) / (2*h)
  cat(sprintf("n=%8d  h=%.5f  f0_powell=%.4f  ratio_to_true=%.4f\n",
      n, h, f0_pow, f0_pow/dnorm(0)))
}
cat("\nIf ratio -> 1 as n grows, estimator is consistent (just finite-sample biased).\n")
cat("If ratio plateaus below 1, there is a normalization bug.\n")
