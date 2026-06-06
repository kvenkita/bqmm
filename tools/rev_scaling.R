source("R/family.R")
source("R/adjust.R")
suppressMessages(library(quantreg))

cat("=== TASK 1: n-scaling test ===\n")
cat("If V = D1^{-1} D0 D1^{-1}/n with D0,D1 already /n, does Var ~ 1/n?\n\n")

# Monte Carlo: empirical Var(betahat) across reps vs our sandwich, across n.
tau <- 0.5
f0_true <- dnorm(0)
for (n in c(500, 1000, 2000, 4000)) {
  reps <- 400
  bhats <- matrix(NA, reps, 3)
  sw_diag <- matrix(NA, reps, 3)
  for (r in 1:reps) {
    x1 <- rnorm(n); x2 <- rnorm(n)
    X <- cbind(1, x1, x2)
    y <- as.numeric(X %*% c(1,2,-1)) + rnorm(n)
    fit <- rq.fit.br(X, y, tau = tau)
    bhats[r,] <- fit$coefficients
    sw <- compute_ywh_sandwich(X, as.numeric(fit$residuals), tau = tau)
    sw_diag[r,] <- diag(sw$vcov)
  }
  emp_var <- apply(bhats, 2, var)
  mean_sw <- colMeans(sw_diag)
  cat(sprintf("n=%5d  empVar=[%.2e %.2e %.2e]  meanSand=[%.2e %.2e %.2e]  ratio=[%.3f %.3f %.3f]\n",
      n, emp_var[1], emp_var[2], emp_var[3],
      mean_sw[1], mean_sw[2], mean_sw[3],
      mean_sw[1]/emp_var[1], mean_sw[2]/emp_var[2], mean_sw[3]/emp_var[3]))
}
cat("\nratio ~ 1.0 across all n  => n-scaling correct.\n")
cat("ratio constant but != 1 => density-bias only (kernel), scaling still fine.\n")
