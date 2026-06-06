## Adversarial review checks for the YWH sandwich
source("R/family.R")
source("R/adjust.R")
suppressMessages(library(quantreg))

cat("==================================================================\n")
cat("TASK 1 & 2: compare compute_ywh_sandwich to quantreg summary.rq\n")
cat("==================================================================\n")

set.seed(42)
n <- 2000
x1 <- rnorm(n)
x2 <- rnorm(n)
X  <- cbind(1, x1, x2)
beta_true <- c(1, 2, -1)
# homoscedastic errors so f_i(0) is constant across i
err <- rnorm(n)                       # standard normal errors
tau <- 0.5
y <- as.numeric(X %*% beta_true) + err

fit <- rq(y ~ x1 + x2, tau = tau)
resid_rq <- residuals(fit)

# quantreg "ker" (Powell kernel) sandwich
s_ker <- summary(fit, se = "ker", covariance = TRUE)
V_ker <- s_ker$cov
# quantreg "nid" sandwich
s_nid <- summary(fit, se = "nid", covariance = TRUE)
V_nid <- s_nid$cov

# Our sandwich on the SAME residuals & X (using rq residuals)
sw <- compute_ywh_sandwich(X = X, resid = resid_rq, tau = tau)
V_ours <- sw$vcov

cat("\nbandwidth (hall_sheather, ours):", sw$bandwidth, "\n")
# quantreg's own bandwidth for comparison:
hs_quantreg <- quantreg::bandwidth.rq(tau, n, hs = TRUE)
cat("bandwidth.rq (hs=TRUE):", hs_quantreg, "\n\n")

cat("Our V diag:    ", diag(V_ours), "\n")
cat("quantreg ker:  ", diag(V_ker), "\n")
cat("quantreg nid:  ", diag(V_nid), "\n\n")

cat("Ratio ours/ker:", diag(V_ours)/diag(V_ker), "\n")
cat("Ratio ours/nid:", diag(V_ours)/diag(V_nid), "\n\n")

# Analytic truth for homoscedastic N(0,1) at tau=0.5:
# f(0) = dnorm(0) = 0.3989; Var = tau(1-tau)/f0^2 * (X'X)^{-1}
f0_true <- dnorm(0)
XtX_inv <- solve(crossprod(X))
V_analytic <- tau*(1-tau)/f0_true^2 * XtX_inv
cat("Analytic V diag (true f0):", diag(V_analytic), "\n")
cat("Ratio ours/analytic:", diag(V_ours)/diag(V_analytic), "\n")
cat("Ratio ker/analytic: ", diag(V_ker)/diag(V_analytic), "\n\n")

# What density does our Powell estimate give vs truth?
h <- sw$bandwidth
f0_powell <- mean(abs(resid_rq) < h) / (2*h)
cat("f0 true:", f0_true, " f0 powell(ours):", f0_powell, "\n")
