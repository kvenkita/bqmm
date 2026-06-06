## Adversarial verification of ALD density and RNG
source("R/family.R")

set.seed(123)

cat("==== TASK 1: density self-consistency & quantile ====\n")

## 1a. integrate to 1
for (tau in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  for (sig in c(0.5, 2.5)) {
    I <- integrate(function(x) dald(x, mu = 1.0, sigma = sig, tau = tau),
                   lower = -Inf, upper = Inf)$value
    cat(sprintf("integral tau=%.1f sigma=%.1f : %.8f\n", tau, sig, I))
  }
}

## 1b. P(Y<=mu) == tau  (CDF at mu)
cat("\n-- P(Y<=mu) via numeric integration of dald --\n")
for (tau in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  P <- integrate(function(x) dald(x, mu = 1.0, sigma = 2.5, tau = tau),
                 lower = -Inf, upper = 1.0)$value
  cat(sprintf("tau=%.1f : P(Y<=mu)=%.6f  (should be %.1f)\n", tau, P, tau))
}

## 1c. R dald vs explicit formula tau(1-tau)/sigma exp(-rho)
cat("\n-- dald vs explicit formula --\n")
xs <- c(-3, -1, 0.5, 1.0, 2.3, 5)
mu <- 1.0; sig <- 2.5; tau <- 0.3
rho <- function(u, t) u * (t - (u < 0))
explicit <- tau*(1-tau)/sig * exp(-rho((xs-mu)/sig, tau))
fromfun  <- dald(xs, mu, sig, tau)
cat("max abs diff R dald vs explicit:", max(abs(explicit - fromfun)), "\n")

## 1d. replicate Stan lpdf scalar and compare to log(dald)
stan_lpdf_scalar <- function(y, mu, sigma, p) {
  u <- (y - mu)/sigma
  log(p) + log1p(-p) - log(sigma) - rho(u, p)
}
diff_stan <- max(abs(stan_lpdf_scalar(xs, mu, sig, tau) - dald(xs, mu, sig, tau, log = TRUE)))
cat("max abs diff Stan-scalar-lpdf vs R dald(log):", diff_stan, "\n")

cat("\n==== TASK 2: RNG (Kozumi-Kobayashi) ====\n")
tau <- 0.3; sig <- 2.5; mu <- 1.0
N <- 1e6
y <- rald(N, mu = mu, sigma = sig, tau = tau)

## (a) empirical P(Y<=mu)
cat(sprintf("(a) empirical P(Y<=mu) = %.5f  (target tau=%.2f)\n", mean(y <= mu), tau))

## (b) variance: closed form
var_cf <- sig^2 * (1 - 2*tau + 2*tau^2) / ((1-tau)^2 * tau^2)
cat(sprintf("(b) empirical var = %.5f ; closed-form ALD var = %.5f ; ratio = %.5f\n",
            var(y), var_cf, var(y)/var_cf))

## also check mean closed form: E[Y] = mu + sigma*(1-2tau)/(tau(1-tau))
mean_cf <- mu + sig*(1-2*tau)/(tau*(1-tau))
cat(sprintf("    empirical mean = %.5f ; closed-form mean = %.5f\n", mean(y), mean_cf))

## (c) KS test against dald numerical CDF
pald <- function(q, mu, sigma, tau) {
  sapply(q, function(qq) integrate(function(x) dald(x, mu, sigma, tau),
                                   lower = -Inf, upper = qq,
                                   rel.tol = 1e-9)$value)
}
samp <- sample(y, 5000)   # KS on a subsample to keep CDF evals tractable
ks <- ks.test(samp, function(q) pald(q, mu, sig, tau))
cat(sprintf("(c) KS D = %.5f, p-value = %.4f\n", ks$statistic, ks$p.value))

## Extra: sweep several tau to make sure P(Y<=mu)=tau holds generally
cat("\n-- sweep empirical P(Y<=mu) over tau (sigma=2.5, N=1e6) --\n")
for (t in c(0.1, 0.5, 0.7, 0.9)) {
  yy <- rald(N, mu = mu, sigma = sig, tau = t)
  vcf <- sig^2 * (1 - 2*t + 2*t^2) / ((1-t)^2 * t^2)
  cat(sprintf("tau=%.1f : P(Y<=mu)=%.5f  var ratio(emp/cf)=%.5f\n",
              t, mean(yy <= mu), var(yy)/vcf))
}

cat("\n==== TASK: sigma-scaling check ====\n")
## Does variance scale as sigma^2? Compare sigma=1 vs sigma=2.5 empirically.
y1 <- rald(N, mu=0, sigma=1,   tau=0.3)
y2 <- rald(N, mu=0, sigma=2.5, tau=0.3)
cat(sprintf("var(sigma=1)=%.4f var(sigma=2.5)=%.4f ratio=%.4f (sigma^2 ratio=%.4f)\n",
            var(y1), var(y2), var(y2)/var(y1), 2.5^2))
