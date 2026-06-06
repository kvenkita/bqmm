# Stand-alone smoke test of the Stan-free logic (no package build / no Stan).
# Sources the relevant R files and exercises formula parsing, the YWH sandwich,
# rearrangement and the ALD family. Run with:
#   Rscript tools/smoke-test.R
suppressMessages({
  library(lme4)
  library(Matrix)
})

src <- function(f) source(file.path("R", f), local = .GlobalEnv)
invisible(lapply(c("family.R", "priors.R", "formula.R", "standata.R",
                   "adjust.R", "noncrossing.R"), src))

ok <- function(cond, msg) {
  if (isTRUE(cond)) cat("PASS:", msg, "\n")
  else stop("FAIL: ", msg, call. = FALSE)
}
approx <- function(a, b, tol = 1e-5) all(abs(a - b) < tol)

## --- formula parsing -------------------------------------------------------
set.seed(1)
n <- 60
d <- data.frame(y = rnorm(n), x = rnorm(n),
                g = factor(rep(1:6, length.out = n)),
                h = factor(rep(1:3, each = n / 3)))

p1 <- bqmm_parse_formula(y ~ x + (1 | g), data = d)
ok(ncol(p1$X) == 2, "random intercept: X has 2 columns")
ok(ncol(p1$Z) == nlevels(d$g), "random intercept: Z has one col per level")
ok(length(p1$re_components) == 1, "random intercept: one variance component")
ok(all(p1$sd_map == 1), "random intercept: all Z cols share the component")

p2 <- bqmm_parse_formula(y ~ x + (1 + x | g), data = d)
ok(ncol(p2$Z) == 2 * nlevels(d$g), "random slope: Z has 2*levels columns")
ok(length(p2$re_components) == 2, "random slope: two variance components")
ok(all(table(p2$sd_map) == nlevels(d$g)), "random slope: even split of columns")

p3 <- bqmm_parse_formula(y ~ x + (1 | g) + (1 | h), data = d)
ok(ncol(p3$Z) == nlevels(d$g) + nlevels(d$h), "crossed: Z column count")
ok(length(p3$re_components) == 2, "crossed: two variance components")

Zref <- as.matrix(Matrix::t(
  lme4::glFormula(y ~ x + (1 + x | g), data = d, family = gaussian())$reTrms$Zt))
ok(approx(unname(p2$Z), unname(Zref)), "Z matches t(Zt) from lme4")

## --- standata --------------------------------------------------------------
pr <- bqmm_default_priors(NULL, p1$y, ncol(p1$X))
sd <- bqmm_standata(p1, tau = 0.5, prior = pr)
ok(sd$N == n && sd$K == 2 && sd$Q == ncol(p1$Z), "standata dimensions")
ok(length(sd$sd_map) == sd$Q, "standata sd_map length")
ok(sd$p == 0.5, "standata carries tau")

## --- ALD family ------------------------------------------------------------
for (tau in c(0.25, 0.5, 0.8)) {
  area <- integrate(function(x) dald(x, 0, 1, tau), -Inf, Inf)$value
  ok(approx(area, 1, 1e-4), sprintf("ALD density integrates to 1 (tau=%.2f)", tau))
  cdf0 <- integrate(function(x) dald(x, 0, 1, tau), -Inf, 0)$value
  ok(approx(cdf0, tau, 1e-4), sprintf("ALD location is the tau-quantile (tau=%.2f)", tau))
}
set.seed(42)
yq <- rald(2e5, 0, 1, 0.75)
ok(approx(mean(yq <= 0), 0.75, 0.01), "rald empirical quantile ~ tau")

## --- Yang-Wang-He sandwich -------------------------------------------------
set.seed(2)
N <- 500; tau <- 0.5
X <- matrix(1, N, 1)
resid <- rald(N, 0, 1, tau)
sw <- compute_ywh_sandwich(X, resid, tau = tau)
h  <- hall_sheather_bandwidth(N, tau)
f0 <- mean(abs(resid) < h) / (2 * h)
ok(approx(as.numeric(sw$vcov), tau * (1 - tau) / (N * f0^2), 1e-6),
   "YWH one-coef closed form matches analytic")

Xk <- cbind(1, rnorm(300))
rk <- rald(300, 0, 1, 0.7)
swk <- compute_ywh_sandwich(Xk, rk, tau = 0.7)
ok(isSymmetric(unname(swk$vcov)), "YWH vcov symmetric")
ok(all(eigen(swk$vcov, only.values = TRUE)$values >= -1e-8), "YWH vcov PSD")

g <- rep(1:12, each = 20)
Xc <- cbind(1, rnorm(240))
base <- rald(12, 0, 1, 0.5)[g]
rc <- base + rald(240, 0, 0.3, 0.5)
indep <- compute_ywh_sandwich(Xc, rc, 0.5)$vcov
clust <- compute_ywh_sandwich(Xc, rc, 0.5, groups = g)$vcov
ok(!isTRUE(all.equal(indep, clust)), "cluster-robust meat differs from independence")

## --- noncrossing -----------------------------------------------------------
m <- rbind(c(1, 0.5, 2), c(0, 1, 0.8), c(-1, -2, 0))
ok(is_noncrossing(rearrange_quantiles(m)), "rearrangement removes crossings")
mono <- rbind(c(0.1, 0.5, 0.9), c(-1, 0, 1))
ok(approx(unname(rearrange_quantiles(mono)), unname(mono)),
   "rearrangement is identity when monotone")

cat("\nAll smoke-test assertions passed.\n")
