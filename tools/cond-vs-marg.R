# Conditional vs marginal quantile estimands in a mixed model.
#
# DGP (two-level, heteroscedastic location-scale):
#   y_ij = b0 + b1 x_ij + u_j + (a + c*x_ij) * eps_ij,
#   x ~ U(0.5, 2.5),  u_j ~ N(0, su^2),  eps ~ N(0,1).
#
# Conditional tau-quantile given the random effect u_j (what bqmm's fixed
# effects target) is LINEAR in x:
#   Q_tau(y | x, u_j) = [b0 + u_j + a*z_tau] + (b1 + c*z_tau) * x,  z_tau=qnorm(tau)
#   => conditional fixed slope = b1 + c*z_tau   (closed form).
#
# Marginal tau-quantile (over u and eps; what a marginal QR + cluster-robust SE,
# e.g. Ji et al. 2024, targets) is NONLINEAR in x:
#   y | x ~ N(b0 + b1 x, su^2 + (a + c x)^2)
#   => Q_tau(y | x) = b0 + b1 x + z_tau * sqrt(su^2 + (a + c x)^2).
# A linear marginal QR estimates the best-linear projection of this curve; its
# slope differs from the conditional slope.
#
# This script (1) tabulates conditional vs (oracle) marginal slopes across tau,
# and (2) checks that bqmm recovers the CONDITIONAL slope with ~nominal coverage,
# while its estimate does NOT match the marginal slope -- i.e. the two methods
# answer different questions.

suppressMessages({ library(bqmm); library(quantreg) })
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else 60L
outdir <- file.path("tools", "validation"); dir.create(outdir, FALSE, TRUE)
logmsg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n"); flush.console() }

b0 <- 1; b1 <- 2; a <- 0.6; c <- 0.9; su <- 0.8
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
cond_slope <- function(tau) b1 + c * qnorm(tau)
cond_int   <- function(tau) b0 + a * qnorm(tau)   # fixed conditional intercept

sim <- function(G, npg, seed) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), each = npg)); x <- runif(G * npg, 0.5, 2.5)
  u <- rnorm(G, 0, su)[g]
  y <- b0 + b1 * x + u + (a + c * x) * rnorm(G * npg)
  data.frame(y = y, x = x, g = g)
}

## (1) Oracle marginal slope (best linear projection) via a huge rq fit -------
logmsg("computing oracle marginal slopes ...")
big <- sim(G = 4000, npg = 30, seed = 1)   # ~120k obs
marg_slope <- sapply(taus, function(t) coef(rq(y ~ x, tau = t, data = big))[2])
tab <- data.frame(tau = taus, conditional = cond_slope(taus), marginal = marg_slope,
                  gap = cond_slope(taus) - marg_slope)
logmsg("conditional vs marginal slopes:")
print(round(tab, 3))
utils::write.csv(tab, file.path(outdir, "cond_vs_marg_slopes.csv"), row.names = FALSE)

## (2) bqmm targets the CONDITIONAL estimand; IJ vs YWH coverage --------------
logmsg("coverage of CONDITIONAL fixed effects (bqmm), reps =", REPS)
covd <- function(ci, truth) truth >= ci[1] && truth <= ci[2]
rows <- list()
for (tau in c(0.25, 0.5, 0.75)) {
  ti <- cond_int(tau); ts <- cond_slope(tau)
  for (r in seq_len(REPS)) {
    d <- sim(G = 40, npg = 12, seed = 1000 + r)
    fit <- tryCatch(suppressWarnings(suppressMessages(
      bqmm(y ~ x + (1 | g), d, tau = tau, chains = 2, iter = 800, warmup = 400,
           cores = 1, seed = 1000 + r, refresh = 0))), error = function(e) NULL)
    if (is.null(fit)) next
    cy <- confint(fit, method = "ywh"); ci_ij <- confint(fit, method = "ij")
    cn <- confint(fit, adjusted = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      tau = tau, rep = r, slope_est = fixef(fit)["x"],
      ywh_int = covd(cy["(Intercept)", ], ti),  ywh_slope = covd(cy["x", ], ts),
      ij_int  = covd(ci_ij["(Intercept)", ], ti), ij_slope = covd(ci_ij["x", ], ts),
      nai_int = covd(cn["(Intercept)", ], ti),  nai_slope = covd(cn["x", ], ts))
  }
  df <- do.call(rbind, rows); utils::write.csv(df, file.path(outdir, "cond_vs_marg_cov.csv"), row.names = FALSE)
  s <- df[df$tau == tau, ]
  ms <- marg_slope[match(tau, taus)]
  logmsg(sprintf("  tau=%.2f cond(int/slope)=%.2f/%.2f marg_slope=%.2f | cover int(YWH/IJ/naive)=%.2f/%.2f/%.2f slope=%.2f/%.2f/%.2f | bqmm slope mean=%.3f",
    tau, ti, ts, ms,
    mean(s$ywh_int), mean(s$ij_int), mean(s$nai_int),
    mean(s$ywh_slope), mean(s$ij_slope), mean(s$nai_slope), mean(s$slope_est)))
}
logmsg("COND VS MARG DONE")
