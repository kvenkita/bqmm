# Comprehensive empirical validation of bqmm (long-running). Writes incremental
# CSVs to tools/validation/ so partial results survive interruption.
#
# Studies:
#  (1) RECOVERY   - correct ALD multilevel DGP; bias of fixed effects + variance
#                   components, and coverage of BOTH interval types (~nominal).
#  (2) COVERAGE   - MISSPECIFIED normal-error multilevel DGPs (homoscedastic and
#                   heteroscedastic) with KNOWN tau-varying true coefficients;
#                   frequentist coverage of naive vs Yang-Wang-He intervals.
#                   This is the headline test of the package's central claim.
#  (3) COMPARE    - bqmm random-intercept fixef vs lqmm on identical data.
#
# Usage: Rscript tools/validate.R [reps]   (default reps below)

suppressMessages({
  lib <- file.path(Sys.getenv("TEMP"), "bqmm_lib")
  library(bqmm, lib.loc = lib)
})
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else 100L

outdir <- file.path("tools", "validation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
flush_csv <- function(df, file) utils::write.csv(df, file.path(outdir, file), row.names = FALSE)
logmsg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n"); flush.console() }

# Sampler settings tuned for many fits: the YWH interval depends on the
# posterior MEDIAN + a frequentist sandwich, which stabilise quickly.
CHAINS <- 2L; ITER <- 700L; WARMUP <- 350L
fit_q <- function(form, data, tau, seed) {
  suppressWarnings(suppressMessages(
    bqmm(form, data, tau = tau, chains = CHAINS, iter = ITER, warmup = WARMUP,
         cores = 1L, seed = seed, refresh = 0,
         control = list(adapt_delta = 0.9))
  ))
}

# ---- shared DGP: 2-level location-scale model ------------------------------
# y_ij = b0 + b1 x_ij + u_j + (a + c|x_ij|) * z_ij     (z standard normal)
# Marginal fixed tau-quantile coefficients (u_j has mean 0):
#   intercept(tau) = b0 + a * qnorm(tau)
#   slope(tau)     = b1 + c * qnorm(tau)
sim_lsm <- function(G = 25, npg = 8, b0 = 1, b1 = 2, a = 1, c = 0, su = 0.7,
                    errfun = function(n) rnorm(n)) {
  n <- G * npg
  g <- factor(rep(seq_len(G), each = npg))
  x <- runif(n, 0.5, 2.5)
  u <- rnorm(G, 0, su)[g]
  y <- b0 + b1 * x + u + (a + c * abs(x)) * errfun(n)
  data.frame(y = y, x = x, g = g)
}
true_coef <- function(tau, b0 = 1, b1 = 2, a = 1, c = 0) {
  c(`(Intercept)` = b0 + a * stats::qnorm(tau), x = b1 + c * stats::qnorm(tau))
}

covered <- function(ci, truth) truth >= ci[, 1] & truth <= ci[, 2]

# ============================================================ (1) RECOVERY ===
logmsg("STUDY 1: parameter recovery under correct ALD DGP, reps =", REPS)
rec <- list()
sigma_true <- 1.0; su_true <- 0.7; tau_rec <- 0.5
for (r in seq_len(REPS)) {
  set.seed(1000 + r)
  d <- sim_lsm(su = su_true, a = 1, c = 0,
               errfun = function(n) bqmm:::rald(n, 0, sigma_true, tau_rec))
  # under correct spec, error tau-quantile is 0 => true fixed = (b0, b1)
  truth <- c(`(Intercept)` = 1, x = 2)
  fit <- tryCatch(fit_q(y ~ x + (1 | g), d, tau_rec, seed = 1000 + r),
                  error = function(e) NULL)
  if (is.null(fit)) next
  beta <- fixef(fit)
  vc   <- tryCatch(VarCorr(fit), error = function(e) NA)
  sig  <- stats::median(rstan::extract(fit$stanfit, "sigma")$sigma)
  ci_a <- confint(fit, adjusted = TRUE); ci_n <- confint(fit, adjusted = FALSE)
  rec[[length(rec) + 1]] <- data.frame(
    rep = r,
    b0 = beta[1], b1 = beta[2],
    su_hat = as.numeric(vc[1]), sigma_hat = sig,
    cov_a_int = covered(ci_a, truth)[1], cov_a_slope = covered(ci_a, truth)[2],
    cov_n_int = covered(ci_n, truth)[1], cov_n_slope = covered(ci_n, truth)[2]
  )
  if (r %% 10 == 0) { flush_csv(do.call(rbind, rec), "recovery.csv"); logmsg("  recovery rep", r) }
}
rec_df <- do.call(rbind, rec); flush_csv(rec_df, "recovery.csv")
logmsg("STUDY 1 done. bias b0 =", round(mean(rec_df$b0) - 1, 3),
       " bias b1 =", round(mean(rec_df$b1) - 2, 3),
       " mean sigma_hat =", round(mean(rec_df$sigma_hat), 3),
       " mean su_hat =", round(mean(rec_df$su_hat), 3))

# ============================================================ (2) COVERAGE ===
logmsg("STUDY 2: coverage under MISSPECIFIED normal-error DGPs, reps =", REPS)
dgps <- list(
  homoscedastic = list(a = 1.0, c = 0.0),
  heteroscedastic = list(a = 0.6, c = 0.9)
)
taus <- c(0.25, 0.5, 0.75)
cov_rows <- list()
for (dg in names(dgps)) {
  pp <- dgps[[dg]]
  for (tau in taus) {
    for (r in seq_len(REPS)) {
      set.seed(2000 + r)
      d <- sim_lsm(a = pp$a, c = pp$c, su = 0.7, errfun = function(n) rnorm(n))
      truth <- true_coef(tau, a = pp$a, c = pp$c)
      fit <- tryCatch(fit_q(y ~ x + (1 | g), d, tau, seed = 2000 + r),
                      error = function(e) NULL)
      if (is.null(fit)) next
      ci_a <- confint(fit, adjusted = TRUE); ci_n <- confint(fit, adjusted = FALSE)
      ca <- covered(ci_a, truth); cn <- covered(ci_n, truth)
      cov_rows[[length(cov_rows) + 1]] <- data.frame(
        dgp = dg, tau = tau, rep = r,
        cov_adj_int = ca[1], cov_adj_slope = ca[2],
        cov_nai_int = cn[1], cov_nai_slope = cn[2],
        wid_adj_slope = ci_a[2, 2] - ci_a[2, 1],
        wid_nai_slope = ci_n[2, 2] - ci_n[2, 1]
      )
    }
    cdf <- do.call(rbind, cov_rows); flush_csv(cdf, "coverage.csv")
    sub <- cdf[cdf$dgp == dg & cdf$tau == tau, ]
    logmsg(sprintf("  %s tau=%.2f: adj cover (int/slope)=%.2f/%.2f  naive=%.2f/%.2f  (n=%d)",
                   dg, tau, mean(sub$cov_adj_int), mean(sub$cov_adj_slope),
                   mean(sub$cov_nai_int), mean(sub$cov_nai_slope), nrow(sub)))
  }
}
logmsg("STUDY 2 done.")

# ============================================================= (3) COMPARE ===
if (requireNamespace("lqmm", quietly = TRUE)) {
  logmsg("STUDY 3: bqmm vs lqmm (random intercept), reps =", min(REPS, 40))
  cmp <- list()
  for (r in seq_len(min(REPS, 40L))) {
    set.seed(3000 + r)
    d <- sim_lsm(a = 1, c = 0, su = 0.7, errfun = function(n) rnorm(n))
    fb <- tryCatch(fit_q(y ~ x + (1 | g), d, 0.5, seed = 3000 + r), error = function(e) NULL)
    fl <- tryCatch(lqmm::lqmm(y ~ x, random = ~ 1, group = g, tau = 0.5, data = d),
                   error = function(e) NULL)
    if (is.null(fb) || is.null(fl)) next
    bb <- fixef(fb); bl <- as.numeric(lqmm::coef.lqmm(fl))
    cmp[[length(cmp) + 1]] <- data.frame(rep = r,
      bqmm_int = bb[1], bqmm_slope = bb[2], lqmm_int = bl[1], lqmm_slope = bl[2])
  }
  cmp_df <- do.call(rbind, cmp); flush_csv(cmp_df, "compare.csv")
  logmsg("STUDY 3 done. mean |bqmm-lqmm| int =",
         round(mean(abs(cmp_df$bqmm_int - cmp_df$lqmm_int)), 3),
         " slope =", round(mean(abs(cmp_df$bqmm_slope - cmp_df$lqmm_slope)), 3))
}

logmsg("ALL VALIDATION STUDIES COMPLETE")
