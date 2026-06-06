# Recovery of the random-effect CORRELATION (and SDs, betas) for the correlated
# model. Simulates a random intercept + slope with a known correlation rho and
# ALD errors (correct spec at tau=0.5), fits cov="unstructured", and checks that
# the posterior recovers rho. Uses the separate library bqmm_lib2.
suppressMessages({
  library(bqmm, lib.loc = file.path(Sys.getenv("TEMP"), "bqmm_lib2"))
  library(MASS)
})
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else 40L
outdir <- file.path("tools", "validation"); dir.create(outdir, FALSE, TRUE)
logmsg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n"); flush.console() }

# Well-identified DGP (RE SDs comparable to noise; enough groups) so the RE
# covariance and correlation are estimable. See tools/diag-corr.R.
b0 <- 1; b1 <- 2; s0 <- 1.5; s1 <- 1.0; rho <- 0.5; sigma <- 0.5; tau <- 0.5
G <- 60; npg <- 12; n <- G * npg
Su <- matrix(c(s0^2, rho*s0*s1, rho*s0*s1, s1^2), 2)

rows <- list()
for (r in seq_len(REPS)) {
  set.seed(6000 + r)
  g <- factor(rep(seq_len(G), each = npg)); x <- rnorm(n)
  U <- MASS::mvrnorm(G, c(0, 0), Su)
  y <- b0 + b1 * x + U[g, 1] + U[g, 2] * x + bqmm:::rald(n, 0, sigma, tau)
  d <- data.frame(y = y, x = x, g = g)
  fit <- tryCatch(suppressWarnings(suppressMessages(
    bqmm(y ~ x + (1 + x | g), d, tau = tau, cov = "unstructured",
         chains = 2, iter = 2000, warmup = 1000, cores = 1, seed = 6000 + r,
         refresh = 0, control = list(adapt_delta = 0.97)))),
    error = function(e) NULL)
  if (is.null(fit)) next
  bb <- fixef(fit); vc <- VarCorr(fit); cm <- attr(vc, "correlation")
  rows[[length(rows) + 1]] <- data.frame(
    rep = r, b0 = bb[1], b1 = bb[2],
    sd0 = vc[1], sd1 = vc[2], rho = cm[1, 2])
  if (r %% 5 == 0) {
    df <- do.call(rbind, rows); utils::write.csv(df, file.path(outdir, "corr_recovery.csv"), row.names = FALSE)
    logmsg("rep", r, "of", REPS)
  }
}
df <- do.call(rbind, rows); utils::write.csv(df, file.path(outdir, "corr_recovery.csv"), row.names = FALSE)
logmsg(sprintf("CORR RECOVERY (%d reps): b0=%.3f (t %.1f) b1=%.3f (t %.1f) sd0=%.3f (t %.2f) sd1=%.3f (t %.2f) rho=%.3f (t %.2f)",
  nrow(df), mean(df$b0), b0, mean(df$b1), b1, mean(df$sd0), s0, mean(df$sd1), s1, mean(df$rho), rho))
logmsg("CORR RECOVERY COMPLETE")
