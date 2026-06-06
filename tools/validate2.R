# Confirmation validation against the INSTALLED package after the YWH fix.
# Verifies that the package wiring (ywh_adjust -> vcov -> confint, default
# meat = "cluster") reproduces the bake-off coverage end-to-end, and compares
# fixed effects to lqmm.
suppressMessages(library(bqmm, lib.loc = file.path(Sys.getenv("TEMP"), "bqmm_lib")))
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else 80L
outdir <- file.path("tools", "validation"); dir.create(outdir, FALSE, TRUE)
logmsg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n"); flush.console() }

sim_lsm <- function(G = 25, npg = 8, b0 = 1, b1 = 2, a = 1, c = 0, su = 0.7) {
  n <- G * npg; g <- factor(rep(seq_len(G), each = npg)); x <- runif(n, 0.5, 2.5)
  u <- rnorm(G, 0, su)[g]; y <- b0 + b1 * x + u + (a + c * abs(x)) * rnorm(n)
  data.frame(y = y, x = x, g = g)
}
true_coef <- function(tau, b0 = 1, b1 = 2, a = 1, c = 0)
  c(b0 + a * qnorm(tau), b1 + c * qnorm(tau))

dgps <- list(homoscedastic = c(a = 1, c = 0), heteroscedastic = c(a = 0.6, c = 0.9))
taus <- c(0.25, 0.5)
rows <- list()
for (dg in names(dgps)) {
  pp <- dgps[[dg]]
  for (tau in taus) {
    for (r in seq_len(REPS)) {
      set.seed(9000 + r)
      d <- sim_lsm(a = pp["a"], c = pp["c"])
      truth <- true_coef(tau, a = pp["a"], c = pp["c"])
      fit <- tryCatch(suppressWarnings(suppressMessages(
        bqmm(y ~ x + (1 | g), d, tau = tau, chains = 2, iter = 700, warmup = 350,
             cores = 1, seed = 9000 + r, refresh = 0, control = list(adapt_delta = 0.9)))),
        error = function(e) NULL)
      if (is.null(fit)) next
      ca <- confint(fit, adjusted = TRUE); cn <- confint(fit, adjusted = FALSE)
      rows[[length(rows) + 1]] <- data.frame(dgp = dg, tau = tau, rep = r,
        adj_int = truth[1] >= ca[1,1] && truth[1] <= ca[1,2],
        adj_slope = truth[2] >= ca[2,1] && truth[2] <= ca[2,2],
        nai_int = truth[1] >= cn[1,1] && truth[1] <= cn[1,2],
        nai_slope = truth[2] >= cn[2,1] && truth[2] <= cn[2,2])
    }
    df <- do.call(rbind, rows); utils::write.csv(df, file.path(outdir, "validate2.csv"), row.names = FALSE)
    s <- df[df$dgp == dg & df$tau == tau, ]
    logmsg(sprintf("%s tau=%.2f (n=%d): adjusted int/slope=%.2f/%.2f  naive=%.2f/%.2f",
      dg, tau, nrow(s), mean(s$adj_int), mean(s$adj_slope), mean(s$nai_int), mean(s$nai_slope)))
  }
}

# lqmm comparison
if (requireNamespace("lqmm", quietly = TRUE)) {
  logmsg("lqmm comparison (random intercept, tau=0.5), reps =", min(REPS, 40))
  cmp <- list()
  for (r in seq_len(min(REPS, 40L))) {
    set.seed(9500 + r)
    d <- sim_lsm(a = 1, c = 0)
    fb <- tryCatch(suppressWarnings(bqmm(y ~ x + (1 | g), d, 0.5, chains = 2, iter = 700,
      warmup = 350, cores = 1, seed = 9500 + r, refresh = 0)), error = function(e) NULL)
    fl <- tryCatch(lqmm::lqmm(y ~ x, random = ~ 1, group = g, tau = 0.5, data = d),
      error = function(e) NULL)
    if (is.null(fb) || is.null(fl)) next
    bb <- fixef(fb); bl <- as.numeric(lqmm::coef.lqmm(fl))
    cmp[[length(cmp) + 1]] <- data.frame(rep = r, d_int = bb[1] - bl[1], d_slope = bb[2] - bl[2])
  }
  cdf <- do.call(rbind, cmp); utils::write.csv(cdf, file.path(outdir, "compare.csv"), row.names = FALSE)
  logmsg(sprintf("lqmm: mean|bqmm-lqmm| int=%.3f slope=%.3f (sd slope of estimates ~ ref)",
    mean(abs(cdf$d_int)), mean(abs(cdf$d_slope))))
}
logmsg("VALIDATE2 COMPLETE")
