# Variance-estimator bake-off for the mixed-model fixed effects.
# Fits ONE bqmm model per replicate and computes four fixed-effect covariance
# estimators from the same fit, then reports frequentist coverage of each:
#   (1) naive       : posterior covariance Sigma_post
#   (2) koenker     : D1^{-1} G_cluster D1^{-1}/n  (the current package impl)
#   (3) ywh_indep   : Sigma_post G_indep   Sigma_post,  G = (1/s^2) sum x x' psi^2
#   (4) ywh_cluster : Sigma_post G_cluster Sigma_post,  G = (1/s^2) sum_g s_g s_g'
# Theory: the YWH multiplicative form uses Sigma_post (which encodes the
# multilevel structure) as the bread, so it should retain mixed-model variance
# while correcting ALD misspecification, and reduce to ~Sigma_post under correct
# specification. The bake-off decides which achieves nominal coverage.

suppressMessages(library(bqmm, lib.loc = file.path(Sys.getenv("TEMP"), "bqmm_lib")))
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else 60L
outdir <- file.path("tools", "validation"); dir.create(outdir, FALSE, TRUE)
logmsg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n"); flush.console() }

sim_lsm <- function(G = 25, npg = 8, b0 = 1, b1 = 2, a = 1, c = 0, su = 0.7) {
  n <- G * npg; g <- factor(rep(seq_len(G), each = npg)); x <- runif(n, 0.5, 2.5)
  u <- rnorm(G, 0, su)[g]; y <- b0 + b1 * x + u + (a + c * abs(x)) * rnorm(n)
  data.frame(y = y, x = x, g = g)
}
true_coef <- function(tau, b0 = 1, b1 = 2, a = 1, c = 0)
  c(b0 + a * qnorm(tau), b1 + c * qnorm(tau))

# all four estimators from one fit
estimators <- function(fit) {
  db <- rstan::extract(fit$stanfit, pars = "beta", permuted = TRUE)$beta
  db <- as.matrix(db); K <- ncol(db)
  betahat <- apply(db, 2, median)
  Sig <- cov(db)
  s <- median(rstan::extract(fit$stanfit, "sigma")$sigma)
  X <- fit$parsed$X
  bdr <- as.matrix(rstan::extract(fit$stanfit, "b", permuted = TRUE)$b)
  bhat <- apply(bdr, 2, median)
  resid <- as.numeric(fit$parsed$y - X %*% betahat - fit$parsed$Z %*% bhat)
  g <- as.integer(fit$parsed$flist[[1]])
  tau <- fit$tau; n <- nrow(X)
  psi <- tau - (resid < 0)

  # meats
  M_ind <- crossprod(X * psi) / s^2
  s_g <- rowsum(X * psi, group = g)
  M_clu <- crossprod(s_g) / s^2

  # koenker bread (Powell) + cluster meat, current package form
  h <- bqmm:::hall_sheather_bandwidth(n, tau)
  kd <- (abs(resid) < h) / (2 * h)
  D1 <- crossprod(X * sqrt(kd)) / n
  if (qr(D1)$rank < K) { f0 <- mean(dnorm(resid/sd(resid)))/sd(resid); D1 <- f0*crossprod(X)/n }
  D0 <- crossprod(s_g) / n
  V_koen <- solve(D1) %*% D0 %*% solve(D1) / n

  list(
    betahat = betahat,
    naive       = Sig,
    koenker     = V_koen,
    ywh_indep   = Sig %*% M_ind %*% Sig,
    ywh_cluster = Sig %*% M_clu %*% Sig
  )
}

dgps <- list(homoscedastic = c(a = 1.0, c = 0.0),
             heteroscedastic = c(a = 0.6, c = 0.9))
taus <- c(0.25, 0.5, 0.75)
methods <- c("naive", "koenker", "ywh_indep", "ywh_cluster")
rows <- list()

for (dg in names(dgps)) {
  pp <- dgps[[dg]]
  for (tau in taus) {
    for (r in seq_len(REPS)) {
      set.seed(7000 + r)
      d <- sim_lsm(a = pp["a"], c = pp["c"])
      truth <- true_coef(tau, a = pp["a"], c = pp["c"])
      fit <- tryCatch(suppressWarnings(suppressMessages(
        bqmm(y ~ x + (1 | g), d, tau = tau, chains = 2, iter = 700, warmup = 350,
             cores = 1, seed = 7000 + r, refresh = 0, control = list(adapt_delta = 0.9),
             adjust = FALSE))), error = function(e) NULL)
      if (is.null(fit)) next
      es <- estimators(fit)
      for (m in methods) {
        se <- sqrt(diag(es[[m]]))
        lo <- es$betahat - 1.96 * se; hi <- es$betahat + 1.96 * se
        rows[[length(rows) + 1]] <- data.frame(
          dgp = dg, tau = tau, rep = r, method = m,
          cov_int = truth[1] >= lo[1] && truth[1] <= hi[1],
          cov_slope = truth[2] >= lo[2] && truth[2] <= hi[2],
          se_int = se[1], se_slope = se[2])
      }
    }
    df <- do.call(rbind, rows); utils::write.csv(df, file.path(outdir, "bakeoff.csv"), row.names = FALSE)
    sub <- df[df$dgp == dg & df$tau == tau, ]
    agg <- aggregate(cbind(cov_int, cov_slope) ~ method, sub, mean)
    logmsg(sprintf("%s tau=%.2f (n=%d):", dg, tau, sum(sub$method=="naive")))
    for (i in seq_len(nrow(agg)))
      logmsg(sprintf("    %-12s int=%.2f slope=%.2f", agg$method[i], agg$cov_int[i], agg$cov_slope[i]))
  }
}
logmsg("BAKEOFF COMPLETE")
