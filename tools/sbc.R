# Simulation-Based Calibration (Talts et al. 2018) for the diagonal bqmm model.
# For each simulation: draw theta* from the model's (fixed) priors, simulate ALD
# multilevel data, fit with the SAME prior, and rank theta* among thinned
# posterior draws. Under a correctly-implemented posterior the ranks are
# Uniform; systematic deviation flags a sampler/likelihood/prior bug.
#
# Usage: Rscript tools/sbc.R [nsims]
suppressMessages(library(bqmm, lib.loc = file.path(Sys.getenv("TEMP"), "bqmm_lib")))
args <- commandArgs(trailingOnly = TRUE)
NSIM <- if (length(args) >= 1) as.integer(args[1]) else 300L
outdir <- file.path("tools", "validation"); dir.create(outdir, FALSE, TRUE)
logmsg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n"); flush.console() }

tau <- 0.5
# FIXED priors (not data-scaled) so we can draw theta* from them exactly.
PB_MEAN <- 0; PB_SD <- 2; P_SIGMA <- 1; P_RE <- 1
pri <- bqmm_prior(beta_mean = PB_MEAN, beta_sd = PB_SD,
                  sigma_scale = P_SIGMA, re_scale = P_RE)

# Fixed design across simulations (standard SBC).
set.seed(2024)
G <- 20; npg <- 10; n <- G * npg
gf <- factor(rep(seq_len(G), each = npg))
xx <- rnorm(n)
Xd <- cbind(1, xx)                      # intercept + slope (K = 2)
Zlev <- as.integer(gf)

L <- 256L                               # thinned posterior draws per sim
rank_rows <- list()
params <- c("b0", "b1", "sigma", "sd_re")

for (s in seq_len(NSIM)) {
  set.seed(5000 + s)
  # draw theta* from priors
  b_star  <- rnorm(2, PB_MEAN, PB_SD)
  sig_star <- abs(rnorm(1, 0, P_SIGMA))          # half-normal
  sdre_star <- abs(rnorm(1, 0, P_RE))            # half-normal
  u_star <- rnorm(G, 0, sdre_star)               # random intercepts
  mu <- as.numeric(Xd %*% b_star) + u_star[Zlev]
  y <- mu + bqmm:::rald(n, 0, sig_star, tau)
  d <- data.frame(y = y, x = xx, g = gf)

  fit <- tryCatch(suppressWarnings(suppressMessages(
    bqmm(y ~ x + (1 | g), d, tau = tau, prior = pri, adjust = FALSE,
         chains = 2, iter = 1200, warmup = 400, cores = 1, seed = 5000 + s,
         refresh = 0, control = list(adapt_delta = 0.95)))),
    error = function(e) NULL)
  if (is.null(fit)) next

  ex <- rstan::extract(fit$stanfit, pars = c("beta", "sigma", "sd_re"),
                       permuted = TRUE)
  S <- length(ex$sigma)
  idx <- round(seq(1, S, length.out = L))        # thin to L draws
  draws <- cbind(
    b0 = ex$beta[idx, 1], b1 = ex$beta[idx, 2],
    sigma = ex$sigma[idx], sd_re = as.matrix(ex$sd_re)[idx, 1]
  )
  truth <- c(b0 = b_star[1], b1 = b_star[2], sigma = sig_star, sd_re = sdre_star)
  ranks <- vapply(params, function(p) sum(draws[, p] < truth[p]), integer(1))
  rank_rows[[length(rank_rows) + 1]] <- as.data.frame(c(list(sim = s), as.list(ranks)))

  if (s %% 20 == 0) {
    rdf <- do.call(rbind, rank_rows)
    utils::write.csv(rdf, file.path(outdir, "sbc_ranks.csv"), row.names = FALSE)
    logmsg("sim", s, "of", NSIM)
  }
}

rdf <- do.call(rbind, rank_rows)
utils::write.csv(rdf, file.path(outdir, "sbc_ranks.csv"), row.names = FALSE)

# Uniformity test: bin ranks (0..L) into B bins, chi-square GOF vs uniform.
B <- 16L
logmsg("SBC uniformity (chi-square GOF of rank histogram, ", nrow(rdf), " sims):")
for (p in params) {
  r <- rdf[[p]]
  bins <- cut(r, breaks = seq(0, L + 1, length.out = B + 1), include.lowest = TRUE)
  obs <- as.numeric(table(bins))
  chi <- suppressWarnings(stats::chisq.test(obs))
  logmsg(sprintf("   %-7s chisq p = %.3f  (mean rank %.1f, expected %.1f)",
                 p, chi$p.value, mean(r), L / 2))
}
logmsg("SBC COMPLETE")
