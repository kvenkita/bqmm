# Diagnostic single fit for the correlated model: higher RE signal-to-noise,
# more groups, long well-converged chains. Prints convergence + recovery to
# distinguish weak-identification/under-convergence from a model bug.
suppressMessages({ library(bqmm, lib.loc = file.path(Sys.getenv("TEMP"), "bqmm_lib2")); library(MASS) })

b0 <- 1; b1 <- 2; s0 <- 1.5; s1 <- 1.0; rho <- 0.5; sigma <- 0.5; tau <- 0.5
G <- 60; npg <- 12; n <- G * npg
set.seed(11)
g <- factor(rep(seq_len(G), each = npg)); x <- rnorm(n)
Su <- matrix(c(s0^2, rho*s0*s1, rho*s0*s1, s1^2), 2)
U <- MASS::mvrnorm(G, c(0, 0), Su)
y <- b0 + b1 * x + U[g, 1] + U[g, 2] * x + bqmm:::rald(n, 0, sigma, tau)
d <- data.frame(y = y, x = x, g = g)

fit <- bqmm(y ~ x + (1 + x | g), d, tau = tau, cov = "unstructured",
            chains = 4, iter = 3000, warmup = 1500, cores = 1, seed = 11,
            refresh = 0, control = list(adapt_delta = 0.99, max_treedepth = 12))

sm <- rstan::summary(fit$stanfit, pars = c("beta", "sigma", "sd_re", "Corr"))$summary
cat("=== convergence (Rhat / n_eff) ===\n")
print(round(sm[, c("mean", "2.5%", "97.5%", "n_eff", "Rhat")], 3))
cat("\ndivergences:", rstan::get_num_divergent(fit$stanfit), "\n")
cat("\n=== recovery (posterior median) ===\n")
cat("fixef:", round(fixef(fit), 3), " (true", b0, b1, ")\n")
vc <- VarCorr(fit); cat("sd_re:", round(vc, 3), " (true", s0, s1, ")\n")
cat("rho:", round(attr(vc, "correlation")[1, 2], 3), " (true", rho, ")\n")
cat("sigma:", round(median(rstan::extract(fit$stanfit, 'sigma')$sigma), 3), " (true", sigma, ")\n")
cat("\nDIAG DONE\n")
