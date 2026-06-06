# Does rstan::extract(permuted=TRUE) use a FIXED permutation across calls?
# If two separate "beta" extracts are identical, the permutation is stored in
# the fit and the old separate-call code was aligned (fragile but correct).
# If they differ, separate beta/b extracts were genuinely scrambled.
suppressMessages(library(bqmm, lib.loc = file.path(Sys.getenv("TEMP"), "bqmm_lib")))
set.seed(1)
g <- factor(rep(1:6, each = 6)); x <- rnorm(36)
d <- data.frame(y = 1 + 2 * x + rnorm(6)[g] + rnorm(36), x = x, g = g)
fit <- suppressWarnings(bqmm(y ~ x + (1 | g), d, tau = 0.5,
                             chains = 1, iter = 400, warmup = 200, seed = 7, refresh = 0))

b1 <- rstan::extract(fit$stanfit, pars = "beta", permuted = TRUE)$beta
b2 <- rstan::extract(fit$stanfit, pars = "beta", permuted = TRUE)$beta
cat("two separate 'beta' extracts identical:", identical(b1, b2), "\n")

# within a single call, do beta/b/sigma share ordering? (always yes) — sanity
ex <- rstan::extract(fit$stanfit, pars = c("beta", "b", "sigma"), permuted = TRUE)
cat("single-call dims  beta:", paste(dim(ex$beta), collapse = "x"),
    " b:", paste(dim(ex$b), collapse = "x"),
    " sigma:", length(ex$sigma), "\n")

# new posterior_epred coherence: epred should equal X*beta + Z*b per draw
ep <- posterior_epred(fit)
cat("posterior_epred dims:", paste(dim(ep), collapse = "x"),
    " (expect S x N =", nrow(b1), "x", nrow(d), ")\n")
# reconstruct draw 1 by hand from a single joint extract and compare
loc_manual <- as.numeric(ex$beta[1, ] %*% t(fit$parsed$X) + ex$b[1, ] %*% t(fit$parsed$Z))
ep2 <- bqmm:::bqmm_location_draws(fit)$loc
cat("hand-rebuilt draw1 matches bqmm_location_draws draw1:",
    isTRUE(all.equal(loc_manual, as.numeric(ep2[1, ]))), "\n")
cat("DONE\n")
