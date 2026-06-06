# Permutation-behavior check only (works against the currently-installed pkg).
suppressMessages(library(bqmm, lib.loc = file.path(Sys.getenv("TEMP"), "bqmm_lib")))
set.seed(1)
g <- factor(rep(1:6, each = 6)); x <- rnorm(36)
d <- data.frame(y = 1 + 2 * x + rnorm(6)[g] + rnorm(36), x = x, g = g)
fit <- suppressWarnings(bqmm(y ~ x + (1 | g), d, tau = 0.5,
                             chains = 1, iter = 400, warmup = 200, seed = 7, refresh = 0))
b1 <- rstan::extract(fit$stanfit, pars = "beta", permuted = TRUE)$beta
b2 <- rstan::extract(fit$stanfit, pars = "beta", permuted = TRUE)$beta
cat("two separate 'beta' extracts identical:", identical(b1, b2), "\n")
# Also: are beta and b from separate calls aligned? Check a stored permutation.
cat("has stored sim$permutation:",
    !is.null(fit$stanfit@sim$permutation), "\n")
