source("R/family.R")
set.seed(7)

## Confirm closed-form variance used in check is the TRUE ALD variance,
## by computing Var analytically from the density via numeric integration
## (independent of the RNG and independent of the formula).
cat("== independent check: variance from density integration ==\n")
mom <- function(k, mu, sigma, tau)
  integrate(function(x) x^k * dald(x, mu, sigma, tau), -Inf, Inf, rel.tol=1e-10)$value
for (tau in c(0.3, 0.7)) {
  m1 <- mom(1, 0, 2.5, tau)
  m2 <- mom(2, 0, 2.5, tau)
  v_int <- m2 - m1^2
  v_cf  <- 2.5^2 * (1 - 2*tau + 2*tau^2) / ((1-tau)^2 * tau^2)
  cat(sprintf("tau=%.1f: var(integration)=%.5f var(closed-form)=%.5f\n", tau, v_int, v_cf))
}

## Analytic derivation of the mixture variance to check sigma factoring.
## y = mu + sigma*(theta*w + k*sqrt(w)*z), w~Exp(1), z~N(0,1), indep.
## Var contribution (drop mu, sigma): Var(theta*w) + E[k^2 * w] (since z indep, mean 0)
## = theta^2 * Var(w) + k^2 * E[w] = theta^2 * 1 + k^2 * 1  (Exp(1): mean=var=1)
## with theta=(1-2tau)/(tau(1-tau)), k^2 = 2/(tau(1-tau)).
## Then Var(y) = sigma^2 * (theta^2 + k^2).
cat("\n== analytic mixture var vs ALD closed form ==\n")
for (tau in c(0.1,0.3,0.5,0.7,0.9)) {
  theta <- (1-2*tau)/(tau*(1-tau)); k2 <- 2/(tau*(1-tau))
  v_mix <- 2.5^2 * (theta^2 + k2)
  v_cf  <- 2.5^2 * (1 - 2*tau + 2*tau^2) / ((1-tau)^2 * tau^2)
  cat(sprintf("tau=%.1f: sigma^2(theta^2+k^2)=%.6f  ALD var=%.6f  diff=%.2e\n",
              tau, v_mix, v_cf, v_mix - v_cf))
}
