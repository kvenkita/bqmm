source("R/family.R"); source("R/adjust.R")
options(width=120)
cat("=== Dimensional analysis of Sigma_post fallback (line 100) ===\n\n")
cat("Real sandwich:    V = D1^{-1} D0 D1^{-1} / n\n")
cat("  D0 = tau(1-tau) (X'X)/n        [O(1) in n, since X'X ~ n]\n")
cat("  D1 = f0 (X'X)/n                [O(1) in n]\n")
cat("  => D1^{-1} ~ O(1/n) * (avg gram)^{-1};  V ~ O(1/n). correct.\n\n")
cat("Fallback:         V = Sigma_post %*% D0 %*% Sigma_post\n")
cat("  Sigma_post ~ posterior cov of beta ~ O(1/n)\n")
cat("  D0 ~ O(1) (it is the AVERAGED gram times tau(1-tau))\n")
cat("  => V_fallback ~ O(1/n) * O(1) * O(1/n) = O(1/n^2).  WRONG ORDER!\n\n")

cat("Numeric demonstration: build a case where D1 is nonsingular,\n")
cat("compute the TRUE sandwich, then ALSO compute the fallback expression\n")
cat("on the same D0 and a realistic Sigma_post, and compare magnitudes.\n\n")

set.seed(11)
n <- 1500
X <- cbind(1, rnorm(n), rnorm(n))
tau <- 0.5
# realistic residuals + a 'posterior cov' that mimics ALD naive cov.
resid <- rald(n, 0, 1, tau)
sw <- compute_ywh_sandwich(X, resid, tau=tau)           # real sandwich
V_real <- sw$vcov
D0 <- sw$D0

# A plausible naive posterior cov: ALD posterior var ~ (some const) * (X'X)^{-1}
# Use the actual asymptotic naive ALD cov: sigma^2 * (2/(tau(1-tau))) ...
# Simpler: just take Sigma_post = c*(X'X)^{-1} for c making it O(1/n)-correct.
Sigma_post <- 1.0 * solve(crossprod(X))   # O(1/n) magnitude
V_fallback <- Sigma_post %*% D0 %*% Sigma_post

cat("diag(V_real)    :", diag(V_real), "\n")
cat("diag(V_fallback):", diag(V_fallback), "\n")
cat("ratio fallback/real:", diag(V_fallback)/diag(V_real), "\n\n")

cat("Now scale n x10 (n=15000) and see if the ratio changes by ~10x\n")
cat("(it should, proving fallback is off by a factor ~n):\n")
n2 <- 15000
X2 <- cbind(1, rnorm(n2), rnorm(n2))
resid2 <- rald(n2,0,1,tau)
sw2 <- compute_ywh_sandwich(X2, resid2, tau=tau)
Sp2 <- 1.0*solve(crossprod(X2))
Vf2 <- Sp2 %*% sw2$D0 %*% Sp2
cat("n=1500  ratio:", mean(diag(V_fallback)/diag(V_real)), "\n")
cat("n=15000 ratio:", mean(diag(Vf2)/diag(sw2$vcov)), "\n")
