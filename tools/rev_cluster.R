source("R/family.R"); source("R/adjust.R")
suppressMessages(library(quantreg))
options(width=120)

cat("=== TASK 3: cluster-robust meat formula ===\n\n")
cat("Code: psi = tau - (resid<0); s_g = rowsum(X*psi, groups); D0 = crossprod(s_g)/n\n")
cat("Standard cluster meat for QR estimating eqn sum_i psi_i x_i = 0:\n")
cat("  J = (1/n) sum_g (sum_{i in g} psi_i x_i)(...)'  = crossprod(s_g)/n.  MATCHES.\n")
cat("psi already encodes tau(1-tau): E[psi^2]=tau(1-tau) under correct quantile,\n")
cat("so NOT multiplying by tau(1-tau) again is CORRECT (would be double counting).\n\n")

# numeric: independence meat should ~equal cluster meat with 1-obs clusters
set.seed(1); n<-300; X<-cbind(1,rnorm(n)); resid<-rald(n,0,1,0.5); tau<-0.5
D0_ind <- tau*(1-tau)*crossprod(X)/n
psi <- tau-(resid<0)
# singleton groups
g <- 1:n
s_g <- rowsum(X*psi, g)
D0_clust_singleton <- crossprod(s_g)/n
cat("Singleton-cluster meat vs independence meat (should be close, psi^2~tau(1-tau)):\n")
cat(" indep D0[1,1]:", D0_ind[1,1], "  singleton-cluster D0[1,1]:", D0_clust_singleton[1,1],"\n")
cat(" psi^2 mean:", mean(psi^2), " tau(1-tau):", tau*(1-tau), "\n")
cat(" (independence meat uses tau(1-tau) exactly; cluster uses empirical psi^2,\n")
cat("  these differ by sampling noise + cross terms = expected.)\n\n")

cat("Cross-check cluster meat against sandwich::vcovCL logic on a linear proxy?\n")
cat("Direct: build correlated clusters, compare to block-bootstrap Var(betahat).\n\n")

# Block bootstrap validation of cluster sandwich
set.seed(99)
G <- 60; m <- 20; n <- G*m
grp <- rep(1:G, each=m)
ranef <- rnorm(G, 0, 1.2)[grp]          # cluster random intercept
x <- rnorm(n)
X <- cbind(1, x)
tau <- 0.5
y <- 1 + 2*x + ranef + rnorm(n)         # within-cluster correlated errors
fit <- rq.fit.br(X, y, tau=tau)
resid <- as.numeric(fit$residuals)
sw_ind <- compute_ywh_sandwich(X, resid, tau=tau)$vcov
sw_cl  <- compute_ywh_sandwich(X, resid, tau=tau, groups=grp)$vcov

# cluster (block) bootstrap truth
B<-400; bb<-matrix(NA,B,2)
for (b in 1:B){
  gid <- sample(1:G, G, replace=TRUE)
  idx <- unlist(lapply(gid, function(g) which(grp==g)))
  fb <- rq.fit.br(X[idx,], y[idx], tau=tau)
  bb[b,]<-fb$coefficients
}
V_boot <- var(bb)
cat("diag independence sandwich:", diag(sw_ind), "\n")
cat("diag cluster   sandwich   :", diag(sw_cl), "\n")
cat("diag cluster bootstrap    :", diag(V_boot), "\n")
cat("ratio clusterSand/boot    :", diag(sw_cl)/diag(V_boot), "\n")
cat("ratio indepSand/boot      :", diag(sw_ind)/diag(V_boot), "\n")

cat("\n=== TASK 5: analytic single-coefficient sandwich ===\n")
# X = intercept only, errors N(0,s). V = tau(1-tau)/(n f0^2). f0=dnorm(0)/s.
set.seed(7); n<-5000; s<-2; tau<-0.5
X<-matrix(1,n,1); resid<-rnorm(n,0,s)   # already centered at quantile
sw<-compute_ywh_sandwich(X,resid,tau=tau)
f0_true<-dnorm(0)/s
V_analytic<-tau*(1-tau)/(n*f0_true^2)
h<-sw$bandwidth; f0_pow<-mean(abs(resid)<h)/(2*h)
V_with_powell<-tau*(1-tau)/(n*f0_pow^2)
cat("sandwich vcov:", as.numeric(sw$vcov), "\n")
cat("analytic (powell f0):", V_with_powell, " -> match:", isTRUE(all.equal(as.numeric(sw$vcov),V_with_powell,tol=1e-8)),"\n")
cat("analytic (true f0)  :", V_analytic, "  ratio sand/true:", as.numeric(sw$vcov)/V_analytic,"\n")
