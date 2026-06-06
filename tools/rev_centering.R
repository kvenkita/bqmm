source("R/family.R"); source("R/adjust.R")
options(width=120)
cat("=== Residual centering / score-mean assumption ===\n\n")
cat("The sandwich assumes residuals are evaluated where the tau-quantile is 0,\n")
cat("i.e. (1/n)sum 1{resid<0} ~ tau, so the score sum_i psi_i x_i ~ 0.\n")
cat("In bqmm, resid = y-(Xbeta+Zb) at POSTERIOR MEDIANS, not the QR optimum.\n")
cat("The ALD posterior median need not satisfy the tau-quantile estimating eqn.\n\n")

# Demonstrate: ALD/Bayesian fit median vs rq solution -> different residual sign balance
set.seed(3); n<-400; X<-cbind(1,rnorm(n)); tau<-0.7
y <- 1 + 2*X[,2] + (rnorm(n) - qnorm(tau))  # tau-quantile of error is 0
# pretend 'posterior median beta' is slightly off the rq solution:
suppressMessages(library(quantreg))
fit <- rq.fit.br(X,y,tau=tau); beta_rq<-fit$coefficients
resid_rq <- as.numeric(fit$residuals)
cat("At rq optimum: frac(resid<0)=", mean(resid_rq<0), " (target tau=",tau,")\n")

beta_off <- beta_rq + c(0.15,-0.1)   # mimic posterior-median offset
resid_off <- y - as.numeric(X%*%beta_off)
cat("At offset beta: frac(resid<0)=", mean(resid_off<0), "\n\n")

# Effect on the band-density (D1) and on the score (matters for cluster meat)
h<-hall_sheather_bandwidth(n,tau)
cat("f0 at rq resid :", mean(abs(resid_rq)<h)/(2*h), "\n")
cat("f0 at off resid:", mean(abs(resid_off)<h)/(2*h), "\n")
cat("--> band density is anchored at 0; if resid not centered at the tau-quantile,\n")
cat("    f0 estimates the density at the WRONG point. Mild for small offsets.\n\n")
cat("Also: independence D0 = tau(1-tau)X'X/n ignores resid entirely, so it is\n")
cat("    robust to centering; only D1 (and cluster D0 via psi) depend on it.\n")
