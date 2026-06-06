source("R/family.R")

cat("== rho at u=0 ==\n")
cat("rho_tau(0, 0.3) =", rho_tau(0, 0.3), "\n")  # (0<0)=FALSE -> 0*(0.3-0)=0

cat("\n== edge tau enforcement ==\n")
print(tryCatch(dald(1, tau = 0),   error = function(e) conditionMessage(e)))
print(tryCatch(dald(1, tau = 1),   error = function(e) conditionMessage(e)))
print(tryCatch(dald(1, sigma = 0), error = function(e) conditionMessage(e)))
print(tryCatch(rald(5, tau = 1),   error = function(e) conditionMessage(e)))
print(tryCatch(rald(5, sigma = -1),error = function(e) conditionMessage(e)))

cat("\n== tau near 0/1 numerical behavior of rald ==\n")
set.seed(1)
for (t in c(1e-3, 1e-6, 1-1e-6)) {
  y <- rald(1e5, 0, 1, t)
  cat(sprintf("tau=%.1e: any NaN/Inf? %s ; P(Y<=0)=%.4f\n",
              t, any(!is.finite(y)), mean(y <= 0)))
}

cat("\n== dald log handling at extreme x ==\n")
cat("dald(1e6, 0,1,0.5, log=TRUE) =", dald(1e6,0,1,0.5,log=TRUE), "\n")
cat("dald(-1e6,0,1,0.5, log=TRUE) =", dald(-1e6,0,1,0.5,log=TRUE), "\n")
cat("dald(1e6,0,1,0.5) =", dald(1e6,0,1,0.5), "(underflow to 0, expected)\n")

cat("\n== vectorized sigma check in dald: dald allows vector sigma? ==\n")
## family.R dald: any(sigma<=0) guard, u=(x-mu)/sigma vectorizes
print(dald(c(0,1,2), mu=0, sigma=c(1,2,3), tau=0.5))
