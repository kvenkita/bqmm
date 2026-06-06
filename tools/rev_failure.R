source("R/family.R")
source("R/adjust.R")
suppressMessages(library(quantreg))
options(width=120)

cat("=== TASK 4: failure modes ===\n\n")

## (a) No residuals in band -> normal-density fallback
cat("(a) Empty band fallback:\n")
set.seed(1); n <- 50
X <- cbind(1, rnorm(n))
resid <- rald(n, 0, 1, 0.5)
# force tiny bandwidth so band is empty
sw <- compute_ywh_sandwich(X, resid, tau = 0.5, bandwidth = 1e-12)
cat("  all kdens zero -> fallback used. vcov finite?", all(is.finite(sw$vcov)),
    " PSD?", all(eigen(sw$vcov,only.values=TRUE)$values > -1e-8), "\n")
# BUT: is the fallback triggered correctly? D1 is computed BEFORE the all-zero check.
# Check: when SOME but the matrix still ends up rank-deficient? See (b).

## (b) Partial band: only 1 residual in band with K=2 -> D1 rank 1 -> singular
cat("\n(b) Few-in-band -> singular bread -> Sigma_post fallback:\n")
set.seed(7); n <- 40
X <- cbind(1, rnorm(n))
resid <- rald(n, 0, 1, 0.5)
# bandwidth that lets exactly 1 residual in
ord <- sort(abs(resid))
h <- (ord[1]+ord[2])/2     # exactly 1 inside
nin <- sum(abs(resid) < h)
cat("  residuals in band:", nin, "\n")
Sp <- diag(c(0.02,0.02))
sw <- compute_ywh_sandwich(X, resid, tau=0.5, bandwidth=h, Sigma_post=Sp)
cat("  D1 rank:", qr(sw$D1)$rank, " used Sigma_post fallback (vcov!=sandwich form)\n")
cat("  Sigma_post fallback gives V = Sp D0 Sp; D0 here = tau(1-tau)X'X/n (units mismatch!)\n")
print(sw$vcov)
cat("  Note: D0 has scale ~O(n)*... while Sp ~ posterior. Dimensions of fallback:\n")
cat("  V_fallback = Sp %*% D0 %*% Sp ; D0=tau(1-tau)crossprod(X)/n\n\n")

## (c) collinear X
cat("(c) Collinear X (exact duplicate column):\n")
set.seed(3); n <- 200
x <- rnorm(n)
X <- cbind(1, x, x)   # x repeated
resid <- rald(n,0,1,0.5)
sw <- tryCatch(compute_ywh_sandwich(X, resid, tau=0.5),
               error=function(e) {cat("  ERROR:",conditionMessage(e),"\n"); NULL})
if (!is.null(sw)) {
  cat("  D1 rank:", qr(sw$D1)$rank, "/", ncol(X), "; vcov finite?", all(is.finite(sw$vcov)),"\n")
  cat("  (no Sigma_post supplied -> should it stop? )\n")
}

## (d) tiny n
cat("\n(d) Tiny n=3, K=2:\n")
set.seed(5)
X <- cbind(1, rnorm(3)); resid <- rald(3,0,1,0.5)
sw <- tryCatch(compute_ywh_sandwich(X, resid, tau=0.5, Sigma_post=diag(2)),
               error=function(e){cat("  ERROR:",conditionMessage(e),"\n");NULL})
if(!is.null(sw)) cat("  vcov finite?", all(is.finite(sw$vcov)), "\n")

## (e) Does the all-zero guard actually fire? It recomputes D1 but check the cond.
cat("\n(e) all(kdens==0) guard: is 'in_band' strict '<' consistent? boundary?\n")
cat("    uniform kernel uses abs(resid) < h (strict). fine.\n")
