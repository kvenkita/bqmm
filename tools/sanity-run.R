# End-to-end sanity run against the installed package. Fits the Orthodont
# random-intercept model at three quantiles and exercises the S3 surface.
lib <- file.path(Sys.getenv("TEMP"), "bqmm_lib")
library(bqmm, lib.loc = lib)
suppressMessages(library(nlme))

cat("== single quantile (median) ==\n")
fit <- bqmm(distance ~ age + (1 | Subject), data = Orthodont,
            tau = 0.5, chains = 2, iter = 1000, seed = 1,
            cores = 1, refresh = 0)
print(fit)

cat("\n== summary (YWH-adjusted intervals) ==\n")
print(summary(fit))

cat("\n== vcov: adjusted vs naive ==\n")
va <- vcov(fit, adjusted = TRUE)
vn <- vcov(fit, adjusted = FALSE)
cat("adjusted SEs:", round(sqrt(diag(va)), 4), "\n")
cat("naive    SEs:", round(sqrt(diag(vn)), 4), "\n")

cat("\n== VarCorr ==\n")
print(VarCorr(fit))

cat("\n== compare to lqmm point estimates ==\n")
if (requireNamespace("lqmm", quietly = TRUE)) {
  lq <- lqmm::lqmm(fixed = distance ~ age, random = ~ 1, group = Subject,
                   tau = 0.5, data = Orthodont)
  cat("bqmm fixef:", round(fixef(fit), 4), "\n")
  cat("lqmm fixef:", round(as.numeric(lqmm::coef.lqmm(lq)), 4), "\n")
}

cat("\n== multiple quantiles + non-crossing ==\n")
fit3 <- bqmm(distance ~ age + (1 | Subject), data = Orthodont,
             tau = c(0.1, 0.5, 0.9), chains = 2, iter = 1000, seed = 2,
             cores = 1, refresh = 0)
print(class(fit3))
print(coef(fit3))

# predictions across quantiles for a grid of ages, then rearrange
nd <- data.frame(age = seq(8, 14, by = 2))
preds <- sapply(fit3$fits, function(f) predict(f, newdata = nd, re.form = NA))
cat("\npopulation-level predicted quantiles (rows = ages 8..14):\n")
print(round(preds, 3))
cat("non-crossing after rearrangement:\n")
print(round(rearrange_quantiles(preds), 3))

cat("\nALL SANITY CHECKS COMPLETED\n")
