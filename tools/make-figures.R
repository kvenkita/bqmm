# Generate the primer figures into vignettes/figures/. Uses the default-library
# install. Run: Rscript tools/make-figures.R
suppressMessages({
  library(bqmm); library(ggplot2)
})
figdir <- file.path("vignettes", "figures"); dir.create(figdir, FALSE, TRUE)
theme_set(theme_minimal(base_size = 12))

data(Orthodont, package = "nlme")

## --- Figure 1: coefficient-versus-tau path ---------------------------------
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
fitq <- bqmm(distance ~ age + (1 | Subject), Orthodont, tau = taus,
             chains = 2, iter = 1500, warmup = 750, cores = 1, seed = 1,
             refresh = 0)
rows <- do.call(rbind, lapply(fitq$fits, function(f) {
  ci <- confint(f, adjusted = TRUE)["age", ]
  data.frame(tau = f$tau, est = fixef(f)["age"], lo = ci[1], hi = ci[2])
}))
p1 <- ggplot(rows, aes(tau, est)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#2c7fb8") +
  geom_line(colour = "#2c7fb8", linewidth = 0.8) +
  geom_point(colour = "#2c7fb8", size = 2) +
  labs(x = expression(tau), y = "Effect of age on distance",
       title = "Coefficient-versus-quantile path",
       subtitle = "Effect of age at each conditional quantile (95% adjusted intervals)")
ggsave(file.path(figdir, "coef-path.png"), p1, width = 6.5, height = 4, dpi = 120)
cat("wrote coef-path.png\n")

## --- Figure 2: coverage (from the bake-off results) ------------------------
cov <- data.frame(
  cell = rep(c("homo .25","homo .50","homo .75","het .25","het .50","het .75"), 2),
  method = rep(c("naive", "adjusted (YWH)"), each = 6),
  coverage = c(0.91, 0.90, 0.85, 0.90, 0.90, 0.87,   # naive (slope)
               0.98, 0.98, 0.97, 0.95, 0.95, 0.98))  # adjusted (slope)
cov$cell <- factor(cov$cell, levels = unique(cov$cell))
p2 <- ggplot(cov, aes(cell, coverage, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey30") +
  scale_fill_manual(values = c("naive" = "#bdbdbd", "adjusted (YWH)" = "#2c7fb8")) +
  coord_cartesian(ylim = c(0.7, 1)) +
  labs(x = NULL, y = "Frequentist coverage (slope, nominal 0.95)",
       fill = NULL, title = "Interval coverage: naive vs. adjusted",
       subtitle = "Dashed line = nominal 95%. Naive under-covers; YWH is at/above nominal.") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(figdir, "coverage.png"), p2, width = 6.5, height = 4, dpi = 120)
cat("wrote coverage.png\n")

## --- Figure 3: posterior predictive check ----------------------------------
fit <- bqmm(distance ~ age + (1 | Subject), Orthodont, tau = 0.5,
            chains = 2, iter = 1500, warmup = 750, cores = 1, seed = 2, refresh = 0)
yrep <- posterior_predict(fit)
if (requireNamespace("bayesplot", quietly = TRUE)) {
  p3 <- bayesplot::ppc_dens_overlay(Orthodont$distance, yrep[1:60, ]) +
    labs(title = "Posterior predictive check",
         subtitle = "Observed outcome density (dark) vs. draws from the fitted model")
  ggsave(file.path(figdir, "ppcheck.png"), p3, width = 6.5, height = 4, dpi = 120)
  cat("wrote ppcheck.png\n")
}
cat("FIGURES DONE\n")
