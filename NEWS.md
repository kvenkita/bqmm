# bqmm 0.1.0

First release.

* `bqmm()` — Bayesian mixed-effects quantile regression via the asymmetric
  Laplace likelihood and Stan, with an `lme4`-style formula interface
  (`y ~ x + (1 + x | group)`) supporting nested and crossed random effects.
* `cov = "diagonal"` (default) or `cov = "unstructured"` — the latter adds
  LKJ-correlated random effects for a single grouping factor; `VarCorr()`
  reports the posterior-median correlation matrix.
* Fitting one or several quantiles in a single call (a vector `tau`).
* Post-hoc non-crossing rearrangement of fitted quantiles
  (Chernozhukov, Fernandez-Val and Galichon, 2010).
* Yang, Wang and He (2016) posterior-variance correction for valid
  frequentist inference, exposed through `vcov(object, adjusted = TRUE)`.
* `lme4`/`rstanarm`-style S3 methods and integration with the `posterior`
  and `bayesplot` ecosystems.
