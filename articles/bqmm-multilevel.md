# Multilevel structure in bqmm

`bqmm` uses `lme4`’s formula grammar, so random effects are written
inline and nested *or* crossed structures come for free.

## Random intercepts

``` r

library(bqmm)
bqmm(y ~ x + (1 | group), data, tau = 0.5)
```

Each group gets its own intercept deviation `u_j ~ N(0, σ_u²)`.
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) returns
the posterior-median deviations;
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) returns `σ_u`.

## Random slopes

``` r

bqmm(y ~ x + (1 + x | group), data, tau = 0.5)            # diagonal
bqmm(y ~ x + (1 + x | group), data, tau = 0.5,
     cov = "unstructured")                                # correlated
```

With `cov = "diagonal"` (the default) the intercept and slope deviations
are independent. With `cov = "unstructured"` they share an
LKJ-correlated covariance and
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) carries the
correlation matrix:

``` r

fit <- bqmm(y ~ x + (1 + x | group), data, tau = 0.5, cov = "unstructured")
VarCorr(fit)
attr(VarCorr(fit), "correlation")
```

`cov = "unstructured"` currently supports a **single** grouping factor.
Use the default diagonal covariance for multiple or crossed terms.

## Nested and crossed grouping

``` r

bqmm(y ~ x + (1 | school/classroom), data, tau = 0.5)    # nested
bqmm(y ~ x + (1 | school) + (1 | neighbourhood), data, tau = 0.5)  # crossed
```

Both are parsed by `lme4` and handled by the diagonal model — no special
syntax is needed. The variance-component mapping (which random-effect
column belongs to which `(term, coefficient)`) is built directly from
`lme4::mkReTrms()` and is verified against `lme4`’s own design matrices
in the package tests.

## Practical notes

- Variance components and correlations need **many groups** to be
  estimable (≥ 20-30). With few groups the estimates shrink toward the
  prior; this is a feature of the data, not a defect.
- When random-effect SDs are small relative to the asymmetric-Laplace
  scale (which has SD ≈ 2.83·σ at τ = 0.5), the random structure is
  weakly identified.
- See
  [`vignette("bqmm-inference")`](https://kvenkita.github.io/bqmm/articles/bqmm-inference.md)
  for how the multilevel structure feeds into the fixed-effect variance
  correction.
