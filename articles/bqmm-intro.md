# Introduction to bqmm

`bqmm` fits Bayesian mixed-effects quantile regression using the
asymmetric Laplace working likelihood and Stan.

## A first model

``` r

library(bqmm)

fit <- bqmm(distance ~ age + (1 | Subject),
            data = nlme::Orthodont,
            tau  = 0.5)

summary(fit)
```

## Several quantiles at once

Passing a vector of quantiles fits each one and returns a `bqmm_multi`:

``` r

fit3 <- bqmm(distance ~ age + (1 | Subject),
             data = nlme::Orthodont,
             tau  = c(0.1, 0.5, 0.9))

coef(fit3)      # tau-by-coefficient matrix
plot(fit3)      # coefficient-versus-tau paths
```

## Valid inference

By default `bqmm` applies the Yang, Wang and He (2016) correction so
that fixed-effect intervals are asymptotically valid despite the
misspecified asymmetric Laplace likelihood:

``` r

vcov(fit, adjusted = TRUE)
vcov(fit, adjusted = FALSE)   # naive posterior covariance
```

See
[`vignette("bqmm-inference")`](https://kvenkita.github.io/bqmm/articles/bqmm-inference.md)
for details.
