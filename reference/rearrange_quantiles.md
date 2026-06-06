# Rearrange fitted quantiles to remove crossing

Post-hoc monotonisation of estimated quantile curves (Chernozhukov,
Fernandez-Val and Galichon, 2010): at any covariate point, sorting the
fitted values across quantile levels into increasing order yields a
valid, non-crossing set of quantiles and never increases estimation
error. This is the v0.1 non-crossing strategy; joint constrained
estimation is deferred.

## Usage

``` r
rearrange_quantiles(preds)
```

## Arguments

- preds:

  A numeric matrix of fitted quantiles with one column per quantile
  level, ordered by increasing `tau` (rows = observations).

## Value

A matrix of the same shape with each row sorted increasingly.

## Examples

``` r
m <- rbind(c(1, 0.5, 2), c(0, 1, 0.8))   # some crossings
rearrange_quantiles(m)
#>      [,1] [,2] [,3]
#> [1,]  0.5  1.0    2
#> [2,]  0.0  0.8    1
```
