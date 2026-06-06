# Check whether quantile predictions are monotone in tau

Check whether quantile predictions are monotone in tau

## Usage

``` r
is_noncrossing(preds)
```

## Arguments

- preds:

  Matrix of fitted quantiles (columns ordered by increasing tau).

## Value

Logical scalar: `TRUE` if every row is non-decreasing.
