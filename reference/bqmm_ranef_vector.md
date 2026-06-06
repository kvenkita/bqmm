# Posterior-median random effects as a Q-vector aligned with the columns of Z

Handles both the diagonal model (`b` drawn as an S x Q matrix) and the
correlated model (`b` drawn as an S x M x L array). For the correlated
model the column order `(level - 1) * M + coef` matches the columns of
`Z`.

## Usage

``` r
bqmm_ranef_vector(object)
```
