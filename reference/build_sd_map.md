# Map random-effect columns to variance components

lme4 orders the columns of `t(Zt)` block by term (delimited by `Gp`),
and within a term block the coefficients vary fastest within each level.
This function turns that layout into an integer `sd_map` assigning each
column to a `(term, coefficient)` variance component, plus
human-readable labels.

## Usage

``` r
build_sd_map(cnms, Gp)
```

## Arguments

- cnms:

  Named list of coefficient names per grouping factor (from
  `lme4::mkReTrms()`).

- Gp:

  Integer vector of block boundaries (from `mkReTrms`).

## Value

List with `sd_map`, `components` (labels) and `terms` (metadata).
