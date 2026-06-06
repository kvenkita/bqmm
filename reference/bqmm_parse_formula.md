# Parse an lme4-style mixed-model formula

Thin wrapper around
[`lme4::glFormula()`](https://rdrr.io/pkg/lme4/man/modular.html) that
extracts everything `bqmm` needs: the fixed-effect design matrix `X`,
the random-effect design matrix `Z` (dense), the response `y`, and a
mapping from each column of `Z` to a variance component. Reusing lme4's
parser means nested *and* crossed random effects are handled for free.

## Usage

``` r
bqmm_parse_formula(formula, data, na.action = stats::na.omit, contrasts = NULL)
```

## Arguments

- formula:

  A model formula such as `y ~ x + (1 + x | group)`.

- data:

  A data frame.

- na.action, contrasts:

  Passed through to model-frame construction.

## Value

A list with elements:

- y:

  numeric response vector.

- X:

  fixed-effect design matrix (N x K).

- Z:

  random-effect design matrix (N x Q), dense.

- sd_map:

  integer vector (length Q) mapping each Z column to a variance
  component in `1:G`.

- re_components:

  character labels for the `G` variance components.

- re_terms:

  per-term metadata: grouping factor, coefficient names, number of
  levels.

- cnms, flist, Gp:

  the raw lme4 random-effect structures.

- fixed_names:

  column names of `X`.
