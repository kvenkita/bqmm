## Submission

This is a new submission of bqmm, providing Bayesian multilevel quantile
regression via the asymmetric Laplace working likelihood and Stan.

## R CMD check results

0 errors | 0 warnings | 1 note

* Note (CRAN incoming feasibility): "New submission", and the package URLs
  (https://github.com/kvenkita/bqmm and https://kvenkita.github.io/bqmm/) are
  reported as not yet reachable. These become live when the repository and its
  GitHub Pages site are published alongside this submission.
* Packages that pre-compile Stan models may also raise an installed-size note;
  this is expected and matches rstanarm and brms.
* GNU make is a SystemRequirement, declared in DESCRIPTION (inherited from the
  rstan/StanHeaders toolchain).
* The package compiles with the standard rstantools-generated `src/Makevars`
  (verified by a clean `R CMD INSTALL` with no user Makevars); no non-portable
  compiler flags are shipped.

## Test environments

* Windows 11, R 4.4.3 (local)
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), macOS, Windows

## Downstream dependencies

There are currently no downstream dependencies.
