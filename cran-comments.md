## Submission

This is a new submission of bqmm, providing Bayesian multilevel quantile
regression via the asymmetric Laplace working likelihood and Stan.

## R CMD check results

0 errors | 0 warnings | 3 notes

* Note (CRAN incoming feasibility): "New submission". The package URLs
  (https://github.com/kvenkita/bqmm and https://kvenkita.github.io/bqmm/) become
  live when the repository and its GitHub Pages site are published alongside this
  submission.
* Note: "unable to verify current time" -- a transient network/clock check on
  the build machine, unrelated to the package.
* Note: a few examples run for more than 5 seconds. These are the multi-quantile
  methods (`coef`/`summary`/`plot` for `bqmm_multi`), whose examples each fit
  more than one Markov chain Monte Carlo model with Stan; the cost is inherent to
  the method and comparable to other Stan-based packages (e.g. rstanarm, brms).
  All such examples are wrapped in `\donttest{}`.

A note on installed size may also appear: the package pre-compiles Stan models,
so the compiled `libs/` directory dominates, as for rstanarm and brms. The
package builds with the standard rstantools-generated `src/Makevars` (verified by
a clean `R CMD INSTALL` with no user Makevars); no non-portable compiler flags
are shipped. GNU make is declared in SystemRequirements (inherited from the
rstan/StanHeaders toolchain).

## Test environments

* Windows 11, R 4.4.3 (local)
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), macOS, Windows

## Downstream dependencies

There are currently no downstream dependencies.
