# Changelog

## vinereg 0.13.0

DEPENDS

- Adapted the C++ interface to remain compatible with newer
  `rvinecopulib` releases and their C++17 requirements.

BUG FIXES

- Stabilized a platform-sensitive test for fixed variable orders.

- Forwarded `cores` when predicting conditional means.

- Improved validation of variable orders and prediction data.

IMPROVEMENTS

- Formula transformations and interactions now produce an informative
  error; transformed predictors can be computed in `data` before
  fitting.

- Observation weights are aligned with rows omitted because of missing
  data.

- Original factor names can be used in fixed variable orders. Expanded
  dummy names remain supported for backwards compatibility.

NEW FEATURES

- Added [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`formula()`](https://rdrr.io/r/stats/formula.html), and
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) methods
  for fitted models. Standard
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) methods are now available.

DOCUMENTATION

- Clarified the supported formula and prediction interfaces, conditional
  density outputs, and the interpretation of marginal effect plots.

## vinereg 0.12.0

DEPENDS

- Adapts to new and now required version of rvinecopulib (\>=
  0.7.1.1.0).

## vinereg 0.11.0

CRAN release: 2025-01-11

DEPENDS

- Adapts to new and now required version of kde1d (\>= 1.1.0).

BUG FIX

- Fixed conditional density computations for discrete variables
  ([\#17](https://github.com/tnagler/vinereg/issues/17), thanks
  [@FLBuchner](https://github.com/FLBuchner).)

## vinereg 0.10.0

CRAN release: 2024-02-09

NEW FEATURE

- New function
  [`cpdf()`](https://tnagler.github.io/vinereg/reference/cpdf.md) to
  compute the conditional density or probability mass.

## vinereg 0.9.2

CRAN release: 2023-11-02

BUG FIX

- add compiler flag to prevent boost/functional from using
  `unary_function`.

## vinereg 0.9.1

CRAN release: 2023-10-23

BUG FIX

- fix unnecessary error when calling
  [`vinereg()`](https://tnagler.github.io/vinereg/reference/vinereg.md)
  with weights.

## vinereg 0.9.0

NEW FEATURE

- New function
  [`cll()`](https://tnagler.github.io/vinereg/reference/cll.md) to
  compute the conditional log-likelihood.

## vinereg 0.8.3

CRAN release: 2022-03-23

BUG FIX

- avoid bit-wise operations on boolean variables (fixes warnings for
  clang\>=14).

## vinereg 0.8.2

CRAN release: 2022-03-15

DEPENDS

- require recent version of vinecopulib (\>= 0.6.1.1.2) to ensure
  compatible RcppThread versions.

## vinereg 0.8.1

CRAN release: 2022-01-22

BUG FIXES

- [`require()`](https://rdrr.io/r/base/library.html) calls with single
  argument in vignettes.

## vinereg 0.8.0

CRAN release: 2021-11-23

BUG FIXES

- fix [`cpit()`](https://tnagler.github.io/vinereg/reference/cpit.md)
  (last conditioning was sometimes omitted).

- prevent `rvinecopulib` from spawning own threads.

NEW FEATURES

- add `uscale` option to allow for external marginal modeling.

## vinereg 0.7.3

CRAN release: 2021-02-24

BUG FIXES

- fix simulated data size in documentation examples.

- properly handle case where no covariates are selected.

- conditional use of packages in Suggests.

## vinereg 0.7.2

CRAN release: 2020-11-17

BUG FIXES

- remove bias in quantiles for discrete variables.

## vinereg 0.7.1

CRAN release: 2020-10-12

This is a maintenance release following an update in rvinecopulib.

DEPENDS

- requires rvinecopulib (\>= 0.5.4.1.0) to fix an uninitialized value
  issue.

NEW FEATURES

- variables generated from `factor`s are now named with the
  corresponding factor level.

## vinereg 0.7.0

CRAN release: 2019-12-02

DEPENDENCIES

- removed dependence on future and furrr packages.

NEW FEATURES

- faster runtimes, especially for parallelized code.

- handle discrete variables properly with both parametric and
  nonparametric pair-copulas.

REMOVED FEATURES

- removed support for `uscale` argument.

## vinereg 0.6.0

NEW FEATURES

- new function
  [`cpit()`](https://tnagler.github.io/vinereg/reference/cpit.md) to
  compute the conditional cdf.

## vinereg 0.5.0

CRAN release: 2018-08-09

DEPENDS

- require rvinecopulib (\>= 0.3.0) due to breaking changes in this
  package.

BUG FIXES

- prevent nan errors in loglik calculation.

- allow for empty and bivariate models.

- properly pass degree parameter for marginal estimation.

## vinereg 0.4.0

BUG FIXES

- Fix handling of `uscale` in
  [`fitted.vinereg()`](https://tnagler.github.io/vinereg/reference/predict.vinereg.md).

- Fix handling of `mult` parameter for pair-copula fits in
  [`vinereg()`](https://tnagler.github.io/vinereg/reference/vinereg.md).

- Fix orientation of asymmetric pair-copulas.

## vinereg 0.3.0

CRAN release: 2018-06-02

DEPENDS

- Use `furrr` and `future` packages instead of `parallel`, `doParallel`,
  and `foreach` for parallelization.

NEW FEATURES

- New [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) generics for
  `vinereg` objects.

- New
  [`plot_effects()`](https://tnagler.github.io/vinereg/reference/plot_effects.md)
  method to show the marginal effects of variables.

- Allow to predict the mean with `predict(object, alpha = NA)`.

## vinereg 0.2.0

CRAN release: 2018-05-29

- First official release.
