# Standard methods for D-vine regression models

These methods extract the fitted model's call information, data,
effective degrees of freedom, conditional log-likelihood, and
variable-wise summary.

## Usage

``` r
# S3 method for class 'vinereg'
print(x, ...)

# S3 method for class 'vinereg'
summary(object, ...)

# S3 method for class 'vinereg'
logLik(object, ...)

# S3 method for class 'vinereg'
nobs(object, use.fallback = TRUE, ...)

# S3 method for class 'vinereg'
formula(x, ...)

# S3 method for class 'vinereg'
model.frame(formula, ...)
```

## Arguments

- x, object, formula:

  a `vinereg` object.

- ...:

  unused.

- use.fallback:

  unused; included for compatibility with
  [`nobs()`](https://rdrr.io/r/stats/nobs.html).

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with one row for the response and each selected predictor.
[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns an object of
class `logLik`; its `df` attribute contains the effective degrees of
freedom. [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the
number of observations used for fitting.
[`formula()`](https://rdrr.io/r/stats/formula.html) and
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) return the
model formula and frame.
