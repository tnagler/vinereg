# Conditional probability integral transform

Evaluates the fitted conditional distribution of the response given the
covariates. For a discrete response, this is the conditional CDF at the
observed category, not a randomized probability integral transform.

## Usage

``` r
cpit(object, newdata, cores = 1)
```

## Arguments

- object:

  an object of class `vinereg`.

- newdata:

  a data frame containing the response and covariates from the original
  formula, with matching classes and factor levels.

- cores:

  integer; the number of cores to use for computations.

## Value

A numeric vector containing one conditional CDF value per row of
`newdata`.

## Examples

``` r
# simulate data
x <- matrix(rnorm(100), 50, 2)
y <- x %*% c(1, -2)
dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))

# fit vine regression model
fit <- vinereg(y ~ ., dat)

hist(cpit(fit, dat)) # should be approximately uniform
```
