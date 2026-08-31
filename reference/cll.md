# Conditional log-likelihood

Calculates the conditional log-likelihood of the response given the
covariates.

## Usage

``` r
cll(object, newdata, cores = 1)
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

The scalar conditional log-likelihood evaluated on `newdata`.

## Examples

``` r
# simulate data
x <- matrix(rnorm(100), 50, 2)
y <- x %*% c(1, -2)
dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))

# fit vine regression model
fit <- vinereg(y ~ ., dat)

cll(fit, dat)
#> [1] 15.8034
fit$stats$cll
#> [1] 15.8034
```
