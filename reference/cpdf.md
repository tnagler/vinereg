# Conditional density or probability mass

Calculates the conditional density of a continuous response or
conditional probability mass of a discrete response given the
covariates.

## Usage

``` r
cpdf(object, newdata, cores = 1)
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

A numeric vector containing one conditional density or probability mass
per row of `newdata`.

## Examples

``` r
# simulate data
x <- matrix(rnorm(100), 50, 2)
y <- x %*% c(1, -2)
dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))

# fit vine regression model
fit <- vinereg(y ~ ., dat)

cpdf(fit, dat)
#>  [1] 1.5899905 2.9633043 1.8086191 1.4032837 3.7545453 0.1809221 2.9756595
#>  [8] 1.5163262 1.7251766 1.9226761 1.0518828 2.5011815 1.9666177 0.2826336
#> [15] 1.4040541 1.1666493 0.2007590 2.0568615 2.1251529 0.3286866 1.8159674
#> [22] 1.6989435 0.9833907 0.2306088 1.3958616 1.0258675 3.0231066 1.0156668
#> [29] 1.5222920 2.6585333 1.2492579 2.0367193 2.6061828 0.2785767 1.1123120
#> [36] 1.2254587 2.8384611 2.7749277 1.4788255 2.1267000 3.0413040 3.3521151
#> [43] 2.4249209 1.6870749 0.8124165 1.8739937 0.8799563 1.8888287 0.6037526
#> [50] 1.7563609
```
