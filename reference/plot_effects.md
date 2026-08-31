# Plot marginal effects of a D-vine regression model

The points show fitted conditional quantiles against each requested
variable. For variable \\X_k\\, the smooth curve estimates \\E\[\hat
Q\_\alpha(Y \mid X) \mid X_k = x\]\\. It therefore averages over the
conditional distribution of the other variables. A curve for an
unselected variable can vary when that variable is associated with
selected predictors. The curve is descriptive and is not a
partial-dependence or causal effect.

## Usage

``` r
plot_effects(object, alpha = c(0.1, 0.5, 0.9), vars = object$order)
```

## Arguments

- object:

  a `vinereg` object

- alpha:

  vector of quantile levels.

- vars:

  vector of expanded variable names to display.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
# simulate data
x <- matrix(rnorm(100), 50, 2)
y <- x %*% c(1, -2)
dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))

# fit vine regression model
fit <- vinereg(y ~ ., dat)
plot_effects(fit)
#> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'
```
