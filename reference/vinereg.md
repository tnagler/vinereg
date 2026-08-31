# D-vine regression models

Sequential estimation of a regression D-vine for the purpose of quantile
prediction as described in Kraus and Czado (2017).

## Usage

``` r
vinereg(
  formula,
  data,
  family_set = "parametric",
  selcrit = "aic",
  order = NA,
  par_1d = list(),
  weights = numeric(),
  cores = 1,
  ...,
  uscale = FALSE
)
```

## Arguments

- formula:

  a two-sided formula containing untransformed variable names. Compute
  transformations and interactions in `data` before fitting.

- data:

  data frame (or object coercible by
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html))
  containing the variables in the model.

- family_set:

  see `family_set` argument of
  [`rvinecopulib::bicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop.html).

- selcrit:

  selection criterion based on conditional log-likelihood. `"aic"`
  (default) and `"bic"` penalize model complexity; `"loglik"` imposes no
  penalty.

- order:

  the order of covariates in the D-vine, provided as vector of variable
  names. An unordered factor name expands to all of its dummy variables
  in level order. Expanded dummy names remain accepted. The order is
  selected automatically if `order = NA` (default).

- par_1d:

  list of options passed to
  [`kde1d::kde1d()`](https://tnagler.github.io/kde1d/reference/kde1d.html),
  must be one value for each margin, e.g. `list(xmin = c(0, 0, NaN))` if
  the response and first covariate have non-negative support.

- weights:

  optional numeric vector of nonnegative observation weights. Supply one
  value per row of `data` or per row of the model frame after missing
  values have been omitted. Missing weights are omitted as well.

- cores:

  integer; the number of cores to use for computations.

- ...:

  further arguments passed to
  [`rvinecopulib::bicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop.html).

- uscale:

  if TRUE, vinereg assumes that marginal distributions have been taken
  care of in a preliminary step.

## Value

An object of class vinereg. It is a list containing the elements

- formula:

  the formula used for the fit.

- selcrit:

  criterion used for variable selection.

- model_frame:

  the data used to fit the regression model.

- margins:

  list of marginal models fitted by
  [`kde1d::kde1d()`](https://tnagler.github.io/kde1d/reference/kde1d.html).

- vine:

  an
  [`rvinecopulib::vinecop_dist()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_dist.html)
  object containing the fitted D-vine.

- stats:

  fit statistics such as conditional log-likelihood/AIC/BIC and p-values
  for each variable's contribution.

- order:

  order of the covariates chosen by the variable selection algorithm.

- selected_vars:

  indices of selected variables.

- factor_map:

  mapping from original variables to expanded variables.

Use
[`predict.vinereg()`](https://tnagler.github.io/vinereg/reference/predict.vinereg.md)
to predict conditional quantiles.
[`summary.vinereg()`](https://tnagler.github.io/vinereg/reference/vinereg-methods.md)
shows the contribution of each selected variable with the associated
p-value derived from a likelihood ratio test.

## Details

If discrete variables are declared as
[`ordered()`](https://rdrr.io/r/base/factor.html) or
[`factor()`](https://rdrr.io/r/base/factor.html), they are handled as
described in Panagiotelis et al. (2012).

## References

Kraus and Czado (2017), D-vine copula based quantile regression,
Computational Statistics and Data Analysis, 110, 1-18

Panagiotelis, A., Czado, C., & Joe, H. (2012). Pair copula constructions
for multivariate discrete data. Journal of the American Statistical
Association, 107(499), 1063-1072.

## See also

[`predict.vinereg`](https://tnagler.github.io/vinereg/reference/predict.vinereg.md)

## Examples

``` r
# simulate data
x <- matrix(rnorm(100), 50, 2)
y <- x %*% c(1, -2)
dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))

# fit vine regression model
(fit <- vinereg(y ~ ., dat))
#> D-vine regression model: y | x.2, x.1
#> nobs = 50, edf = 6.01, cll = 9.3, caic = -6.56, cbic = 4.94

# inspect model
summary(fit)
#>   var      edf        cll       caic      cbic      p_value
#> 1   y 4.014548 -114.12266  236.27441  243.9503           NA
#> 2 x.2 1.000000   43.56331  -85.12662  -83.2146 1.017911e-20
#> 3 x.1 1.000000   79.85597 -157.71193 -155.7999 1.307941e-36
plot_effects(fit)
#> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'


# model predictions
mu_hat <- predict(fit, newdata = dat, alpha = NA) # mean
med_hat <- predict(fit, newdata = dat, alpha = 0.5) # median

# observed vs predicted
plot(cbind(y, mu_hat))


## fixed variable order (no selection)
(fit <- vinereg(y ~ ., dat, order = c("x.2", "x.1", "z")))
#> D-vine regression model: y | x.2, x.1, z.1, z.2
#> nobs = 50, edf = 9.01, cll = 12.85, caic = -7.68, cbic = 9.56
```
