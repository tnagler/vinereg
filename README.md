
<!-- README.md is generated from README.Rmd. Please edit that file -->
vinereg
=======

[![Travis-CI Build Status](https://travis-ci.org/tnagler/vinereg.svg?branch=master)](https://travis-ci.org/tnagler/vinereg) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/tnagler/vinereg?branch=master&svg=true)](https://ci.appveyor.com/project/tnagler/vinereg) [![Coverage status](https://codecov.io/gh/tnagler/vinereg/branch/master/graph/badge.svg)](https://codecov.io/github/tnagler/vinereg?branch=master) [![CRAN status](https://www.r-pkg.org/badges/version/vinereg)](https://cran.r-project.org/package=vinereg)

An R package for D-vine copula based mean and quantile regression.

How to install
--------------

-   the stable release from CRAN:

``` r
install.packages("vinereg")
```

-   the latest development version:

``` r
# install.packages("devtools")
devtools::install_github("tnagler/vinereg", build_vignettes = TRUE)
```

Example
-------

``` r
library(vinereg)
data(mtcars)

# declare factors and discrete variables
for (var in c("cyl", "vs", "gear", "carb"))
    mtcars[[var]] <- as.ordered(mtcars[[var]])
mtcars[["am"]] <- as.factor(mtcars[["am"]])

# fit model
(fit <- vinereg(mpg ~ ., data = mtcars))
#> D-vine regression model: mpg | disp, cyl, wt, hp, am.1, gear, drat, vs, qsec, carb 
#> nobs = 32, edf = 92, cll = -56.4, caic = 296.8, cbic = 431.65

summary(fit)
#>     var     edf         cll        caic       cbic      p_value
#> 1   mpg 3.2e-11 -98.5927195 197.1854390 197.185439           NA
#> 2  disp 2.0e+00  29.5342816 -55.0685632 -52.137091 1.490817e-13
#> 3   cyl 3.0e+00   3.2125562  -0.4251125   3.972095 9.266313e-02
#> 4    wt 4.0e+00   1.8040933   4.3918134  10.254757 4.616201e-01
#> 5    hp 8.0e+00   1.9098413  12.1803174  23.906205 8.730147e-01
#> 6  am.1 9.0e+00   1.4655935  15.0688130  28.260436 9.669607e-01
#> 7  gear 1.1e+01   1.0700775  19.8598449  35.982940 9.979398e-01
#> 8  drat 1.2e+01   0.5924805  22.8150390  40.403870 9.999637e-01
#> 9    vs 1.0e+01   0.7494010  18.5011979  33.158557 9.989390e-01
#> 10 qsec 1.6e+01   1.4835459  29.0329081  52.484683 9.998425e-01
#> 11 carb 1.7e+01   0.3701619  33.2596762  58.177187 1.000000e+00

# show marginal effects for all selected variables
plot_effects(fit)
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

``` r

# predict mean and median
head(predict(fit, mtcars, alpha = c(NA, 0.5)), 4)
#>       mean      0.5
#> 1 21.05510 20.88661
#> 2 20.89661 20.77428
#> 3 28.01396 28.49054
#> 4 21.26247 21.35846
```

Vignettes
---------

For more examples, have a look at the vignettes with

``` r
vignette("abalone-example", package = "vinereg")
vignette("bike-rental", package = "vinereg")
```

### References

Kraus and Czado (2017). D-vine copula based quantile regression. *Computational Statistics & Data Analysis*, 110, 1-18. [link](https://www.sciencedirect.com/science/article/pii/S0167947316303073), [preprint](https://arxiv.org/abs/1510.04161)

Schallhorn, N., Kraus, D., Nagler, T., Czado, C. (2017). D-vine quantile regression with discrete variables. Working paper, [preprint](https://arxiv.org/abs/1705.08310).
