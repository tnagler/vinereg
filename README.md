
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

Functionality
-------------

See the [package website](https://tnagler.github.io/vinereg).

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
#> D-vine regression model: mpg | disp, hp, gear, carb, wt, vs, drat, qsec, cyl, am.1 
#> nobs = 32, edf = 92, cll = -51.59, caic = 287.18, cbic = 422.02

summary(fit)
#>     var     edf         cll       caic       cbic      p_value
#> 1   mpg 3.2e-11 -98.5927195 197.185439 197.185439           NA
#> 2  disp 2.0e+00  29.5342816 -55.068563 -52.137091 1.490817e-13
#> 3    hp 3.0e+00   2.3323113   1.335377   5.732585 1.980680e-01
#> 4  gear 4.0e+00   2.6417293   2.716541   8.579485 2.594294e-01
#> 5  carb 8.0e+00   4.3416393   7.316721  19.042609 3.697145e-01
#> 6    wt 7.0e+00   1.6074265  10.785147  21.045298 8.644408e-01
#> 7    vs 7.0e+00   1.3034954  11.393009  21.653161 9.188274e-01
#> 8  drat 1.2e+01   1.4321371  21.135726  38.724557 9.964287e-01
#> 9  qsec 1.5e+01   1.2227564  27.554487  49.540526 9.998896e-01
#> 10  cyl 1.7e+01   1.8528062  30.294388  55.211898 9.996929e-01
#> 11 am.1 1.7e+01   0.7366123  32.526775  57.444286 9.999997e-01

# show marginal effects for all selected variables
plot_effects(fit)
#> `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

``` r

# predict mean and median
head(predict(fit, mtcars, alpha = c(NA, 0.5)), 4)
#>       mean      0.5
#> 1 20.68862 20.66071
#> 2 20.50411 20.49184
#> 3 25.70142 25.66380
#> 4 20.99347 20.97258
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
