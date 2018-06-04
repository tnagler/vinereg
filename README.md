
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
set.seed(5)
library(vinereg)
data(mtcars)

# declare factors and discrete variables
for (var in c("cyl", "vs", "gear", "carb"))
    mtcars[[var]] <- as.ordered(mtcars[[var]])
mtcars[["am"]] <- as.factor(mtcars[["am"]])

# fit model
(fit <- vinereg(mpg ~ ., data = mtcars))
#> D-vine regression model: mpg | disp, hp, gear, carb, cyl, am.1, wt, vs, qsec, drat 
#> nobs = 32, edf = 95, cll = -56.09, caic = 302.19, cbic = 441.43

summary(fit)
#>     var     edf          cll       caic       cbic      p_value
#> 1   mpg 3.2e-11 -98.59271949 197.185439 197.185439           NA
#> 2  disp 2.0e+00  29.53428159 -55.068563 -52.137091 1.490817e-13
#> 3    hp 3.0e+00   2.33231128   1.335377   5.732585 1.980680e-01
#> 4  gear 5.0e+00   2.16379041   5.672419  13.001099 5.032784e-01
#> 5  carb 7.0e+00   2.25663902   9.486722  19.746873 7.191178e-01
#> 6   cyl 9.0e+00   1.57187334  14.856253  28.047876 9.583219e-01
#> 7  am.1 1.0e+01   1.75059329  16.498813  31.156172 9.670581e-01
#> 8    wt 1.3e+01   1.62623809  22.747524  41.802091 9.968597e-01
#> 9    vs 1.3e+01   0.52958111  24.940838  43.995405 9.999946e-01
#> 10 qsec 1.6e+01   0.70430954  30.591381  54.043155 9.999992e-01
#> 11 drat 1.7e+01   0.02886845  33.942263  58.859773 1.000000e+00

# show marginal effects for all selected variables
plot_effects(fit)
#> `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

``` r

# predict mean and median
head(predict(fit, mtcars, alpha = c(NA, 0.5)), 4)
#>       mean      0.5
#> 1 19.34836 19.36129
#> 2 19.17641 19.19810
#> 3 25.28064 25.13942
#> 4 19.70841 19.67779
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
