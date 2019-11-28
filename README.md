
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vinereg

[![R build
status](https://github.com/tnagler/vinereg/workflows/R-CMD-check/badge.svg)](https://github.com/tnagler/vinereg)
[![Coverage
status](https://codecov.io/gh/tnagler/vinereg/branch/master/graph/badge.svg)](https://codecov.io/github/tnagler/vinereg?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/vinereg)](https://cran.r-project.org/package=vinereg)

An R package for D-vine copula based mean and quantile regression.

## How to install

  - the stable release from CRAN:
    
    ``` r
    install.packages("vinereg")
    ```

  - the latest development version:
    
    ``` r
    # install.packages("devtools")
    devtools::install_github("tnagler/vinereg", build_vignettes = TRUE)
    ```

## Functionality

See the [package website](https://tnagler.github.io/vinereg).

## Example

``` r
set.seed(5)
library(vinereg)
data(mtcars)

# declare factors and discrete variables
for (var in c("cyl", "vs", "gear", "carb"))
    mtcars[[var]] <- as.ordered(mtcars[[var]])
mtcars[["am"]] <- as.factor(mtcars[["am"]])

# fit model
(fit <- vinereg(mpg ~ ., family = "nonpar", data = mtcars))
#> D-vine regression model: mpg | wt, qsec, drat 
#> nobs = 32, edf = 19.33, cll = -56.94, caic = 152.55, cbic = 180.88

summary(fit)
#>    var       edf         cll       caic       cbic      p_value
#> 1  mpg  0.000000 -100.189867 200.379733 200.379733           NA
#> 2   wt 10.597257   29.409422 -37.624331 -22.091551 1.062440e-08
#> 3 qsec  5.739523    7.867765  -4.256484   4.156141 1.286392e-02
#> 4 drat  2.996463    5.973303  -5.953681  -1.561657 7.542848e-03

# show marginal effects for all selected variables
plot_effects(fit)
#> `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

``` r

# predict mean and median
head(predict(fit, mtcars, alpha = c(NA, 0.5)), 4)
#>       mean      0.5
#> 1 23.33594 22.56025
#> 2 22.30219 21.68133
#> 3 25.71678 25.58609
#> 4 20.21699 20.36290
```

## Vignettes

For more examples, have a look at the vignettes with

``` r
vignette("abalone-example", package = "vinereg")
vignette("bike-rental", package = "vinereg")
```

### References

Kraus and Czado (2017). D-vine copula based quantile regression.
*Computational Statistics & Data Analysis*, 110, 1-18.
[link](https://www.sciencedirect.com/science/article/pii/S0167947316303073),
[preprint](https://arxiv.org/abs/1510.04161)

Schallhorn, N., Kraus, D., Nagler, T., Czado, C. (2017). D-vine quantile
regression with discrete variables. Working paper,
[preprint](https://arxiv.org/abs/1705.08310).
