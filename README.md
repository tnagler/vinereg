
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
(fit <- vinereg(mpg ~ ., data = mtcars))
#> D-vine regression model: mpg | disp, wt, hp, gear, cyl, vs, qsec 
#> nobs = 32, edf = 9, cll = -58.41, caic = 134.82, cbic = 148.01

summary(fit)
#>    var edf         cll        caic        cbic      p_value
#> 1  mpg   0 -100.189867 200.3797334 200.3797334           NA
#> 2 disp   1   27.086917 -52.1738350 -50.7080991 1.835143e-13
#> 3   wt   1    2.676766  -3.3535326  -1.8877967 2.068033e-02
#> 4   hp   1    3.983133  -5.9662654  -4.5005295 4.765716e-03
#> 5 gear   1    1.392314  -0.7846281   0.6811078 9.517278e-02
#> 6  cyl   2    3.116818  -2.2336361   0.6978357 4.429790e-02
#> 7   vs   2    2.458009  -0.9160183   2.0154535 8.560521e-02
#> 8 qsec   1    1.065405  -0.1308095   1.3349264 1.443645e-01

# show marginal effects for all selected variables
plot_effects(fit)
#> `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

``` r

# predict mean and median
head(predict(fit, mtcars, alpha = c(NA, 0.5)), 4)
#>       mean      0.5
#> 1 15.76021 15.76021
#> 2 15.33784 15.33784
#> 3 21.06066 21.06066
#> 4 14.64041 14.64041
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
