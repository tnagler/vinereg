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

-   the stable release from CRAN:

        install.packages("vinereg")

-   the latest development version:

        # install.packages("devtools")
        devtools::install_github("tnagler/vinereg", build_vignettes = TRUE)

## Functionality

See the [package website](https://tnagler.github.io/vinereg/).

## Example

    set.seed(5)
    library(vinereg)
    data(mtcars)

    # declare factors and discrete variables
    for (var in c("cyl", "vs", "gear", "carb"))
        mtcars[[var]] <- as.ordered(mtcars[[var]])
    mtcars[["am"]] <- as.factor(mtcars[["am"]])

    # fit model
    (fit <- vinereg(mpg ~ ., family = "nonpar", data = mtcars))
    #> D-vine regression model: mpg | disp, qsec, hp 
    #> nobs = 32, edf = 21.86, cll = -55.94, caic = 155.59, cbic = 187.63

    summary(fit)
    #>    var       edf         cll       caic        cbic      p_value
    #> 1  mpg  0.000000 -100.189867 200.379733 200.3797334           NA
    #> 2 disp 11.177711   29.363521 -36.371618 -19.9880453 1.873313e-08
    #> 3 qsec  2.328636    4.167727  -3.678182  -0.2650159 2.180106e-02
    #> 4   hp  8.353178   10.723533  -4.740711   7.5028411 7.480400e-03

    # show marginal effects for all selected variables
    plot_effects(fit)
    #> `geom_smooth()` using method = 'loess' and formula 'y ~ x'

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />


    # predict mean and median
    head(predict(fit, mtcars, alpha = c(NA, 0.5)), 4)
    #>       mean      0.5
    #> 1 22.36600 22.27170
    #> 2 22.18247 22.01755
    #> 3 25.33357 24.90170
    #> 4 20.24950 20.03959

## Vignettes

For more examples, have a look at the vignettes with

    vignette("abalone-example", package = "vinereg")
    vignette("bike-rental", package = "vinereg")

### References

Kraus and Czado (2017). D-vine copula based quantile regression.
*Computational Statistics & Data Analysis*, 110, 1-18.
[link](https://www.sciencedirect.com/science/article/pii/S0167947316303073),
[preprint](https://arxiv.org/abs/1510.04161)

Schallhorn, N., Kraus, D., Nagler, T., Czado, C. (2017). D-vine quantile
regression with discrete variables. Working paper,
[preprint](https://arxiv.org/abs/1705.08310).
