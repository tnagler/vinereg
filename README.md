<!-- README.md is generated from README.Rmd. Please edit that file -->

# vinereg

[![R build
status](https://github.com/tnagler/vinereg/workflows/R-CMD-check/badge.svg)](https://github.com/tnagler/vinereg)
<!-- [![Coverage status](https://codecov.io/gh/tnagler/vinereg/branch/main/graph/badge.svg)](https://app.codecov.io/github/tnagler/vinereg?branch=main) -->
[![CRAN
status](https://www.r-pkg.org/badges/version/vinereg)](https://cran.r-project.org/package=vinereg)

An R package for D-vine copula based mean and quantile regression.

## How to install

-   the stable release from CRAN:

        install.packages("vinereg")

-   the latest development version:

        # install.packages("remotes")
        remotes::install_github("tnagler/vinereg", build_vignettes = TRUE)

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
    #> D-vine regression model: mpg | wt, qsec, drat, gear 
    #> nobs = 32, edf = 23.63, cll = -55.86, caic = 158.98, cbic = 193.62

    summary(fit)
    #>    var       edf         cll       caic        cbic      p_value
    #> 1  mpg  0.000000 -100.135440 200.270879 200.2708794           NA
    #> 2   wt 11.452248   28.706110 -34.507723 -17.7217520 4.161832e-08
    #> 3 qsec  6.091637    7.596924  -3.010573   5.9181583 1.990142e-02
    #> 4 drat  5.089693    5.742895  -1.306405   6.1537401 4.494112e-02
    #> 5 gear  1.000000    2.232423  -2.464845  -0.9991094 3.459922e-02

    # show marginal effects for all selected variables
    plot_effects(fit)
    #> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />


    # predict mean and median
    head(predict(fit, mtcars, alpha = c(NA, 0.5)), 4)
    #>       mean      0.5
    #> 1 23.38467 23.04676
    #> 2 22.69125 22.36638
    #> 3 26.29842 26.10553
    #> 4 20.62143 20.63283

## Vignettes

For more examples, have a look at the vignettes with

    vignette("abalone-example", package = "vinereg")
    vignette("bike-rental", package = "vinereg")

### References

Kraus and Czado (2017). D-vine copula based quantile regression.
*Computational Statistics & Data Analysis*, 110, 1-18.
[link](https://doi.org/10.1016/j.csda.2016.12.009),
[preprint](https://arxiv.org/abs/1510.04161)

Schallhorn, N., Kraus, D., Nagler, T., Czado, C. (2017). D-vine quantile
regression with discrete variables. Working paper,
[preprint](https://arxiv.org/abs/1705.08310).
