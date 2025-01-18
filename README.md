<!-- README.md is generated from README.Rmd. Please edit that file -->

# vinereg

[![R-CMD-check](https://github.com/tnagler/vinereg/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/tnagler/vinereg/actions/workflows/R-CMD-check.yaml)
[![Coverage
status](https://codecov.io/gh/tnagler/vinereg/branch/main/graph/badge.svg)](https://app.codecov.io/github/tnagler/vinereg?branch=main)
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
    #> D-vine regression model: mpg | disp, qsec, hp, drat 
    #> nobs = 32, edf = 18.39, cll = -50.08, caic = 136.93, cbic = 163.88

    summary(fit)
    #>    var      edf         cll       caic       cbic      p_value
    #> 1  mpg 0.000000 -100.135440 200.270879 200.270879           NA
    #> 2 disp 8.391335   31.185601 -45.588532 -33.289052 2.446502e-10
    #> 3 qsec 1.624310    3.907191  -4.565762  -2.184953 1.300182e-02
    #> 4   hp 7.371096   11.928452  -9.114713   1.689367 1.576882e-03
    #> 5 drat 1.000000    3.038036  -4.076071  -2.610335 1.370252e-02

    # show marginal effects for all selected variables
    plot_effects(fit)
    #> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />


    # predict mean and median
    head(predict(fit, mtcars, alpha = c(NA, 0.5)), 4)
    #>       mean      0.5
    #> 1 22.46158 22.35297
    #> 2 22.45410 22.35836
    #> 3 24.89114 24.60640
    #> 4 20.44469 20.44982

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
