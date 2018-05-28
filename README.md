# vinereg
[![Travis-CI Build Status](https://travis-ci.org/tnagler/vinereg.svg?branch=master)](https://travis-ci.org/tnagler/vinereg)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/tnagler/vinereg?branch=master&svg=true)](https://ci.appveyor.com/project/tnagler/vinereg)
[![Coverage status](https://codecov.io/gh/tnagler/vinereg/branch/master/graph/badge.svg)](https://codecov.io/github/tnagler/vinereg?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/vinereg)](https://cran.r-project.org/package=vinereg)

An R package for D-vine quantile regression.

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

For examples of its usage, have a look at the vignettes with

``` r
vignette("abalone-example", package = "vinereg")
vignette("bike-rental", package = "vinereg")
```
