# vinereg
[![Travis-CI Build Status](https://travis-ci.org/tnagler/vinereg.svg?branch=master)](https://travis-ci.org/tnagler/vinereg)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/tnagler/vinereg?branch=master&svg=true)](https://ci.appveyor.com/project/tnagler/vinereg)
[![Coverage status](https://codecov.io/gh/tnagler/vinereg/branch/master/graph/badge.svg)](https://codecov.io/github/tnagler/vinereg?branch=master)

An R package for D-vine quantile regression.

## How to install

``` r
# install.packages("devtools")
devtools::install_github("tnagler/vinereg", build_vignettes = TRUE)
```

For an example of its usage, have a look at the vignettes with

``` r
vignette("abalone-example", package = "vinereg")
vignette("bike-rental", package = "vinereg")
```
