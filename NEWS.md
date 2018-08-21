# vinereg 0.5.0

DEPENDS

* require rvinecopulib (>= 0.3.0) due to breaking changes in this package.

BUG FIXES

* prevent nan errors in loglik calculation.

* allow for empty and bivariate models.

* properly pass degree parameter for margin al estimation.
  
  
# vinereg 0.4.0

BUG FIXES

* Fix handling of `uscale` in `fitted.vinereg()`.

* Fix handling of `mult` parameter for pair-copula fits in `vinereg()`.

* Fix orientation of asymmetric pair-copulas.
  

# vinereg 0.3.0

DEPENDS

* Use `furrr` and `fututre` packages instead of `parallel`, `doParallel`, and 
  `foreach` for parallelization.

NEW FEATURES

* New `print()` and `summary()` generics for `vinereg` objects.

* New `plot_effects()` method to show the marginal effects of variables.

* Allow to predict the mean with `predict(object, alpha = NA)`.


# vinereg 0.2.0

* First official release.
