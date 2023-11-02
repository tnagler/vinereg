# vinereg 0.9.2

BUG FIX

* add compiler flag to prevent boost/functional from using `unary_function`.


# vinereg 0.9.1

BUG FIX

* fix unnecessary error when calling `vinereg()` with weights.
 

# vinereg 0.9.0

NEW FEATURE

* New function `cll()` to compute the conditional log-likelihood.
 
 
# vinereg 0.8.3

BUG FIX

* avoid bit-wise operations on boolean variables (fixes warnings for clang>=14).
 

# vinereg 0.8.2

DEPENDS

* require recent version of vinecopulib (>= 0.6.1.1.2) to ensure compatible
 RcppThread versions.
 

# vinereg 0.8.1

BUG FIXES

* `require()` calls with single argument in vignettes.


vinereg 0.8.0

BUG FIXES

* fix `cpit()` (last conditioning was sometimes omitted).

* prevent `rvinecopulib` from spawning own threads.

NEW FEATURES

* add `uscale` option to allow for external marginal modeling.


# vinereg 0.7.3

BUG FIXES

* fix simulated data size in documentation examples.

# vinereg 0.7.3

BUG FIXES

* properly handle case where no covariates are selected.
* conditional use of packages in Suggests.


# vinereg 0.7.2

BUG FIXES

* remove bias in quantiles for discrete variables.


# vinereg 0.7.1

This is a maintenance release following an update in rvinecopulib.

DEPENDS

* requires rvinecopulib (>= 0.5.4.1.0) to fix an unitialized value issue.

NEW FEATURES

* variables generated from `factor`s are now named with the corresponding factor level.


# vinereg 0.7.0

DEPENDENCIES

* removed dependence on future and furrr packages.

NEW FEATURES

* faster runtimes, especially for parallelized code.

* handle discrete variables properly with both parametric and  nonparametric 
  pair-copulas.

REMOVED FEATURES

* removed support for `uscale` argument.



# vinereg 0.6.0

NEW FEATURES

* new function `cpit()` to compute the conditional cdf.


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

* Use `furrr` and `future` packages instead of `parallel`, `doParallel`, and 
  `foreach` for parallelization.

NEW FEATURES

* New `print()` and `summary()` generics for `vinereg` objects.

* New `plot_effects()` method to show the marginal effects of variables.

* Allow to predict the mean with `predict(object, alpha = NA)`.


# vinereg 0.2.0

* First official release.
