#' Conditional probability integral transform
#'
#' Calculates the conditional distribution of the response given the covariates.
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the conditional distribution.
#' @param cores integer; the number of cores to use for computations.
#'
#' @export
#'
#' @examples
#' \dontshow{
#' set.seed(1)
#' }
#' # simulate data
#' x <- matrix(rnorm(500), 250, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(250, 2, 0.5)))
#'
#' # fit vine regression model
#' fit <- vinereg(y ~ ., dat)
#'
#' hist(cpit(fit, dat)) # should be approximately uniform
cpit <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  newdata <- to_uscale(newdata, object$margins)
  cond_dist_cpp(newdata, object$vine, cores)
}

#' Conditional log-likelihood
#'
#' Calculates the conditional log-likelihood of the response given the covariates.
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the conditional distribution.
#' @param cores integer; the number of cores to use for computations.
#'
#' @export
#'
#' @examples
#' \dontshow{
#' set.seed(1)
#' }
#' # simulate data
#' x <- matrix(rnorm(500), 250, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(250, 2, 0.5)))
#'
#' # fit vine regression model
#' fit <- vinereg(y ~ ., dat)
#'
#' cll(fit, dat)
#' fit$stats$cll
cll <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  ll_marg <- if (inherits(object$margins[[1]], "kde1d")) {
    sum(log(kde1d::dkde1d(newdata[, 1], object$margins[[1]])))
  } else {
    0
  }
  newdata <- to_uscale(newdata, object$margins)
  ll_cop <- cond_loglik_cpp(newdata, object$vine, cores)
  ll_cop + ll_marg
}
