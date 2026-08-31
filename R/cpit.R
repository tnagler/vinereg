#' Conditional probability integral transform
#'
#' Evaluates the fitted conditional distribution of the response given the
#' covariates. For a discrete response, this is the conditional CDF at the
#' observed category, not a randomized probability integral transform.
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata a data frame containing the response and covariates from the
#'   original formula, with matching classes and factor levels.
#' @param cores integer; the number of cores to use for computations.
#'
#' @return A numeric vector containing one conditional CDF value per row of
#'   `newdata`.
#'
#' @export
#'
#' @examples
#' \dontshow{
#' set.seed(1)
#' }
#' # simulate data
#' x <- matrix(rnorm(100), 50, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))
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
#' @param newdata a data frame containing the response and covariates from the
#'   original formula, with matching classes and factor levels.
#' @param cores integer; the number of cores to use for computations.
#'
#' @return The scalar conditional log-likelihood evaluated on `newdata`.
#'
#' @export
#'
#' @examples
#' \dontshow{
#' set.seed(1)
#' }
#' # simulate data
#' x <- matrix(rnorm(100), 50, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))
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

#' Conditional density or probability mass
#'
#' Calculates the conditional density of a continuous response or conditional
#' probability mass of a discrete response given the covariates.
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata a data frame containing the response and covariates from the
#'   original formula, with matching classes and factor levels.
#' @param cores integer; the number of cores to use for computations.
#'
#' @return A numeric vector containing one conditional density or probability
#'   mass per row of `newdata`.
#'
#' @export
#'
#' @examples
#' \dontshow{
#' set.seed(1)
#' }
#' # simulate data
#' x <- matrix(rnorm(100), 50, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))
#'
#' # fit vine regression model
#' fit <- vinereg(y ~ ., dat)
#'
#' cpdf(fit, dat)
cpdf <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  dens_marg <- if (inherits(object$margins[[1]], "kde1d")) {
    kde1d::dkde1d(newdata[, 1], object$margins[[1]])
  } else {
    1
  }
  newdata <- to_uscale(newdata, object$margins)
  cond_dens_cpp(newdata, object$vine, cores) * dens_marg
}
