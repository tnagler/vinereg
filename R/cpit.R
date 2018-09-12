#' Conditional probability integral transform
#'
#' Calculates the conditional distribution of the response given the covariates.
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the conditional distribution.
#' @param uscale if \code{TRUE} input (newdata) and output is on copula scale.
#'
#' @export
#'
#' @examples
#' \dontshow{set.seed(5)}
#' # simulate data
#' x <- matrix(rnorm(200), 100, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(100, 2, 0.5)))
#'
#' # fit vine regression model
#' fit <- vinereg(y ~ ., dat)
#'
#' hist(cpit(fit, dat))  # should be approximately uniform
cpit <- function(object, newdata, uscale = FALSE) {
    # check if margins were estimated
    if (!is.null(object$margins[[1]]$u) & (!uscale)) {
        warning("no margins have been estimated, setting uscale = TRUE")
        uscale <- TRUE
    }

    ## preprocessing
    uscale <- adjust_uscale(object, uscale)
    newdata <- check_and_sort_newdata(newdata, object)
    x <- cctools::expand_as_numeric(newdata)

    rosenblatt(to_uscale(x, object), object$vine)[, 1]
}
