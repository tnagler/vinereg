#' Predict conditional mean and quantiles from a D-vine regression model
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of covariate values for which to predict the quantile.
#' @param alpha vector of quantile levels; `NA` predicts the mean based on an
#'   average of the `1:10 / 11`-quantiles.
#' @param cores integer; the number of cores to use for computations.
#' @param ... unused.
#'
#' @return A data.frame of quantiles where each column corresponds to one
#' value of `alpha`.
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(100), 50, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))
#'
#' ## fixed variable order (no selection)
#' (fit <- vinereg(y ~ ., dat, order = c("x.2", "x.1", "z.1")))
#'
#' # model predictions
#' mu_hat <- predict(fit, newdata = dat, alpha = NA) # mean
#' med_hat <- predict(fit, newdata = dat, alpha = 0.5) # median
#'
#' # observed vs predicted
#' plot(cbind(y, mu_hat))
#'
#' @seealso \code{\link{vinereg}}
#'
#' @export
#'
#' @importFrom kde1d pkde1d qkde1d
#' @importFrom stats predict
predict.vinereg <- function(object, newdata, alpha = 0.5, cores = 1, ...) {
  if (missing(newdata)) {
    return(fitted.vinereg(object, alpha = alpha))
  }

  stopifnot(length(alpha) > 0)
  if (any(is.na(alpha)) & inherits(object$model_frame[[1]], "ordered")) {
    stop("cannot predict mean for ordered response.")
  }

  # predict the conditional mean if alpha contains NA
  if (any(is.na(alpha))) {
    alpha <- alpha[!is.na(alpha)] # remove NA for quantile estimation
    preds_mean <- predict_mean(object, newdata)
  } else {
    preds_mean <- NULL
  }

  ## computation of conditional quantiles
  if (length(alpha) > 0) {
    stopifnot(is.numeric(alpha), all(alpha > 0), all(alpha < 1))

    ## preprocessing
    newdata <- prepare_newdata(newdata, object)
    newdata <- to_uscale(newdata, object$margins[-1], add_response = TRUE)
    preds <- qdvine(newdata, alpha, vine = object$vine, cores)

    ## actual predictions on original scale
    preds <- to_yscale(preds, object)
    if (!is.null(preds_mean))
      preds <- cbind(preds_mean, preds)
  } else {
    preds <- preds_mean
  }

  preds
}

#' @rdname predict.vinereg
#' @importFrom stats fitted
#' @export
fitted.vinereg <- function(object, alpha = 0.5, ...) {
  predict.vinereg(object, newdata = object$model_frame, alpha = alpha, ...)
}


#' predicts the conditional mean as the average of quantiles.
#' @noRd
predict_mean <- function(object, newdata) {
  preds <- predict.vinereg(object, newdata, alpha = 1:10 / 11)
  data.frame(mean = rowMeans(preds))
}

#' @importFrom rvinecopulib rosenblatt inverse_rosenblatt
qdvine <- function(u, alpha, vine, cores) {
  vine$var_types[1] <- "c"
  q_hat <- as.data.frame(cond_quantile_cpp(alpha, as.matrix(u), vine, cores))
  names(q_hat) <- alpha
  q_hat
}
