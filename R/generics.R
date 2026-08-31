#' Standard methods for D-vine regression models
#'
#' These methods extract the fitted model's call information, data, effective
#' degrees of freedom, conditional log-likelihood, and variable-wise summary.
#'
#' @param x,object,formula a `vinereg` object.
#' @param use.fallback unused; included for compatibility with [nobs()].
#' @param ... unused.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns a data frame with
#'   one row for the response and each selected predictor. `logLik()` returns an
#'   object of class `logLik`; its `df` attribute contains the effective degrees
#'   of freedom. `nobs()` returns the number of observations used for fitting.
#'   `formula()` and `model.frame()` return the model formula and frame.
#'
#' @name vinereg-methods
#' @importFrom stats formula logLik model.frame nobs
NULL

#' @rdname vinereg-methods
#' @export
print.vinereg <- function(x, ...) {
  cat("D-vine regression model: ")
  n_predictors <- length(x$order)
  if (n_predictors <= 10) {
    predictors <- paste(x$order, collapse = ", ")
  } else {
    predictors <- paste(x$order[1:10], collapse = ", ")
    predictors <- paste0(predictors, ", ... (", n_predictors - 10, " more)")
  }
  cat(names(x$model_frame)[1], "|", predictors, "\n")
  stats <- unlist(x$stats[1:5])
  stats <- paste(names(stats), round(stats, 2), sep = " = ")
  cat(paste(stats, collapse = ", "), "\n")
  invisible(x)
}

#' @rdname vinereg-methods
#' @export
summary.vinereg <- function(object, ...) {
  data.frame(
    var = c(names(object$model_frame)[1], object$order),
    edf = object$stats$var_edf,
    cll = object$stats$var_cll,
    caic = object$stats$var_caic,
    cbic = object$stats$var_cbic,
    p_value = object$stats$var_p_value
  )
}

#' @rdname vinereg-methods
#' @export
logLik.vinereg <- function(object, ...) {
  structure(
    object$stats$cll,
    df = object$stats$edf,
    nobs = object$stats$nobs,
    class = "logLik"
  )
}

#' @rdname vinereg-methods
#' @export
nobs.vinereg <- function(object, use.fallback = TRUE, ...) {
  object$stats$nobs
}

#' @rdname vinereg-methods
#' @export
formula.vinereg <- function(x, ...) {
  x$formula
}

#' @rdname vinereg-methods
#' @export
model.frame.vinereg <- function(formula, ...) {
  formula$model_frame
}

#' Plot marginal effects of a D-vine regression model
#'
#' The points show fitted conditional quantiles against each requested variable.
#' For variable \eqn{X_k}, the smooth curve estimates
#' \eqn{E[\hat Q_\alpha(Y \mid X) \mid X_k = x]}. It therefore averages over
#' the conditional distribution of the other variables. A curve for an
#' unselected variable can vary when that variable is associated with selected
#' predictors. The curve is descriptive and is not a partial-dependence or
#' causal effect.
#'
#' @param object a `vinereg` object
#'
#' @param alpha vector of quantile levels.
#' @param vars vector of expanded variable names to display.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @export
#' @examples
#' # simulate data
#' x <- matrix(rnorm(100), 50, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(50, 2, 0.5)))
#'
#' # fit vine regression model
#' fit <- vinereg(y ~ ., dat)
#' plot_effects(fit)
plot_effects <- function(object,
                         alpha = c(0.1, 0.5, 0.9),
                         vars = object$order) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package must be installed to plot.")
  }

  mf <- expand_factors(object$model_frame)
  if (!all(vars %in% colnames(mf)[-1])) {
    stop(
      "unknown variable in 'vars'; allowed values: ",
      paste(colnames(mf)[-1], collapse = ", ")
    )
  }

  preds <- fitted(object, alpha)
  preds <- lapply(seq_along(alpha), function(a)
    cbind(data.frame(alpha = alpha[a], prediction = preds[[a]])))
  preds <- do.call(rbind, preds)

  df <- lapply(vars, function(var)
    cbind(data.frame(var = var, value = as.numeric(unname(mf[, var])), preds)))
  df <- do.call(rbind, df)
  df$value <- as.numeric(df$value)
  df$alpha <- as.factor(df$alpha)

  value <- prediction <- NULL # for CRAN checks
  suppressWarnings(
    ggplot2::ggplot(df, ggplot2::aes(value, prediction, color = alpha)) +
      ggplot2::geom_point(alpha = 0.15) +
      ggplot2::geom_smooth(se = FALSE) +
      ggplot2::facet_wrap(~var, scale = "free_x") +
      ggplot2::ylab(quote(Q(y * "|" * x[1] * ",...," * x[p]))) +
      ggplot2::xlab(quote(x[k])) +
      ggplot2::theme(legend.position = "bottom")
  )
}
