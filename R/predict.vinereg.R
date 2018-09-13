#' Predict conditional mean and quantiles from a D-vine regression model
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of covariate values for which to predict the quantile.
#' @param alpha vector of quantile levels; `NA` predicts the mean based on an
#'   average of the `1:10 / 11`-quantiles.
#' @param uscale if \code{TRUE} input (newdata) and output is on copula scale.
#' @param ... unused.
#'
#' @return A data.frame of quantiles where each column corresponds to one
#' value of `alpha`.
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(200), 100, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(100, 2, 0.5)))
#'
#' # fit vine regression model
#' (fit <- vinereg(y ~ ., dat))
#'
#' # inspect model
#' summary(fit)
#' plot_effects(fit)
#'
#' # model predictions
#' mu_hat  <- predict(fit, newdata = dat, alpha = NA)          # mean
#' med_hat <- predict(fit, newdata = dat, alpha = 0.5)         # median
#'
#' # observed vs predicted
#' plot(cbind(y, mu_hat))
#'
#' ## fixed variable order (no selection)
#' (fit <- vinereg(y ~ ., dat, order = c("x.2", "x.1", "z.1")))
#'
#' @seealso \code{\link{vinereg}}
#'
#' @export
#'
#' @importFrom kde1d pkde1d qkde1d
#' @importFrom stats predict
predict.vinereg <- function(object, newdata, alpha = 0.5, uscale = FALSE, ...) {
    if (missing(newdata))
        return(fitted.vinereg(object, alpha = alpha, uscale = uscale))

    stopifnot(length(alpha) > 0)
    if (any(is.na(alpha)) & inherits(object$model_frame[[1]], "ordered"))
        stop("cannot predict mean for ordered response.")

    # predict the conditional mean if alpha contains NA
    if (any(is.na(alpha))) {
        alpha <- alpha[!is.na(alpha)]  # remove NA for quantile estimation
        if (uscale)
            stop("predicting the mean is not meaningful with `uscale = TRUE.")
        preds_mean <- predict_mean(object, newdata)
    } else {
        preds_mean <- NULL
    }

    ## computation of conditional quantiles
    if (length(alpha) > 0) {
        stopifnot(is.numeric(alpha), all(alpha > 0), all(alpha < 1))

        ## preprocessing
        uscale <- adjust_uscale(object, uscale)
        newdata <- prepare_newdata(newdata, object)

        ## quantile prediction on u-scale
        if (!uscale)
            newdata <- to_uscale(newdata, object)
        preds <- qdvine(newdata, alpha, vine = object$vine)

        ## actual predictions on original scale
        if (!uscale)
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
qdvine <- function(u, alpha, vine) {
    d <- dim(vine)[1]
    if (ncol(u) != d - 1)
        stop("Dimensions of u and vine are not compatible")

    cpits <- rosenblatt(cbind(0.5, u), vine)[, -1, drop = FALSE]
    q_hat <- lapply(
        alpha,
        function(a) inverse_rosenblatt(cbind(a, cpits), vine)[, 1]
    )

    q_hat <- as.data.frame(q_hat)
    names(q_hat) <- alpha
    q_hat
}
