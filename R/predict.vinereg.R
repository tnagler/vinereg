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
#' x <- matrix(rnorm(300), 100, 2)
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
        object$model_frame <- object$model_frame[-1]  # remove response (unused)
        newdata <- check_and_sort_newdata(newdata, object)

        # factors must be expanded to dummy numeric variables
        x <- cctools::expand_as_numeric(newdata)

        # x must be sorted in the order of the data used for fitting
        x <- x[, names(sort(object$selected_vars)), drop = FALSE]

        # quantile prediction on u-scale
        if (!uscale)
            x <- to_uscale(x, object)
        preds <- qdvine(x, alpha, vine = object$vine)

        # actual predictions on original scale
        if (!uscale)
            preds <- to_yscale(preds, object)

        if (!is.null(preds_mean))
            preds <- cbind(preds_mean, preds)
    } else {
        preds <- preds_mean
    }

    preds
}

#' checks if newdata has appropriate columns and sorts according to the order
#' used for fitting.
#' @noRd
check_and_sort_newdata <- function(newdata, object) {
    newdata <- as.data.frame(newdata)
    check_var_availability(newdata, names(object$model_frame))

    # the check_x functions expect variables in newdata and model_frame in
    # the same order
    newdata <- newdata[names(object$model_frame)]
    check_types(newdata, object$model_frame)
    check_levels(newdata, object$model_frame)

    newdata
}

#' checks if all *selected* covariates are in newdata.
#' @noRd
check_var_availability <- function(newdata, vars) {
    vars_avail <- match(vars, colnames(newdata))
    if (any(is.na(vars_avail))) {
        vars_missing <- paste(vars[is.na(vars_avail)], collapse = ", ")
        stop("'newdata' is missing variables ", vars_missing)
    }
}

#' checks if margins were estimated.
#' @noRd
adjust_uscale <- function(object, uscale) {
    if (!is.null(object$margins[[1]]$u) & (!uscale)) {
        warning("no margins have been estimated, setting uscale = TRUE")
        uscale <- TRUE
    }
    uscale
}

#' transforms data to uniform scale with probability integral transform.
#' @noRd
to_uscale <- function(x, object) {
    for (var in colnames(x))
        x[, var] <- truncate_u(pkde1d(x[, var], object$margins[[var]]))
    x
}

truncate_u <- function(u) {
    pmin(pmax(u, 1e-10), 1 - 1e-10)
}

#' transforms predicted response back to orginal variable scale.
#' @noRd
to_yscale <- function(u, object) {
    nms <- colnames(u)
    u <- lapply(u, qkde1d, obj = object$margins[[1]])
    if (inherits(object$model_frame[[1]], "ordered")) {
        # when response is discrete, we need to adjust the quantiles
        lvls <- levels(object$model_frame[[1]])
        u <- lapply(u, with_levels, lvls = lvls)
    }

    u <- as.data.frame(u)
    names(u) <- nms
    u
}

#' checks if variable types are equal in original data and new data.
#' @noRd
check_types <- function(actual, expected) {
    different_type <- sapply(
        seq_along(actual),
        function(i) !identical(class(actual[[i]])[1], class(expected[[i]])[1])
    )
    if (any(different_type)) {
        errors <- data.frame(
            expected = sapply(actual[different_type], function(x) class(x)[1]),
            actual = sapply(expected[different_type], function(x) class(x)[1])
        )
        errors <- paste(capture.output(print(errors)), collapse = "\n")
        stop("some columns have incorrect type:\n", errors, call. = FALSE)
    }
}

#' checks if factor levels are equal in original data and new data.
#' @noRd
check_levels <- function(actual, expected) {
    # only check factors
    actual   <- actual[sapply(actual, is.factor)]
    expected <- expected[sapply(expected, is.factor)]
    if (length(expected) == 0)
        return(TRUE)

    different_levels <- sapply(
        seq_along(actual),
        function(i) !identical(levels(actual[[i]]), levels(expected[[i]]))
    )
    if (any(different_levels)) {
        errors <- data.frame(
            expected = sapply(actual[different_levels],
                              function(x) paste(levels(x), collapse = ",")),
            actual = sapply(expected[different_levels],
                            function(x) paste(levels(x), collapse = ","))
        )
        errors <- paste(capture.output(print(errors)), collapse = "\n")
        stop("some factors have incorrect levels\n", errors, call. = FALSE)
    }

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

#' returns the quantile predictions as order variable with appropriate levels.
#' @noRd
with_levels <- function(q, lvls) {
    q <- ceiling(q)
    q <- pmax(q, 1)
    q <- pmin(q, length(lvls))
    ordered(lvls[q], levels = lvls)
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
