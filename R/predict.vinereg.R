#' Predict conditional mean and quantiles from a D-vine regression model
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of covariate values for which to predict the quantile.
#' @param alpha vector of quantile levels; `NA` predicts the mean based on an
#'   average of the `1:20 / 21`-quantiles.
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
#' @importFrom rvinecopulib rvinecop
#' @importFrom stats predict
predict.vinereg <- function(object, newdata, alpha = 0.5, uscale = FALSE, ...) {
    if (missing(newdata))
        return(fitted.vinereg(object, alpha, uscale))
    stopifnot(length(alpha) > 0)
    if (any(is.na(alpha)) & inherits(object$model_frame[[1]], "ordered"))
        stop("cannot predict mean for ordered response.")

    # check if margins were estimated
    if (!is.null(object$margins[[1]]$u) & (!uscale)) {
        warning("no margins have been estimated, setting uscale = TRUE")
        uscale <- TRUE
    }

    # check if all variables in the model are in newdata
    if (is.matrix(newdata))
        newdata <- as.data.frame(newdata)
    missing_vars <- setdiff(colnames(object$model_frame)[-1], colnames(newdata))
    if (length(missing_vars) > 0)
        stop("'newdata' is missing variables ",
             paste(missing_vars, collapse = ", "))

    # predict the mean if alpha contains NA
    if (any(is.na(alpha))) {
        alpha <- alpha[!is.na(alpha)]
        preds <- predict_mean(object, newdata, uscale)
    } else {
        preds <- NULL
    }

    if (length(alpha) > 0) {
        stopifnot(is.numeric(alpha))
        stopifnot(all(alpha > 0) & all(alpha < 1))

        # expand factors and make ordered variables numeric
        x <- cctools::expand_as_numeric(newdata[colnames(object$model_frame)[-1]])

        # remove unused variables
        selected_vars <- match(object$order, colnames(x))
        u <- x <- x[, selected_vars, drop = FALSE]

        # transform to uniform scale
        if (!uscale) {
            for (j in seq_len(ncol(x)))
                u[, j] <- pkde1d(x[, j], object$margins[[selected_vars[j] + 1]])
        }

        # calculate quantile on uniform scale
        q_hat <- qdvine(u, alpha, vine = object$vine)

        # transform to original scale
        if (!uscale) {
            q_hat <- lapply(q_hat, qkde1d, obj = object$margins[[1]])
            if (inherits(object$model_frame[[1]], "ordered")) {
                # when response is discrete, we need to adjust the quantiles
                lvls <- levels(object$model_frame[[1]])
                q_hat <- lapply(q_hat, with_levels, lvls = lvls)
            }
        }

        ## always return as data frame
        q_hat <- as.data.frame(q_hat)
        names(q_hat) <- alpha
        if (!is.null(preds)) {
            preds <- cbind(preds, q_hat)
            alpha <- c("mean", alpha)
        } else {
            preds <- q_hat
        }
    } else {
        alpha <- "mean"
    }


    preds <- as.data.frame(preds)
    names(preds) <- alpha
    preds
}

#' @rdname predict.vinereg
#' @importFrom stats fitted
#' @export
fitted.vinereg <- function(object, alpha = 0.5, ...) {
    predict.vinereg(object, newdata = object$model_frame, alpha = alpha)
}


predict_mean <- function(object, newdata, uscale) {
    preds <- predict.vinereg(object, newdata, alpha = 1:20 / 21, uscale)
    data.frame(mean = rowMeans(preds))
}

with_levels <- function(q, lvls) {
    q <- ceiling(q)
    q <- pmax(q, 1)
    q <- pmin(q, length(lvls))
    q <- ordered(lvls[q], levels = lvls)
    q
}

qdvine <- function(u, alpha, vine) {
    d <- ncol(vine$matrix)
    if (ncol(u) != d - 1)
        stop("Dimensions of u and vine are not compatible")
    vine$matrix <- gen_dvine_mat(d)

    ## obtain diagonal entries in V matrix
    n <- nrow(u)
    V <- array(NA, dim = c(d, d, n))
    V[d, -1, ] <- t(u)
    V2 <- V
    if (d > 2) {
        for (j in (d - 1):2) {
            for (k in (d - 1):j) {
                tmp <- cbind(V2[k + 1, j, ], V[k + 1, j + 1, ])
                V[k, j, ]  <- hbicop(tmp, 1, vine$pair_copulas[[d - k]][[j]])
                V2[k, j, ] <- hbicop(tmp, 2, vine$pair_copulas[[d - k]][[j]])
            }
        }
        tmp <- t(apply(V2, 3, diag)[-1, ])
    } else {
        tmp <- V[d, 2, ]
    }

    # return as list (will be processed further)
    lapply(alpha,
           function(a)
               matrix(rvinecop(n, vine, U = cbind(a, tmp)), ncol = d)[, 1])
}

