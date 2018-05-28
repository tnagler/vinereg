#' Predict quantiles from a D-vine regression
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of covariate values for which to predict the quantile.
#' @param alpha vector of quantile levels.
#' @param uscale if \code{TRUE} input (newdata) and output is on copula scale.
#' @param ... unused.
#'
#' @return A data.frame of quantiles where each column corresponds to one
#' value of `alpha`.
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(300), 100, 3)
#' y <- x %*% c(1, -1, 2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(100, 3, 0.5)))
#'
#' # fit vine regression model (parametric)
#' fit <- vinereg(y ~ ., dat, family_set = "parametric")
#'
#' # selected variables and their order in the D-vine
#' fit$order
#'
#' # model predictions (median)
#' pred <- predict(fit, newdata = dat, alpha = 0.5)
#'
#' # since we evaluate at the training data, this is equivalent to
#' pred <- fitted(fit, alpha = 0.5)
#'
#' # observed vs predicted
#' plot(cbind(y, pred))
#'
#'
#' ## fixed variable order (no selection)
#' fit <- vinereg(y ~ ., dat, order = c("x.3", "x.1", "x.2", "z.1"))
#' fit$order
#'
#' @seealso \code{\link{vinereg}}
#'
#' @export
#'
#' @importFrom kde1d pkde1d qkde1d
#' @importFrom rvinecopulib rvinecop
#'
predict.vinereg <- function(object, newdata, alpha = 0.5, uscale = FALSE, ...) {
    # check if margins were estimated
    if (is.null(object$margins[[1]]) & (!uscale)) {
        warning("no margins have been estimated, setting uscale = TRUE")
        uscale <- TRUE
    }

    # check if all variables in the model are in newdata
    missing_vars <- setdiff(colnames(object$model_frame)[-1], colnames(newdata))
    if (length(missing_vars) > 0)
        stop("'newdata' is missing variables '", paste(missing_vars, sep = "', '"), "'")

    # expand factors and make ordered variables numeric
    x <- cctools::expand_as_numeric(newdata[, colnames(object$model_frame)[-1], drop = FALSE])

    # remove unused variables
    selected_vars <- match(object$order, colnames(x))
    u <- x <- x[, selected_vars, drop = FALSE]

    # transform to uniform scale
    if (!uscale) {
        for (j in seq_len(ncol(x)))
            u[, j] <- pkde1d(x[, j], object$margins[[selected_vars[j] + 1]])
    }

    # calculate quantile on uniform scale
    q <- qdvine(u, alpha, vine = object$vine)

    # transform to original scale
    if (!uscale) {
        q <- lapply(q, qkde1d, obj = object$margins[[1]])
        # when response is discrete, we need to adjust the quantiles accordingly
        if (inherits(object$model_frame[[1]], "ordered"))
            q <- lapply(q, with_levels, lvls = levels(object$model_frame[[1]]))
    }

    ## always return as data frame
    q <- as.data.frame(q)
    names(q) <- alpha
    q
}

#' @rdname predict.vinereg
#' @export
fitted.vinereg <- function(object, alpha = 0.5, ...) {
    predict.vinereg(object, newdata = object[["model_frame"]], alpha = alpha)
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

