#' Predict quantiles from a D-vine regression
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of covariate values for which to predict the quantile.
#' @param alpha vector of quantile levels.
#' @param uscale if \code{TRUE} input (newdata) and output is on copula scale.
#' @param ... unused.
#'
#' @return A vector of quantiles if alpha is a single number; a matrix of
#' quantiles with \code{length(alpha)} columns otherwise.
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(300), 100, 3)
#' y <- x %*% c(1, -1, 2)
#' dat <- data.frame(y = y, x = x)
#'
#' # fit vine regression model (parametric)
#' fit <- vinereg(y ~ ., dat, familyset = NA)
#'
#' # model predictions (median)
#' pred <- predict(fit, newdata = dat, alpha = 0.5)
#'
#' # observed vs predicted
#' plot(y, pred)
#'
#' @seealso \code{\link{vinereg}}
#'
#' @export
#'
#' @importFrom kdevine pkde1d qkde1d
#' @importFrom rvinecopulib rvinecop
#'
predict.vinereg <- function(object, newdata, alpha = 0.5, uscale = FALSE, ...) {
    if (missing(newdata))
        newdata <- object$data
    if (is.null(object$margins[[1]]) & (!uscale)) {
        warning("no margins have been estimated, setting uscale = TRUE")
        uscale <- TRUE
    }
    # remove unused variables
    missing_vars <- setdiff(colnames(object$model_frame)[-1], colnames(newdata))
    if (length(missing_vars) > 0)
        stop("'newdata' is missing variables '", paste(missing_vars, sep = "', '"), "'")
    x <- cctools::expand_as_numeric(newdata[, colnames(object$model_frame)[-1]])
    u <- x <- x[, object$selected_vars, drop = FALSE]
    if (!uscale) {
        for (j in seq.int(ncol(x))) {
            u[, j] <- pkde1d(x[, j], object$margins[[object$selected_vars[j] + 1]])
        }
    }

    uq <- qdvine(u, alpha, vine = object$vine)
    if (!uscale) {
        if (NCOL(uq) > 1) {
            q <- apply(uq, 2, qkde1d, obj = object$margins[[1]])
        } else {
            q <- qkde1d(uq, object$margins[[1]])
        }
        if (inherits(object$model_frame[[1]], "ordered")) {
            lvls <- levels(object$model_frame[[1]])
            q <- ceiling(q)
            q <- pmax(q, 1)
            q <- pmin(q, length(lvls))
            q <- ordered(lvls[q], levels = lvls)
        }
    } else {
        q <- uq
    }

    q
}

qdvine <- function(u, alpha, vine) {
    d <- ncol(vine$matrix)
    if (ncol(u) != d - 1)
        stop("Dimensions of u and vine are not compatible")
    M <- DVineMatGen(d)
    vine$matrix <- M[d:1, ]

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
        tmp <- u
    }

    ## predict quantile
    uq <- sapply(alpha,
                 function(a)
                     matrix(rvinecop(n, vine, U = cbind(a, tmp)), ncol = d)[, 1])
    if (length(alpha) > 1)
        colnames(uq) <- alpha
    uq
}

