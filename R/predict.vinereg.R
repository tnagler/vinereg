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
#' @importFrom VineCopula RVineSim BiCopHfunc
#' @importFrom kdevine rkdevinecop qkde1d
#' @importFrom kdecopula hkdecop
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
    u <- x
    if (!uscale) {
        for (j in object$used) {
            u[, j] <- pkde1d(x[, j], object$margins[[j + 1]])
        }
    }
    u <- u[, object$order, drop = FALSE]
    if (object$copula.type == "kernel") {
        object$DVM$matrix <- DVineMatGen(ncol(object$DVM$matrix))
    } else {
        object$DVM$Matrix <- DVineMatGen(ncol(object$DVM$Matrix))
    }

    uq <- get_uq(u, alpha, DVM = object$DVM, object$copula.type)
    if (!uscale) {
        if (NCOL(uq) > 1) {
            q <- apply(uq, 2, qkde1d, obj = object$margins[[1]])
        } else {
            q <- qkde1d(uq, object$margins[[1]])
        }
        if (inherits(object$model_frame[[1]], "ordered")) {
            lvls <- levels(object$model_frame[[1]])
            q <- round(q)
            q <- pmax(q, 1)
            q <- pmin(q, length(lvls))
            q <- ordered(lvls[q], levels = lvls)
        }
    } else {
        q <- uq
    }

    q
}

get_uq <- function(u, alpha, DVM, copula.type) {
    qDvine <- switch(copula.type, "kernel" = qDvine_np, "parametric" = qDvine_p)
    qDvine(u, alpha, DVM)
}

qDvine_np <- function(u, alpha, DVM) {
    d <- ncol(DVM$matrix)
    if (ncol(u) != d - 1)
        stop("Dimensions of u and DVM are not compatible")
    DVM$matrix <- M <- DVineMatGen(d)
    for (i in (d - 1):1) {
        for (j in 1:(d - i)) {
            DVM[[i]][[j]]$name <- kdevine:::naming(M[c(j , (d - i + 1):d), j])
        }
    }

    ## obtain diagonal entries in V matrix
    n <- nrow(u)
    V <- array(NA, dim = c(d, d, n))
    V[d, -1, ] <- t(u)
    V2 <- V

    if (d > 2) {
        for (j in (d-1):2) {
            for (k in (d-1):j) {
                # temp = BiCopHfunc(V2[k+1,j],V[k+1,j+1],Fam[k+1,j],Par[k+1,j],Par2[k+1,j])
                cfit <- DVM[[d - k]][[j]]$c
                cfit$flip <- TRUE
                V2[k, j, ] <- hkdecop(cbind(V2[k + 1, j, ], V[k + 1, j + 1, ]),
                                      cfit, cond.var = 2)
                V[k, j, ]  <- hkdecop(cbind(V2[k + 1, j, ], V[k + 1, j + 1, ]),
                                      cfit, cond.var = 1)
            }
        }
        tmp <- t(apply(V2, 3, diag)[-1, ])
    } else {
        tmp <- u
    }

    ## predict quantile
    uq <- sapply(alpha,
                 function(a) rkdevinecop(n, DVM, U = cbind(a, tmp))[, 1])
    if (length(alpha) > 1)
        colnames(uq) <- alpha
    uq
}

qDvine_p <- function(u, alpha, DVM) {
    d <- ncol(DVM$Matrix)
    if (ncol(u) != d - 1)
        stop("Dimensions of u and DVM are not compatible")
    M <- DVineMatGen(d)
    DVM$Matrix <- M

    ## obtain diagonal entries in V matrix
    n <- nrow(u)
    V <- array(NA, dim = c(d, d, n))
    V[d, -1, ] <- t(u)
    V2 <- V
    if (d > 2) {
        for (j in (d-1):2) {
            for (k in (d-1):j) {
                # temp = BiCopHfunc(V2[k+1,j],V[k+1,j+1],Fam[k+1,j],Par[k+1,j],Par2[k+1,j])
                tmp <- BiCopHfunc(V2[k + 1, j, ], V[k + 1, j + 1, ],
                                  DVM$family[k + 1, j],
                                  DVM$par[k + 1, j],
                                  DVM$par2[k + 1, j],
                                  check.pars = FALSE)
                V2[k, j, ] <- tmp$hfunc2
                V[k, j, ]  <- tmp$hfunc1
            }
        }
        tmp <- t(apply(V2, 3, diag)[-1, ])
    } else {
        tmp <- u
    }

    ## predict quantile
    uq <- sapply(alpha,
                 function(a)
                     matrix(RVineSim(n, DVM, U = cbind(a, tmp)), ncol = d)[, 1])
    if (length(alpha) > 1)
        colnames(uq) <- alpha
    uq
}
