#' Predict quantiles from a D-vine regression
#'
#' @param object an object of class \code{vinereg}.
#' @param newdata matrix of covariate values for which to predict the quantile.
#' @param alpha vector of quantile levels.
#' @param uscale if \code{TRUE} input (newdata) and output is on copula scale.
#'
#' @return A vector of quantiles if alpha is a single number; a matrix of
#' quantiles with \code{length(alpha)} columns otherwise.
#'
#' @export
#' @importFrom VineCopula RVineSim BiCopHfunc
#' @importFrom kdevine rkdevinecop qkde1d
#' @importFrom kdecopula hkdecop
#'
predict.vinereg <- function(object, newdata, alpha = 0.5, uscale = F) {
    d <- length(object$margins) - 1
    X <- matrix(newdata, ncol = d)
    n <- nrow(X)

    if (uscale == T) {
        U <- X
    } else {
        U <- matrix(0, n, d)
        for (j in object$used) {
            U[, j] <- pkde1d(X[, j], object$margins[[j + 1]])
        }
    }
    U <- U[, object$used, drop = F]
    DVM <- object$DVM
    if (object$copula.type == "kernel") {
        my.order <- diag(DVM$matrix)[-1] - 1
        U <- U[, my.order, drop = F]
        DVM$matrix <- DVineMatGen(ncol(DVM$matrix))
    } else {
        my.order <- diag(DVM$Matrix)[-1] - 1
        U <- U[, my.order, drop = F]
        DVM$Matrix <- DVineMatGen(ncol(DVM$Matrix))
    }

    uq <- get_uq(U, alpha, DVM = DVM, object$copula.type)
    if (!uscale) {
        if (ncol(uq) > 1) {
            q <- apply(uq, 2, qkde1d, obj = object$margins[[1]])
        } else {
            q <- qkde1d(uq, object$margins[[1]])
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
