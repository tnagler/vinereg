#' D-vine regressions
#'
#' @param y response vector
#' @param x covariate matrix
#' @param familyset either \code{"kde"} for kernel estimation of the D-vine or
#' a vector of integers (see \code{\link{BiCopSelect}}).
#' @param correction correction criterion for the conditional log-likelihood.
#' \code{NA} (default) imposes no correction; other choices are \code{"AIC"}
#' and \code{"BIC"}.
#' @param cores integer.
#' @param ... further arguments passed to \code{\link{kde1d}},
#' \code{\link{BiCopSelect}} or \code{\link{kdecop}}.
#'
#' @return An object of class \code{vinereg}.
#'
#' @export
#'
#' @importFrom parallel makeCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom VineCopula BiCopSelect BiCopHfunc BiCop
#' @importFrom VineCopula RVineLogLik RVineMatrix
#' @importFrom kdevine kde1d pkde1d
#' @importFrom kdecopula kdecop hkdecop
#'
vinereg <- function(y, x, familyset = "kde", correction = NA, par.1d = list(),
                    cores = 1, ...) {
    ## adjust input
    x <- as.matrix(x)
    d <- ncol(x)
    dat <- udat <- cbind(y, x)
    if (is.null(colnames(x)))
        colnames(x) <- paste0("x", 1:ncol(x))
    n <- nrow(x)
    if (!is.null(par.1d$xmin)) {
        if(length(par.1d$xmin) != d + 1)
            stop("'xmin' has to be of length d + 1")
    }
    if (!is.null(par.1d$xmax)) {
        if(length(par.1d$xmax) != d + 1)
            stop("'xmin' has to be of length d")
    }
    if (!is.null(par.1d$xmax)) {
        if(length(par.1d$xmax) != d + 1)
            stop("'xmin' has to be of length d")
    }
    if (length(par.1d$bw) != d && !is.null(par.1d$bw))
        stop("'bw' hast to be of length d + 1")
    if (is.null(par.1d$mult)) {
        par.1d$mult <- 1
    }
    if (length(par.1d$mult) == 1)
        par.1d$mult <- rep(par.1d$mult, d + 1)
    if (length(par.1d$mult) != d + 1)
        stop("mult.1d has to be of length 1 or d + 1")

    ## register parallel backend
    if (cores != 1 | is.na(cores)) {
        if (is.na(cores))
            cores <- max(1, detectCores() - 1)
        if (cores > 1) {
            cl <- makeCluster(cores)
            registerDoParallel(cl)
            on.exit(try(stopCluster(), silent = TRUE))
            on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
        }
    }

    ## estimation of the marginals and transformation to copula data
    est_margs <- function(k) {
        fit <- kde1d(dat[, k],
                     xmin = par.1d$xmin[k],
                     xmax = par.1d$xmax[k],
                     bw   = par.1d$bw[k],
                     mult = par.1d$mult[k])
        u <- pkde1d(dat[, k], fit)
        list(fit = fit, u = u)
    }
    if (cores > 1) {
        margs <- foreach(k = seq.int(d + 1)) %dopar% est_margs(k)
    } else {
        margs <- lapply(seq.int(d + 1), est_margs)
    }
    udat <- sapply(margs, function(x) x$u)
    V <- udat[, 1]
    U <- udat[, 2:(d + 1), drop = FALSE]

    ## initialize 1-dimensional D-vine
    if (is.na(familyset) || is.integer(familyset)) {
        copula.type <- "parametric"
        vine <- list(RVM = list(Matrix = matrix(1, 1, 1),
                                family = matrix(0, 1, 1),
                                par = matrix(0, 1, 1),
                                par2 = matrix(0, 1, 1)),
                     V = list(direct = array(V, dim = c(1, 1, n)),
                              indirect = array(NA, dim = c(1, 1, n))))
        update <- update_p
    } else {
        copula.type <- "kernel"
        RVM <- lapply(1:d, list)
        names(RVM)[1:d] <- vapply(1:d,
                                  function(x) paste("T", x, sep = ""), "")
        RVM$matrix <- matrix(1, 1, 1)
        RVM$info <- list()
        vine <- list(RVM = RVM,
                     V = list(direct = array(V, dim = c(1, 1, n)),
                              indirect = array(NA, dim = c(1, 1, n))))
        update <- update_np
    }

    ## estimation and variable selection
    remaining.variables <- 1:d
    my.index <- NULL
    global.max.ll <- -Inf
    for (i in 1:d) {
        # check which variable update increases the loglikelihood of the
        # conditional density f_V|U_I the most
        if (cores > 1) {
            newvines <- foreach(k = seq_along(remaining.variables), ...) %dopar%
                update(k,
                       i = i,
                       my.index = my.index,
                       U = U,
                       V = V,
                       remaining.variables = remaining.variables,
                       vine = vine,
                       correction = correction,
                       familyset = familyset,
                       ...)
        } else {
            newvines <- lapply(seq_along(remaining.variables),
                               update,
                               i = i,
                               my.index = my.index,
                               U = U,
                               V = V,
                               remaining.variables = remaining.variables,
                               vine = vine,
                               correction = correction,
                               familyset = familyset,
                               ...)
        }
        # pick the one with the highest cll. If none of the remaining variables
        # increases the overall cll break the loop
        cll <- sapply(newvines, function(x) x$cll)
        maxInd <- which.max(cll)
        if (cll[maxInd] <= global.max.ll)
            break
        my.index <- c(my.index, remaining.variables[maxInd])
        global.max.ll <- max(global.max.ll, cll[maxInd])
        vine <- newvines[[maxInd]]$newvine
        remaining.variables <- setdiff(remaining.variables, my.index[i])
    }

    ## adjust model matrix and names
    reorder <- my.index
    reorder[order(reorder)] <- 1:length(my.index)
    if (is.na(familyset) || is.integer(familyset)) {
        vine$RVM$Matrix <- DVineMatGen(elements = c(1, reorder + 1))
    } else {
        M <- DVineMatGen(elements = c(1, reorder + 1))
        vine$RVM$matrix <- M
        dd <- ncol(vine$RVM$matrix)
        if (dd < d)
            vine$RVM[dd:d] <- NULL
        for (i in (dd - 1):1) {
            for (j in 1:(dd - i)) {
                vine$RVM[[i]][[j]]$name <- kdevine:::naming(M[c(j , (dd - i + 1):dd), j])
            }
        }
    }

    ## return results
    out <- list(margins = lapply(margs, function(x) x$fit),
                DVM = vine$RVM,
                order = colnames(x)[my.index],
                my.index = my.index,
                used = sort(my.index),
                copula.type = copula.type,
                data = list(y = y, x = x))
    class(out) <- "vinereg"
    out
}

update_p <- function(j, i, my.index, U, V, remaining.variables, vine, correction, ...) {
    # update current D-Vine by adding j-th remaining variable
    newvine <- xtnd_vine_p(newcolumn = U[, remaining.variables[j]],
                           currentDvine = vine,
                           ...)
    RVM <- newvine$RVM
    # number of vine's parameters for cll calculation
    npar <- sum(RVM$par != 0) + sum(RVM$par2 != 0)
    tmpdat <- cbind(V,
                    U[, my.index[1:(i-1)]],
                    U[, remaining.variables[j]])

    cll <- sum(RVineLogLik(tmpdat, RVM)$V$value[, 1])
    if (!is.na(correction))
        cll <- cll - switch(correction,
                            "AIC" = npar,
                            "BIC" = npar * log(length(V)) / 2)
    list(newvine = newvine, cll = cll)
}

update_np <- function(j, i, my.index, U, V, remaining.variables, vine, correction, ...) {
    # update current D-Vine by adding j-th remaining variable
    newvine <- xtnd_vine_np(newcolumn = U[, remaining.variables[j]],
                            currentDvine = vine,
                            ...)
    RVM <- newvine$RVM
    # number of vine's parameters for cll calculation
    npar <- sum(RVM$info$pair.effp)
    cll <- sum(RVM$info$pair.loglik[, 1])
    if (!is.na(correction))
        cll <- cll - switch(correction,
                            "AIC" = npar,
                            "BIC" = npar * log(length(V)) / 2)
    list(newvine = newvine, cll = cll)
}
