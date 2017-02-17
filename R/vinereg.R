#' D-vine quantile regression
#'
#' @param y response vector
#' @param x covariate matrix
#' @param familyset either \code{"kde"} for kernel estimation of the D-vine or a
#'   vector of integers (see \code{\link{BiCopSelect}}).
#' @param correction correction criterion for the conditional log-likelihood.
#'   \code{NA} (default) imposes no correction; other choices are \code{"AIC"}
#'   and \code{"BIC"}.
#' @param par_1d list of options passed to [kdevine::kde1d()].
#' @param cores integer.
#' @param uscale logical indicating whether the data are already on copula scale
#'   (no margins have to be fitted).
#' @param ... further arguments passed to \code{\link{kde1d}},
#'   \code{\link{BiCopSelect}} or \code{\link{kdecop}}.
#'
#' @return An object of class \code{vinereg}.
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(300), 100, 3)
#' y <- x %*% c(1, -1, 2)
#'
#' # fit vine regression model (parametric)
#' fit <- vinereg(y = y, x = x, familyset = NA)
#'
#' # model predictions (median)
#' pred <- predict(fit, newdata = x, alpha = 0.5)
#'
#' # observed vs predicted
#' plot(y, pred)
#'
#' @seealso \code{\link{predict.vinereg}}
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom VineCopula BiCopSelect BiCopHfunc BiCop
#' @importFrom VineCopula RVineLogLik RVineMatrix
#' @importFrom kdevine kde1d pkde1d
#' @importFrom kdecopula kdecop hkdecop
#'
vinereg <- function(y, x, familyset = "kde", correction = NA, par_1d = list(),
                    cores = 1, uscale = FALSE, ...) {
    ## adjust input
    x <- as.matrix(x)
    d <- ncol(x)
    if (is.null(colnames(x)))
        colnames(x) <- paste0("x", 1:ncol(x))
    n <- nrow(x)
    if (!is.null(par_1d$xmin)) {
        if (length(par_1d$xmin) != d + 1)
            stop("'xmin' has to be of length d + 1")
    }
    if (!is.null(par_1d$xmax)) {
        if (length(par_1d$xmax) != d + 1)
            stop("'xmin' has to be of length d + 1")
    }
    if (!is.null(par_1d$xmax)) {
        if (length(par_1d$xmax) != d + 1)
            stop("'xmin' has to be of length d + 1")
    }
    if (length(par_1d$bw) != d && !is.null(par_1d$bw))
        stop("'bw' hast to be of length d + 1")
    if (is.null(par_1d$mult)) {
        par_1d$mult <- 1
    }
    if (length(par_1d$mult) == 1)
        par_1d$mult <- rep(par_1d$mult, d + 1)
    if (length(par_1d$mult) != d + 1)
        stop("mult.1d has to be of length 1 or d + 1")

    ## register parallel backend
    if (cores != 1 | is.na(cores)) {
        if (is.na(cores))
            cores <- max(1, parallel::detectCores() - 1)
        if (cores > 1) {
            cl <- makeCluster(cores)
            registerDoParallel(cl)
            on.exit(try(stopCluster(), silent = TRUE))
            on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
        }
    }

    ## estimation of the marginals and transformation to copula data
    dat <- udat <- cbind(y, x)
    if (uscale) {
        V <- y
        U <- x
        margs <- lapply(seq.int(d + 1), function(i) NULL)
    } else {
        est_margs <- function(k) {
            fit <- kde1d(dat[, k],
                         xmin = par_1d$xmin[k],
                         xmax = par_1d$xmax[k],
                         bw   = par_1d$bw[k],
                         mult = par_1d$mult[k])
            u <- pkde1d(dat[, k], fit)
            list(fit = fit, u = u)
        }
        if (cores > 1) {
            k <- 0  # otherwise CRAN check complains
            margs <- foreach::foreach(k = seq.int(d + 1)) %dopar% est_margs(k)
        } else {
            margs <- lapply(seq.int(d + 1), est_margs)
        }
        udat <- sapply(margs, function(x) x$u)
    }
    V <- udat[, 1]
    U <- udat[, 2:(d + 1), drop = FALSE]

    ## initialize 1-dimensional D-vine
    if (is.na(familyset) || is.numeric(familyset)) {
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
    if (is.na(familyset) || is.numeric(familyset)) {
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
                    U[, my.index[1:(i - 1)]],
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
