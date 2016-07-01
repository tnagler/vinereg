vinereg_cf <- function(y, x) {
    dat <- cbind(y, x)
    fit <- kdevine(dat)
    out <- list(fit = fit, y = y, x = x)
    class(out) <- "vinereg_cf"
    out
}

quantile.vinereg_cf <- function(x, q, newdata, cores = 1) {
    obj <- x
    newdata <- as.matrix(newdata)
    if (ncol(newdata) == 1)
        newdata <- t(newdata)

    n <- nrow(newdata)
    rho <- function(z, ...)
        z * (q - as.numeric(z < 0))
    critfun <- function(a, y, w, x, ...) {
        sum(rho(obj$y - a) * w)
    }
    opt_for_k <- function(k) {
        xpand <- t(sapply(1:length(obj$y),
                          function(i) newdata[k, , drop = FALSE]))
        print(k)
        optimize(critfun,
                 range(obj$y),
                 tol = 1e-2,
                 y = obj$y,
                 w = copula_weights(cbind(obj$y, xpand), obj$fit),
                 q = q)$minimum
    }

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

    if (cores > 1) {
        opt <- foreach (k = 1:n) %dopar% opt_for_k(k)
        opt <- unlist(opt)
    } else {
        opt <- vapply(1:n, opt_for_k, 1)
    }

    opt
}

copula_weights <- function(x, obj) {
    x <- as.matrix(x)
    n <- length(obj$marg.dens[[1]]$data)
    if (ncol(x) == 1)
        x <- t(x)
    d <- ncol(x)

    stopifnot(class(obj) == "kdevine")
    if (length(obj$marg.dens) != d)
        stop("'x' has incorrect dimension")

    ## evaluate copula density (if necessary)
    if (!is.null(obj$vine)) {
        u <- x
        # PIT to copula level
        for (i in 1:d)
            u[, i] <- pkde1d(x[, i], obj$marg.dens[[i]])
        vinevals <- dkdevinecop(u, obj = obj$vine, stable = TRUE)
    } else {
        vinevals <- rep(1, nrow(x))
    }

    ## final density estimate is product of marginals and copula density
    vinevals
}



