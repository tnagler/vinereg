xtnd_vine_p <- function (newcolumn, familyset = NA, selectioncrit = "BIC",
                        indeptest = TRUE, level = 0.05, trunclevel = NA, rotations = TRUE,
                        currentDvine, ...) {
    d <- n <- ncol(currentDvine$RVM$Matrix) + 1
    N <- length(newcolumn)
    if (is.na(trunclevel))
        trunclevel <- n
    types <- familyset
    if (trunclevel == 0)
        types <- 0
    Matrix = DVineMatGen(d)
    M <- Matrix
    Mold <- M
    o <- diag(M)
    M <- VineCopula:::reorderRVineMatrix(M)
    MaxMat <- VineCopula:::createMaxMat(M)
    #CondDistr <- neededCondDistr(M)
    Types <- cbind(rbind(0,currentDvine$RVM$family),0)
    Params <- cbind(rbind(0,currentDvine$RVM$par),0)
    Params2 <- cbind(rbind(0,currentDvine$RVM$par2),0)
    V <- list()
    V$direct <- array(NA, dim = c(n, n, N))
    V$indirect <- array(NA, dim = c(n, n, N))
    for (l in 1:N){
        V$direct[ , ,l] <- cbind(rbind(NA,currentDvine$V$direct[ , ,l]),NA)
        V$indirect[ , ,l] <- cbind(rbind(NA,currentDvine$V$indirect[ , ,l]),NA)
    }
    V$direct[n, n, ] <- newcolumn
    doEst <- function(i) {
        m <- MaxMat[k, i]
        zr1 <- V$direct[k, i, ]
        zr2 <- if (m == M[k, i]) {
            V$direct[k, (d - m + 1), ]
        }else {
            V$indirect[k, (d - m + 1), ]
        }
        cfit <- BiCopSelect(zr2, zr1, familyset, selectioncrit,
                            indeptest, level, weights = NA, rotations, ...)
        direct <- indirect <- NULL
        myHfunc = BiCopHfunc(zr2, zr1, cfit, check.pars = FALSE)
        #    if (CondDistr$direct[k - 1, i])
        direct <- myHfunc$hfunc1
        #    if (CondDistr$indirect[k - 1, i])
        indirect <- myHfunc$hfunc2
        list(direct = direct, indirect = indirect, cfit = cfit)
    }
    for (k in d:2) {
        i <- k-1
        res <- doEst(i)
        Types[k, i] <- res$cfit$family
        Params[k, i] <- res$cfit$par
        Params2[k, i] <- res$cfit$par2
        if (!is.null(res$direct))
            V$direct[k - 1, i, ] <- res$direct
        if (!is.null(res$indirect))
            V$indirect[k - 1, i, ] <- res$indirect
    }
    RVM <- RVineMatrix(Mold, family = Types, par = Params, par2 = Params2,
                       names = paste("V", 1:n, sep = ""))
    return(list(RVM=RVM, V=V))
}
