xtnd_vine_np <- function(newcolumn, matrix = NA, method = "TLL2", renorm.iter = 3L,
                         mult = 1, test.level = NA, trunc.level = NA,
                         treecrit = "tau", info = TRUE,
                         currentDvine, ...) {
    ## sanity checks
    if (is.na(trunc.level))
        trunc.level <- Inf
    if (is.na(test.level))
        test.level <- 1
    d <- ncol(currentDvine$RVM$matrix) + 1
    n <- length(newcolumn)

    M <- DVineMatGen(d)

    Mold <- M
    o <- diag(M)
    M <- VineCopula:::reorderRVineMatrix(M)
    MaxMat <- VineCopula:::createMaxMat(M)

    ## initialize objects
    res <- as.list(numeric(d - 1))
    for (i in 1:(d - 1))
        res[[i]] <- as.list(numeric(d - i))
    nms <- matrix("", d - 1, d - 1)
    res <- currentDvine$RVM

    if (is.null(currentDvine$RVM$info$pair.loglik)) {
        llik <- matrix(0, d, d)
        AIC <- matrix(0, d, d)
        cAIC <- matrix(0, d, d)
        BIC <- matrix(0, d, d)
        effp <- matrix(0, d, d)
    } else {
        llik <- cbind(rbind(0, currentDvine$RVM$info$pair.loglik), 0)
        cAIC <- cbind(rbind(0, currentDvine$RVM$info$pair.cAIC), 0)
        AIC <- cbind(rbind(0, currentDvine$RVM$info$pair.AIC), 0)
        BIC <- cbind(rbind(0, currentDvine$RVM$info$pair.BIC), 0)
        effp <- cbind(rbind(0, currentDvine$RVM$info$pair.effp), 0)
    }


    V <- list()
    V$direct <- array(NA, dim = c(d, d, n))
    V$indirect <- array(NA, dim = c(d, d, n))
    for (l in 1:n) {
        V$direct[, , l] <- cbind(rbind(NA, currentDvine$V$direct[, , l]), NA)
        V$indirect[, , l] <- cbind(rbind(NA, currentDvine$V$indirect[, , l]), NA)
    }
    V$direct[d, d, ] <- newcolumn

    for (k in d:2) {
        doEst <- function(i) {
            if (k > i) {
                m <- MaxMat[k, i]
                zr1 <- V$direct[k, i, ]
                zr2 <- if (m == M[k, i]) {
                    V$direct[k, (d - m + 1), ]
                } else {
                    V$indirect[k, (d - m + 1), ]
                }

                samples <- cbind(zr2, zr1)

                indep <- ifelse(test.level < 1,
                                BiCopIndTest(zr2, zr1)$p.value >= test.level,
                                FALSE)
                if (trunc.level <= (d-k))
                    indep <- TRUE

                if (indep) {
                    cfit <- list(effp = 0,
                                 likvalues = rep(1, n),
                                 loglik = 0,
                                 effp = 0,
                                 AIC = 0,
                                 cAIC = 0,
                                 BIC = 0)
                    class(cfit) <- c("kdecopula", "indep.copula")
                } else {
                    cfit <- kdecop(samples,
                                   mult = mult,
                                   method = method,
                                   renorm.iter = renorm.iter,
                                   info = info)
                }

                hfit <- list()
                direct <- hkdecop(samples,
                                  obj = cfit,
                                  cond.var = 1L)

                indirect <- hkdecop(samples,
                                    obj = cfit,
                                    cond.var = 2L)
                names <- kdevine:::naming(Mold[c(i, k:d), i])
                res.ki <- list(c = cfit, name = names)
                return(list(direct = direct,
                            indirect = indirect,
                            res.ki = res.ki))
            } else {
                return(NULL)
            }
        }

        res.k <- doEst(k - 1)

        res[[d - k + 1]][[k - 1]] <- res.k$res.ki
        V$direct[k - 1, k - 1, ] <- res.k$direct
        V$indirect[k - 1, k - 1, ] <- res.k$indirect

        cfit <- res.k$res.ki$c$info
        llik[k, k - 1] <- cfit$loglik
        effp[k, k - 1] <- cfit$effp
        AIC[k, k - 1] <- -2 * cfit$loglik + 2 * effp[k, k - 1]
        cAIC[k, k - 1] <-
            AIC[k, k - 1] + (2 * effp[k, k - 1] * (effp[k, k - 1] + 1)) /
            (n - effp[k, k - 1] - 1)
        BIC[k, k - 1] <- - 2 * cfit$loglik + log(n) * effp[k, k - 1]

        ## clean up and finalize
        # i <- k - 1
        # nums <- Mold[c(i, k:d), i]
        # #nums[1:2] <- nums[2:1]
        # name <- kdevine:::naming(nums)
    } # end k = d:2
    ## finalize results


    res$matrix <- Mold
    res$info <- list(loglik      = sum(llik),
                     pair.loglik = llik,
                     effp        = sum(effp),
                     pair.effp   = effp,
                     AIC         = sum(AIC),
                     pair.AIC    = AIC,
                     cAIC        = sum(AIC) +
                         (2 * sum(effp) * (sum(effp) + 1)) /
                         (n - sum(effp) - 1),
                     pair.cAIC   = cAIC,
                     BIC         = sum(BIC),
                     pair.BIC    = BIC)

    ## return results
    class(res) <- "kdevinecop"
    return(list(RVM = res, V = V))
}


