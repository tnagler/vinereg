xtnd_vine_p <- function(newcolumn, familyset = NA, selectioncrit = "BIC",
                        indeptest = TRUE, level = 0.05, trunclevel = NA, rotations = TRUE,
                        currentDvine, ...) {
    d <- ncol(currentDvine$RVM$Matrix) + 1
    n <- length(newcolumn)

    if (is.na(trunclevel))
        trunclevel <- n
    types <- familyset
    if (trunclevel == 0)
        types <- 0

    Types <- cbind(rbind(0,currentDvine$RVM$family),0)
    Params <- cbind(rbind(0,currentDvine$RVM$par),0)
    Params2 <- cbind(rbind(0,currentDvine$RVM$par2),0)

    V <- list(
        direct = array(NA, dim = c(d, d, n)),
        indirect = array(NA, dim = c(d, d, n))
    )
    V$direct[-1, -d, ] <- currentDvine$V$direct
    V$indirect[-1, -d, ] <- currentDvine$V$indirect
    V$direct[d, d, ] <- newcolumn

    for (i in rev(seq.int(d - 1))) {
        zr1 <- V$direct[i + 1, i, ]
        zr2 <- if (i == d - 1) {
            V$direct[i + 1, i + 1, ]
        } else {
            V$indirect[i + 1, i + 1, ]
        }
        cfit <- BiCopSelect(zr2, zr1, familyset, selectioncrit,
                            indeptest, level, weights = NA, rotations, ...)
        myHfunc <- BiCopHfunc(zr2, zr1, cfit, check.pars = FALSE)

        Types[i + 1, i] <- cfit$family
        Params[i + 1, i] <- cfit$par
        Params2[i + 1, i] <- cfit$par2

        V$direct[i, i, ] <- myHfunc$hfunc1
        V$indirect[i, i, ] <- myHfunc$hfunc2
    }
    RVM <- RVineMatrix(DVineMatGen(d),
                       family = Types,
                       par = Params,
                       par2 = Params2,
                       names = paste0("V", 1:d))
    return(list(RVM = RVM, V = V))
}
