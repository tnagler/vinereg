#' Extend D-vine by one variable
#'
#' @importFrom rvinecopulib bicop hbicop vinecop_dist
xtnd_vine <- function(new_var, currentDvine, ...) {
    d <- ncol(currentDvine$RVM$matrix) + 1
    n <- length(new_var)

    V <- list(
        direct = array(NA, dim = c(d, d, n)),
        indirect = array(NA, dim = c(d, d, n))
    )
    V$direct[-1, -d, ] <- currentDvine$V$direct
    V$indirect[-1, -d, ] <- currentDvine$V$indirect
    V$direct[d, d, ] <- new_var

    currentDvine$RVM$pair_copulas[[d - 1]] <- list()
    for (i in rev(seq.int(d - 1))) {
        zr1 <- V$direct[i + 1, i, ]
        zr2 <- if (i == d - 1) {
            V$direct[i + 1, i + 1, ]
        } else {
            V$indirect[i + 1, i + 1, ]
        }
        pc_fit <- bicop(cbind(zr2, zr1))
        currentDvine$RVM$pair_copulas[[d - i]][[i]] <- pc_fit
        V$direct[i, i, ] <- hbicop(cbind(zr2, zr1), 1, pc_fit)
        V$indirect[i, i, ] <- hbicop(cbind(zr2, zr1), 2, pc_fit)
    }

    RVM <- vinecop_dist(currentDvine$RVM$pair_copulas, DVineMatGen(d)[d:1, ])
    return(list(RVM = RVM, V = V, cll = currentDvine$cll + logLik(pc_fit)))
}
