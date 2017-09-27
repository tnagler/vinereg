#' Extend D-vine by one variable
#'
#' @importFrom rvinecopulib bicop hbicop vinecop_dist
xtnd_vine <- function(new_var, old_fit, ...) {
    d <- ncol(old_fit$vine$matrix) + 1
    n <- length(new_var)

    psobs <- list(
        direct = array(NA, dim = c(d, d, n)),
        indirect = array(NA, dim = c(d, d, n))
    )
    psobs$direct[-1, -d, ] <- old_fit$psobs$direct
    psobs$indirect[-1, -d, ] <- old_fit$psobs$indirect
    psobs$direct[d, d, ] <- new_var

    old_fit$vine$pair_copulas[[d - 1]] <- list()
    npars <- 0
    for (i in rev(seq.int(d - 1))) {
        zr1 <- psobs$direct[i + 1, i, ]
        zr2 <- if (i == d - 1) {
            psobs$direct[i + 1, i + 1, ]
        } else {
            psobs$indirect[i + 1, i + 1, ]
        }
        pc_fit <- bicop(cbind(zr2, zr1))
        old_fit$vine$pair_copulas[[d - i]][[i]] <- pc_fit
        npars <- npars + pc_fit$npars
        psobs$direct[i, i, ] <- hbicop(cbind(zr2, zr1), 1, pc_fit)
        psobs$indirect[i, i, ] <- hbicop(cbind(zr2, zr1), 2, pc_fit)
    }
    vine <- vinecop_dist(old_fit$vine$pair_copulas, DVineMatGen(d)[d:1, ])
    return(list(vine = vine, psobs = psobs, cll = old_fit$cll + logLik(pc_fit)))
}
