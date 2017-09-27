#' D-vine quantile regression
#'
#' @param formula an object of class "formula"; same as [stats::lm()].
#' @param data data frame (or object coercible by
#'   [base::as.data.frame()]) containing the variables in the model.
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
#' @seealso \code{\link{predict.vinereg}}
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom kdevine kde1d pkde1d
#' @importFrom stats model.frame
#'
vinereg <- function(formula, data, familyset = "kde", correction = NA, par_1d = list(),
                    cores = 1, uscale = FALSE, ...) {
    ## pre-processing --------
    # remove unused variables
    mf <- model.frame(formula, data)
    if (!(is.ordered(mf[[1]]) | is.numeric(mf[[1]])))
        stop("response must be numeric or ordered")
    # expand factors and add noise to discrete variable
    x <- cctools::cont_conv(mf)
    d <- ncol(x)

    ## register parallel backend
    if (cores > 1) {
        cl <- makeCluster(cores)
        registerDoParallel(cl)
        on.exit(try(stopCluster(), silent = TRUE))
        on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
    }

    ## estimation of the marginals and transformation to copula data
    margin_models <- fit_margins(x, par_1d, cores, uscale)
    u <- get_pits(x, margin_models, cores)

    ## initialization
    current_fit <- initialize_fit(u)
    status <- initialize_status(d, correction)

    ## estimation and variable selection --------
    for (i in seq.int(d - 1)) {
        # check which variable update increases the conditional log-likelihood
        # of V|U_I the most
        if (cores > 1) {
            new_fits <- foreach(k = status$remaining_vars + 1, ...) %dopar%
                xtnd_vine(u[, k], current_fit, ...)
        } else {
            new_fits <- lapply(
                status$remaining_vars + 1,
                function(k) xtnd_vine(u[, k], current_fit, ...)
            )
        }
        status <- update_status(status, new_fits)
        if (status$optimum_found)
            break
        current_fit <- new_fits[[status$best_ind]]
    }

    ## adjust model matrix and names
    reorder <- status$selected_vars
    reorder[order(reorder)] <- seq_along(status$selected_vars)
    current_fit$vine$matrix <- DVineMatGen(elements = c(1, reorder + 1))

    ## return results
    out <- list(margins = margin_models,
                vine = current_fit$vine,
                order = colnames(x[, -1])[status$selected_vars],
                selected_vars = status$selected_vars,
                formula = formula,
                model_frame = mf)
    class(out) <- "vinereg"
    out
}


fit_margins <- function(x, par_1d, cores, uscale) {
    d <- ncol(x)
    par_1d <- process_par_1d(par_1d, d)
    if (uscale) {
        # data are uniform, no need to estimate margins
        margs <- lapply(seq.int(d), function(i) NULL)
    } else {
        fit_margin <- function(k) {
            kde1d(x[, k],
                  xmin = par_1d$xmin[k],
                  xmax = par_1d$xmax[k],
                  bw   = par_1d$bw[k],
                  mult = par_1d$mult[k])
        }
        if (cores > 1) {
            margs <- foreach::foreach(k = seq.int(d)) %dopar% fit_margin(k)
        } else {
            margs <- lapply(seq.int(d), fit_margin)
        }
    }
}

get_pits <- function(x, margin_models, cores) {
    if (is.null(margin_models[[1]])) {
        # data are uniform, no need for PIT
        u <- x
    } else {
        get_pit <- function(m) pkde1d(m$x_cc, m)
        if (cores > 1) {
            u <- foreach::foreach(m = margin_models) %dopar% get_pit(m)
        } else {
            u <- lapply(margin_models, get_pit)
        }
    }

    do.call(cbind, u)
}


process_par_1d <- function(pars, d) {
    if (!is.null(pars$xmin)) {
        if (length(pars$xmin) != d)
            stop("'xmin'  must be a vector with one value for each variable")
    }
    if (!is.null(pars$xmax)) {
        if (length(pars$xmax) != d)
            stop("'xmin'  must be a vector with one value for each variable")
    }
    if (!is.null(pars$xmax)) {
        if (length(pars$xmax) != d)
            stop("'xmin' must be a vector with one value for each variable")
    }
    if (length(pars$bw) != d && !is.null(pars$bw))
        stop("'bw' must be a vector with one value for each variable")
    if (is.null(pars$mult)) {
        pars$mult <- 1
    }
    if (length(pars$mult) == 1)
        pars$mult <- rep(pars$mult, d)
    if (length(pars$mult) != d)
        stop("mult.1d has to be of length 1 or the number of variables")

    pars
}

# u_k <- pkde1d(x[, k], fit)
# list(fit = fit, u = pkde1d(x[, k], fit))
#
initialize_fit <- function(u) {
    list(
        # 1-dimensional (= empty) vine
        vine = list(pair_copulas = list(list()), matrix = as.matrix(1)),
        # array for storing pseudo-observations
        psobs = list(
            direct = array(u[, 1], dim = c(1, 1, nrow(u))),
            indirect = array(NA, dim = c(1, 1, nrow(u)))
        ),
        # conditional log-likelihood of the model
        cll = 0
    )
}

initialize_status <- function(d, correction) {
    list(
        # remaining variable indices to select from
        remaining_vars = seq.int(d - 1),
        # variables indices included in the model
        selected_vars = NULL,
        # which correction should be used for the selection criterion
        correction = correction,
        # current fit criterion (cll + correction)
        current_crit = -Inf,
        # TRUE when no improvement is possible
        optimum_found = FALSE
    )
}

update_status <- function(status, new_vines) {
    crits <- sapply(new_vines, calculate_crit, status$correction)
    if (max(crits) <= status$current_crit) {
        # optimum found, keep old fit
        status$optimum_found <- TRUE
    } else {
        status$best_ind <- which.max(crits)
        status$selected_vars <- c(status$selected_vars,
                                  status$remaining_vars[status$best_ind])
        status$current_crit <- max(crits)
        status$remaining_vars <- setdiff(status$remaining_vars, status$selected_vars)
    }

    status
}

calculate_crit <- function(fit, correction) {
    crit <- fit$cll
    if (!is.na(correction)) {
        crit <- crit - switch(correction,
                              "AIC" = fit$vine$npars,
                              "BIC" = fit$vine$npars * log(dim(fit$psobs)[3]) / 2)
    }

    crit
}
