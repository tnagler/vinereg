#' D-vine quantile regression
#'
#' Sequential estimation of a regression D-vine for the purpose of quantile
#' prediction as described in Kraus and Czado (2017).
#'
#' @param formula an object of class "formula"; same as [lm()].
#' @param data data frame (or object coercible by
#'   [as.data.frame()]) containing the variables in the model.
#' @param family_set see `family_set` argument of [rvinecopulib::bicop()].
#' @param selcrit selection criterion based on conditional log-likelihood.
#'   \code{"loglik"} (default) imposes no correction; other choices are
#'   \code{"aic"} and \code{"bic"}.
#' @param order the order of covariates in the D-vine, provided as vector
#'   of variable names (after calling `cctools::cont_conv()` on the
#'   `model.frame()`); selected automatically if `order = NA` (default).
#' @param par_1d list of options passed to [kde1d::kde1d()], must be one value
#'   for each margin, e.g. `list(xmin = c(0, 0, -Inf))` if the response and
#'   first covariate have non-negative support.
#' @param cores integer.
#' @param uscale logical indicating whether the data are already on copula scale
#'   (no margins have to be fitted).
#' @param ... further arguments passed to [rvinecopulib::bicop()].
#'
#' @return
#' An object of class vinereg. It is a list containing the elements
#' \describe{
#' \item{margins}{list of marginal models fitted by [kde1d::kde1d()].}
#' \item{vine}{an [rvinecopulib::vinecop_dist()] object containing the fitted
#'   D-vine.}
#' \item{order}{order of the covariates chosen by the variable selection algorithm.}
#' \item{selected_vars}{indices of selected variables.}
#' \item{formula}{the model formula specified by the user.}
#' \item{model_frame}{the data used to fit the D-vine.}
#' }
#' This object can be plugged into [predict.vinereg()] to predict conditional
#' quantiles.
#'
#' @references
#'
#' Kraus and Czado (2017), D-vine copula based quantile regression,
#' Computational Statistics and Data Analysis, 110, 1-18
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(300), 100, 3)
#' y <- x %*% c(1, -1, 2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(100, 3, 0.5)))
#'
#' # fit vine regression model (parametric)
#' fit <- vinereg(y ~ ., dat, family_set = "parametric")
#'
#' # model predictions (median)
#' pred <- predict(fit, newdata = dat, alpha = 0.5)
#'
#' # observed vs predicted
#' plot(cbind(y, pred))
#'
#'
#' ## fixed variable order (no selection)
#' fit <- vinereg(y ~ ., dat, order = c("x.3", "x.1", "x.2", "z.1"))
#' fit$order
#'
#' @seealso \code{\link{predict.vinereg}}
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom kde1d kde1d pkde1d
#' @importFrom stats model.frame logLik
#' @importFrom rvinecopulib bicop dbicop hbicop vinecop_dist
#'
vinereg <- function(formula, data, family_set = "parametric", selcrit = "loglik",
                    order = NA, par_1d = list(), cores = 1, uscale = FALSE, ...) {
    ## pre-processing --------
    # remove unused variables
    if (!missing(data)) {
        mf <- model.frame(formula, data)
    } else {
        mf <- model.frame(formula, parent.frame())
    }

    if (any(sapply(mf, is.factor)) &
        is.na(pmatch(family_set, "nonparametric")) &
        is.na(pmatch(family_set, "tll"))) {
        warning('parametric copula families are misspecified ',
                'for jittered discrete variables; ',
                'use family_set = "nonparametric" to maintain consistency')
    }

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
    var_nms <- colnames(u)

    ## initialization
    current_fit <- initialize_fit(u)
    status <- initialize_status(d, selcrit)

    ## estimation --------
    if (any(is.na(order))) {
        ## automatic variable selection
        for (i in seq_len(d - 1)) {
            if (cores > 1) {
                k <- NULL  # for CRAN checks
                new_fits <- foreach(k = status$remaining_vars + 1, ...) %dopar%
                    xtnd_vine(u[, k], current_fit, family_set, selcrit, ...)
            } else {
                new_fits <- lapply(
                    status$remaining_vars + 1,
                    function(k)
                        xtnd_vine(u[, k], current_fit, family_set, selcrit, ...)
                )
            }
            status <- update_status(status, new_fits)
            if (status$optimum_found)
                break
            current_fit <- new_fits[[status$best_ind]]
        }
        names(status$selected_vars) <- var_nms[status$selected_vars + 1]
    } else {
        if (!all(order %in% var_nms))
            stop("unknown variable name in 'order'; ",
                 "allowed values: '", paste(var_nms[-1], collapse = "', '"), "'.")
        if (any(order == var_nms[1]))
            stop("response variable '", var_nms[1],
                 "' must not appear in 'order'.")
        for (var in order)
            current_fit <- xtnd_vine(u[, var], current_fit, family_set, selcrit, ...)
        status$selected_vars <- sapply(order, function(nm) which(var_nms == nm) - 1)
    }

    ## adjust model matrix and names
    reorder <- status$selected_vars
    reorder[order(reorder)] <- seq_along(status$selected_vars)
    current_fit$vine$matrix <- gen_dvine_mat(elements = c(1, reorder + 1))

    ## return results
    out <- list(margins = margin_models,
                vine = current_fit$vine,
                order = var_nms[status$selected_vars + 1],
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
        margs <- lapply(seq_len(d), function(i) NULL)
    } else {
        fit_margin <- function(k) {
            arg_lst <- list(
                x = x[, k],
                xmin = par_1d$xmin[k],
                xmax = par_1d$xmax[k],
                bw   = par_1d$bw[k],
                mult = par_1d$mult[k]
            )
            arg_lst[sapply(arg_lst, is.null)] <- NULL
            m <- do.call(kde1d, arg_lst)
            m$x_cc <- x[, k]
            m
        }
        if (cores > 1) {
            k <- NULL  #  for CRAN checks
            margs <- foreach::foreach(k = seq_len(d)) %dopar% fit_margin(k)
        } else {
            margs <- lapply(seq_len(d), fit_margin)
        }
    }
}

get_pits <- function(x, margin_models, cores) {
    if (is.null(margin_models[[1]])) {
        # data are uniform, no need for PIT
        u <- x
    } else {
        get_pit <- function(m) pkde1d(m$x_cc, m)
        m <- NULL  # for cran checks
        if (cores > 1) {
            u <- foreach::foreach(m = margin_models) %dopar% get_pit(m)
        } else {
            u <- lapply(margin_models, get_pit)
        }
        u <- do.call(cbind, u)
        colnames(u) <- colnames(x)
    }

    u
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

initialize_status <- function(d, selcrit) {
    list(
        # remaining variable indices to select from
        remaining_vars = seq_len(d - 1),
        # variables indices included in the model
        selected_vars = NULL,
        # selection criterion
        selcrit = selcrit,
        # current fit criterion (cll + correction)
        current_crit = -Inf,
        # TRUE when no improvement is possible
        optimum_found = FALSE
    )
}

update_status <- function(status, new_vines) {
    crits <- sapply(new_vines, calculate_crit, selcrit = status$selcrit)

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

calculate_crit <- function(fit, selcrit) {
    crit <- fit$cll
    crit <- crit - switch(
        selcrit,
        "loglik" = 0,
        "aic" = fit$vine$npars,
        "bic" = fit$vine$npars * log(dim(fit$psobs$direct)[3]) / 2
    )
    crit
}

#' @importFrom utils modifyList
#' @importFrom stats cor
#' @importFrom rvinecopulib bicop_dist
xtnd_vine <- function(new_var, old_fit, family_set, selcrit, ...) {
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
    for (i in rev(seq_len(d - 1))) {
        # get data for current edge
        zr1 <- psobs$direct[i + 1, i, ]
        zr2 <- if (i == d - 1) {
            psobs$direct[i + 1, i + 1, ]
        } else {
            psobs$indirect[i + 1, i + 1, ]
        }

        # correct bandwidth for regression context
        # (optimal rate is n^(-1/5) instead of n^(-1/6))
        n <- length(zr2)
        dots <- modifyList(list(mult = n^(1/6 - 1/5)), list(...))
        args <- modifyList(
            list(
                data = cbind(zr2, zr1),
                family_set = family_set,
                selcrit = selcrit
            ),
            dots
        )
        args$threshold <- NULL

        # fit
        if (!is.null(dots$threshold)) {
            if (abs(cor(args$data, method = "kendall")[1, 2]) < dots$threshold) {
                pc_fit <- bicop_dist()
                class(pc_fit) <- c("bicop", "bicop_dist")
                pc_fit$data <- args$data
            } else {
                pc_fit <- do.call(bicop, args)
            }
        } else {
            pc_fit <- do.call(bicop, args)
        }
        old_fit$vine$pair_copulas[[d - i]][[i]] <- pc_fit
        npars <- npars + pc_fit$npars

        # pseudo observations for next tree
        psobs$direct[i, i, ] <- hbicop(cbind(zr2, zr1), 1, pc_fit)
        psobs$indirect[i, i, ] <- hbicop(cbind(zr2, zr1), 2, pc_fit)
    }
    vine <- vinecop_dist(old_fit$vine$pair_copulas, gen_dvine_mat(d))
    cll <- old_fit$cll + sum(log(dbicop(cbind(zr2, zr1), pc_fit)))
    list(vine = vine, psobs = psobs, cll = cll)
}

