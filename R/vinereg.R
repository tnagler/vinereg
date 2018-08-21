#' D-vine regression models
#'
#' Sequential estimation of a regression D-vine for the purpose of quantile
#' prediction as described in Kraus and Czado (2017). If discrete variables
#' are declared as `ordered()` or `factor()`, jittering is used
#' to make them continuous (see `cctools::cont_conv()`). Although this may
#' make the model estimate inconsistent, predictions are usually still reasonable.
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
#'   for each margin, e.g. `list(xmin = c(0, 0, NaN))` if the response and
#'   first covariate have non-negative support.
#' @param cores integer; the number of cores to use for computations.
#' @param uscale logical indicating whether the data are already on copula scale
#'   (no margins have to be fitted).
#' @param ... further arguments passed to [rvinecopulib::bicop()].
#'
#' @return
#' An object of class vinereg. It is a list containing the elements
#' \describe{
#' \item{formula}{the formula used for the fit.}
#' \item{selcrit}{criterion used for variable selection.}
#' \item{model_frame}{the data used to fit the regression model.}
#' \item{margins}{list of marginal models fitted by [kde1d::kde1d()].}
#' \item{vine}{an [rvinecopulib::vinecop_dist()] object containing the fitted
#'   D-vine.}
#' \item{stats}{fit statistics such as conditional log-likelihood/AIC/BIC and
#' p-values for each variable's contribution.}
#' \item{order}{order of the covariates chosen by the variable selection algorithm.}
#' \item{selected_vars}{indices of selected variables.}
#' }
#' Use [predict.vinereg()] to predict conditional
#' quantiles. `summary.vinereg()` shows the contribution of each selected
#' variable with the associated p-value derived from a likelihood ratio test.
#'
#' @references
#'
#' Kraus and Czado (2017), D-vine copula based quantile regression,
#' Computational Statistics and Data Analysis, 110, 1-18
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(300), 100, 2)
#' y <- x %*% c(1, -2)
#' dat <- data.frame(y = y, x = x, z = as.factor(rbinom(100, 2, 0.5)))
#'
#' # fit vine regression model
#' (fit <- vinereg(y ~ ., dat))
#'
#' # inspect model
#' summary(fit)
#' plot_effects(fit)
#'
#' # model predictions
#' mu_hat  <- predict(fit, newdata = dat, alpha = NA)          # mean
#' med_hat <- predict(fit, newdata = dat, alpha = 0.5)         # median
#'
#' # observed vs predicted
#' plot(cbind(y, mu_hat))
#'
#' ## fixed variable order (no selection)
#' (fit <- vinereg(y ~ ., dat, order = c("x.2", "x.1", "z.1")))
#'
#' @seealso \code{\link{predict.vinereg}}
#'
#' @export
#'
#' @importFrom future plan multiprocess availableCores
#' @importFrom furrr future_map
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
    if (!(is.ordered(mf[[1]]) | is.numeric(mf[[1]])))
        stop("response must be numeric or ordered")

    # expand factors and add noise to discrete variable
    x <- cctools::cont_conv(mf)
    d <- ncol(x)

    ## register parallel backend
    cores <- min(cores, future::availableCores())
    if (cores > 1) {
        future::plan(future::multiprocess, workers = cores)
        on.exit(future::plan(), add = TRUE)
    } else {
        future::plan(future::sequential)
    }

    ## estimation of the marginals and transformation to copula data
    margin_models <- fit_margins(x, par_1d, cores, uscale)
    u <- get_pits(margin_models, cores)
    var_nms <- colnames(u)

    ## initialization
    current_fit <- initialize_fit(u)
    status <- initialize_status(margin_models, selcrit)

    ## estimation --------
    if (any(is.na(order))) {
        ## automatic variable selection
        for (i in seq_len(d - 1)) {
            new_fits <- furrr::future_map(
                status$remaining_vars + 1,
                ~ xtnd_vine(u[, .], current_fit, family_set, selcrit, ...)
            )
            status <- update_status(status, new_fits)
            if (status$optimum_found)
                break
            current_fit <- new_fits[[status$best_ind]]
        }
        if (length(status$selected_vars) > 0)
            names(status$selected_vars) <- var_nms[status$selected_vars + 1]
    } else {
        ## fixed variable order
        check_order(order, var_nms)
        status$selcrit <- "keep_all"
        for (var in order) {
            current_fit <- xtnd_vine(u[, var], current_fit, family_set, selcrit, ...)
            status <- update_status(status, list(current_fit))
        }
        status$selected_vars <- sapply(order, function(nm) which(var_nms == nm) - 1)
        status$remaining_vars <- numeric(0)
    }

    finalize_vinereg_object(
        formula = formula,
        model_frame = mf,
        margins = margin_models,
        vine = current_fit$vine,
        status = status,
        var_nms = var_nms
    )
}

#' @noRd
#' @importFrom stats pchisq
#' @importFrom rvinecopulib as_rvine_structure
finalize_vinereg_object <- function(formula, model_frame, margins, vine,
                                    status, var_nms) {
    ## adjust model matrix and names
    reorder <- status$selected_vars
    reorder[order(reorder)] <- seq_along(status$selected_vars)
    vine$structure <- as_rvine_structure(
        gen_dvine_mat(elements = c(1, reorder + 1))
    )
    vine$names <- c(var_nms[1], names(reorder)[reorder])

    ## compute fit statistics
    nobs <- nrow(model_frame)
    var_edf <- status$edf
    var_cll <- status$cll
    var_caic <- -2 * var_cll + 2 * var_edf
    var_cbic <- -2 * var_cll + log(nobs) * var_edf
    var_p_value <- pchisq(2 * var_cll, var_edf, lower.tail = FALSE)
    var_p_value[1] <- NA
    cll <- sum(var_cll)
    edf <- sum(var_edf)
    caic <- sum(var_caic)
    cbic <- sum(var_cbic)

    stats <- list(
        nobs = nobs,
        edf = edf,
        cll = cll,
        caic = caic,
        cbic = cbic,
        var_edf = var_edf,
        var_cll = var_cll,
        var_caic = var_caic,
        var_cbic = var_cbic,
        var_p_value = var_p_value
    )

    ## return results as S3 object
    out <- list(
        formula = formula,
        selcrit = status$selcrit,
        model_frame = model_frame,
        margins = margins,
        vine = vine,
        stats = stats,
        order = var_nms[status$selected_vars + 1],
        selected_vars = status$selected_vars
    )
    class(out) <- "vinereg"
    out
}


check_order <- function(order, var_nms) {
    if (!all(order %in% var_nms))
        stop("unknown variable name in 'order'; ",
             "allowed values: '", paste(var_nms[-1], collapse = "', '"), "'.")
    if (any(order == var_nms[1]))
        stop("response variable '", var_nms[1],
             "' must not appear in 'order'.")
}


fit_margins <- function(x, par_1d, cores, uscale) {
    d <- ncol(x)
    par_1d <- process_par_1d(par_1d, d)
    if (uscale) {
        # data are uniform, no need to estimate margins
        margs <- lapply(
            seq_len(d),
            function(i) list(u = x[, i], loglik = 0, edf = 0)
        )
    } else {
        fit_margin <- function(k) {
            arg_lst <- list(
                x = x[, k],
                xmin = par_1d$xmin[k],
                xmax = par_1d$xmax[k],
                bw   = par_1d$bw[k],
                mult = par_1d$mult[k],
                deg  = par_1d$deg[k]
            )
            arg_lst[sapply(arg_lst, is.null)] <- NULL
            m <- do.call(kde1d, arg_lst)
            m$x_cc <- x[, k]
            m
        }
        margs <- furrr::future_map(seq_len(d), fit_margin)
    }

    names(margs) <- colnames(x)
    margs
}

get_pits <- function(margin_models, cores) {
    if (!is.null(margin_models[[1]]$u)) {
        # data are uniform, no need for PIT
        u <- sapply(margin_models, function(m) m$u)
    } else {
        get_pit <- function(m) pkde1d(m$x_cc, m)
        u <- furrr::future_map(margin_models, ~ pkde1d(.$x_cc, .))
        u <- do.call(cbind, u)
    }

    colnames(u) <- names(margin_models)
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

    if (is.null(pars$mult))
        pars$mult <- 1
    if (length(pars$mult) == 1)
        pars$mult <- rep(pars$mult, d)
    if (length(pars$mult) != d)
        stop("mult has to be of length 1 or the number of variables")

    if (is.null(pars$deg))
        pars$deg <- 2
    if (length(pars$deg) == 1)
        pars$deg <- rep(pars$deg, d)
    if (length(pars$deg) != d)
        stop("deg has to be of length 1 or the number of variables")

    pars
}

# u_k <- pkde1d(x[, k], fit)
# list(fit = fit, u = pkde1d(x[, k], fit))
#
initialize_fit <- function(u) {
    list(
        # 1-dimensional (= empty) vine
        vine = list(pair_copulas = list(list()), structure = as.matrix(1)),
        # array for storing pseudo-observations
        psobs = list(
            direct = array(u[, 1], dim = c(1, 1, nrow(u))),
            indirect = array(NA, dim = c(1, 1, nrow(u)))
        ),
        # conditional log-likelihood of the model
        cll = 0
    )
}

initialize_status <- function(margin_fits, selcrit) {
    list(
        # remaining variable indices to select from
        remaining_vars = seq_len(length(margin_fits) - 1),
        # variables indices included in the model
        selected_vars = NULL,
        # selection criterion
        selcrit = selcrit,
        # conditional logliklihood (only unconditional margin so for)
        clls = margin_fits[[1]]$loglik,
        # number of parameters in current model
        edf = margin_fits[[1]]$edf,
        # TRUE when no improvement is possible
        optimum_found = FALSE
    )
}

update_status <- function(status, new_vines) {
    clls <- sapply(new_vines, function(vine) vine$cll)
    edf <- sapply(new_vines, function(vine) vine$edf)
    n <- nrow(new_vines[[1]]$psobs[[1]])
    crits <- calculate_crits(clls, edf, n, status$selcrit)

    if (max(crits) < 0) {
        # optimum found, keep old fit
        status$optimum_found <- TRUE
    } else {
        status$best_ind <- which.max(crits)
        status$selected_vars <- c(
            status$selected_vars,
            status$remaining_vars[status$best_ind]
        )
        status$remaining_vars <- setdiff(
            status$remaining_vars,
            status$selected_vars
        )
        status$clls = c(status$clls, clls[status$best_ind])
        status$edf = c(status$edf, edf[status$best_ind])
    }

    status
}

calculate_crits <- function(clls, edf, n, selcrit) {
    clls - switch(
        selcrit,
        "loglik" = 0,
        "aic" = edf,
        "bic" = edf * log(n) / 2,
        "keep_all" = -Inf
    )
}

#' @importFrom utils modifyList
#' @importFrom stats cor
#' @importFrom rvinecopulib bicop_dist
xtnd_vine <- function(new_var, old_fit, family_set, selcrit, ...) {
    d <- dim(old_fit$psobs$direct)[1] + 1
    n <- length(new_var)

    psobs <- list(
        direct = array(NA, dim = c(d, d, n)),
        indirect = array(NA, dim = c(d, d, n))
    )
    psobs$direct[-1, -d, ] <- old_fit$psobs$direct
    psobs$indirect[-1, -d, ] <- old_fit$psobs$indirect
    psobs$direct[d, d, ] <- new_var

    old_fit$vine$pair_copulas[[d - 1]] <- list()
    edf <- 0
    for (i in rev(seq_len(d - 1))) {
        # get data for current edge
        u_e <- matrix(NA, n, 2)
        u_e[, 1] <- psobs$direct[i + 1, i, ]
        u_e[, 2] <- if (i == d - 1) {
            psobs$direct[i + 1, i + 1, ]
        } else {
            psobs$indirect[i + 1, i + 1, ]
        }

        # correct bandwidth for regression context
        # (optimal rate is n^(-1/5) instead of n^(-1/6))
        mult <- ifelse(is.null(list(...)$mult), 1, list(...)$mult)
        dots <- modifyList(list(mult = n^(1/6 - 1/5) * mult), list(...))
        args <- modifyList(
            list(
                data = u_e,
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
        edf <- edf + pc_fit$npars

        # pseudo observations for next tree
        psobs$direct[i, i, ] <- hbicop(u_e, 2, pc_fit)
        psobs$indirect[i, i, ] <- hbicop(u_e, 1, pc_fit)
    }

    list(
        vine = vinecop_dist(old_fit$vine$pair_copulas, gen_dvine_mat(d)),
        psobs = psobs,
        cll = logLik(pc_fit),
        edf = edf
    )
}

