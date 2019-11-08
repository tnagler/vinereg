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
#' @importFrom Rcpp sourceCpp
#' @useDynLib vinereg, .registration = TRUE
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
    mfx <- expand_factors(mf)
    d <- ncol(mfx)
    var_types <- rep("c", d)
    var_types[sapply(mfx, is.ordered)] <- "d"

    ## register parallel backend
    cores <- min(cores, future::availableCores())
    if (cores > 1) {
        suppressWarnings(future::plan(future::multiprocess, workers = cores))
        on.exit(future::plan(), add = TRUE)
    } else {
        future::plan(future::sequential)
    }

    ## estimation of the marginals and transformation to copula data
    margins <- fit_margins(mfx, par_1d, cores, uscale)
    fit <- select_dvine_cpp(to_uscale(mfx, margins), var_types)
browser()
    finalize_vinereg_object(
        formula = formula,
        model_frame = mf,
        margins = margins,
        fit = fit,
        var_nms = colnames(mfx)
    )
}

#' @noRd
#' @importFrom stats pchisq
#' @importFrom rvinecopulib as_rvine_structure
finalize_vinereg_object <- function(formula, model_frame, margins, fit, var_nms) {
    ## adjust model matrix and names
    fit$vine$names <- c(var_nms[1], var_nms[sort(fit$selected_vars + 1)])

    ## compute fit statistics
    nobs <- nrow(model_frame)
    fit$vine$nobs <- nobs
    var_edf <- c(
        margins[[1]]$edf,
        sapply(fit$vine$pair_copulas, function(pcs) pcs[[1]]$npars)
    )
    var_cll <- c(
        margins[[1]]$loglik,
        sapply(fit$vine$pair_copulas, function(pcs) pcs[[1]]$loglik)
    )
    var_caic <- -2 * var_cll + 2 * var_edf
    var_cbic <- -2 * var_cll + log(nobs) * var_edf
    var_p_value <- suppressWarnings(
        pchisq(2 * var_cll, var_edf, lower.tail = FALSE)
    )
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
        selcrit = fit$selcrit,
        model_frame = model_frame,
        margins = margins[1 + c(0, sort(fit$selected_vars))],
        vine = fit$vine,
        stats = stats,
        order = var_nms[fit$selected_vars + 1],
        selected_vars = fit$selected_vars + 1
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
