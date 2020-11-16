#' D-vine regression models
#'
#' Sequential estimation of a regression D-vine for the purpose of quantile
#' prediction as described in Kraus and Czado (2017).
#'
#' If discrete variables are declared as `ordered()` or `factor()`, they are
#' handled as described in Panagiotelis et al. (2012). This is different from
#' previous version where the data was jittered before fitting.
#'
#' @param formula an object of class "formula"; same as [lm()].
#' @param data data frame (or object coercible by [as.data.frame()]) containing
#'   the variables in the model.
#' @param family_set see `family_set` argument of [rvinecopulib::bicop()].
#' @param selcrit selection criterion based on conditional log-likelihood.
#'   \code{"loglik"} (default) imposes no correction; other choices are
#'   \code{"aic"} and \code{"bic"}.
#' @param order the order of covariates in the D-vine, provided as vector of
#'   variable names (after calling
#'   `vinereg:::expand_factors(model.frame(formula, data))`); selected
#'   automatically if `order = NA` (default).
#' @param par_1d list of options passed to [kde1d::kde1d()], must be one value
#'   for each margin, e.g. `list(xmin = c(0, 0, NaN))` if the response and first
#'   covariate have non-negative support.
#' @param weights optional vector of weights for each observation.
#' @param cores integer; the number of cores to use for computations.
#' @param ... further arguments passed to [rvinecopulib::bicop()].
#'
#' @return An object of class vinereg. It is a list containing the elements
#'   \describe{ \item{formula}{the formula used for the fit.}
#'   \item{selcrit}{criterion used for variable selection.}
#'   \item{model_frame}{the data used to fit the regression model.}
#'   \item{margins}{list of marginal models fitted by [kde1d::kde1d()].}
#'   \item{vine}{an [rvinecopulib::vinecop_dist()] object containing the fitted
#'   D-vine.} \item{stats}{fit statistics such as conditional
#'   log-likelihood/AIC/BIC and p-values for each variable's contribution.}
#'   \item{order}{order of the covariates chosen by the variable selection
#'   algorithm.} \item{selected_vars}{indices of selected variables.} } Use
#'   [predict.vinereg()] to predict conditional quantiles. `summary.vinereg()`
#'   shows the contribution of each selected variable with the associated
#'   p-value derived from a likelihood ratio test.
#'
#' @references
#'
#' Kraus and Czado (2017), D-vine copula based quantile regression,
#' Computational Statistics and Data Analysis, 110, 1-18
#'
#' Panagiotelis, A., Czado, C., & Joe, H. (2012). Pair copula constructions for
#' multivariate discrete data. Journal of the American Statistical Association,
#' 107(499), 1063-1072.
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
#' mu_hat <- predict(fit, newdata = dat, alpha = NA) # mean
#' med_hat <- predict(fit, newdata = dat, alpha = 0.5) # median
#'
#' # observed vs predicted
#' plot(cbind(y, mu_hat))
#'
#' ## fixed variable order (no selection)
#' (fit <- vinereg(y ~ ., dat, order = c("x.2", "x.1", "z.1")))
#' @seealso \code{\link{predict.vinereg}}
#'
#' @export
#'
#' @importFrom kde1d kde1d pkde1d
#' @importFrom stats model.frame logLik
#' @importFrom utils modifyList
#' @importFrom rvinecopulib bicop vinecop dvine_structure
#' @importFrom Rcpp sourceCpp
#' @useDynLib vinereg, .registration = TRUE
vinereg <- function(formula, data, family_set = "parametric", selcrit = "aic",
                    order = NA, par_1d = list(), weights = numeric(),
                    cores = 1, ...) {
  # remove unused variables
  if (!missing(data)) {
    mf <- model.frame(formula, data)
  } else {
    mf <- model.frame(formula, parent.frame())
  }
  if (!(is.ordered(mf[[1]]) | is.numeric(mf[[1]]))) {
    stop("response must be numeric or ordered")
  }

  # expand factors and deduce variable types
  mfx <- expand_factors(mf)
  d <- ncol(mfx)
  var_types <- rep("c", d)
  var_types[sapply(mfx, is.ordered)] <- "d"

  ## prepare fit controls (little hack calling bicop() for all checks)
  arg <- list(
    data = t(c(0.5, 0.5)),
    family_set = family_set,
    selcrit = selcrit,
    cores = cores,
    par_method = "mle",
    nonpar_method = "quadratic",
    mult = 1,
    weights = weights,
    psi0 = 0.9,
    presel = TRUE,
    keep_data = FALSE
  )
  ctrl <- do.call(
    bicop,
    modifyList(arg, list(...))
  )$controls

  if (!all(is.na(order))) {
    check_order(order, names(mfx))

    selected_vars <- which(names(mfx) %in% order)
    var_types <- var_types[c(1, selected_vars)]
    mfx <- mfx[, c(1, selected_vars)]

    par_1d <- process_par_1d(mfx, par_1d)
    margins <- fit_margins_cpp(prep_for_kde1d(mfx),
                               sapply(mfx, nlevels),
                               mult = par_1d$mult,
                               xmin = par_1d$xmin,
                               xmax = par_1d$xmax,
                               bw = par_1d$bw,
                               deg = par_1d$deg,
                               weights = weights,
                               cores)
    margins <- finalize_margins(margins, mfx)
    u <- to_uscale(mfx, margins)

    # now we need the correct ordering in selected_vars
    selected_vars <- sapply(order, function(x) which(x == names(mfx)))
    args <- append(
      ctrl,
      list(
        data = u,
        var_types = var_types,
        cores = cores,
        structure = dvine_structure(rank(c(1, selected_vars)))
      )
    )
    fit <- list(
      vine = do.call(vinecop, args),
      selected_vars = selected_vars
    )
  } else {
    par_1d <- process_par_1d(mfx, par_1d)
    margins <- fit_margins_cpp(prep_for_kde1d(mfx),
                               sapply(mfx, nlevels),
                               mult = par_1d$mult,
                               xmin = par_1d$xmin,
                               xmax = par_1d$xmax,
                               bw = par_1d$bw,
                               deg = par_1d$deg,
                               weights = weights,
                               cores)
    margins <- finalize_margins(margins, mfx)
    u <- to_uscale(mfx, margins)
    args <- append(
      ctrl,
      list(data = u, var_types = var_types, cores = cores)
    )
    fit <- do.call(select_dvine_cpp, args)
    margins <- margins[c(1, sort(fit$selected_vars))] # other margins useless
  }

  finalize_vinereg_object(
    formula = formula,
    selcrit = selcrit,
    model_frame = mf,
    margins = margins,
    vine = fit$vine,
    selected_vars = fit$selected_vars,
    var_nms = colnames(mfx)
  )
}

#' @noRd
#' @importFrom stats pchisq
#' @importFrom rvinecopulib as_rvine_structure
finalize_vinereg_object <- function(formula, selcrit, model_frame, margins, vine,
                                    selected_vars, var_nms) {
  vine$names <- c(var_nms[1], var_nms[sort(selected_vars)])
  nobs <- nrow(model_frame)
  vine$nobs <- nobs
  var_edf <- c(
    margins[[1]]$edf,
    sapply(vine$pair_copulas, function(pcs) pcs[[1]]$npars)
  )
  var_cll <- c(
    margins[[1]]$loglik,
    sapply(vine$pair_copulas, function(pcs) pcs[[1]]$loglik)
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

  out <- list(
    formula = formula,
    selcrit = selcrit,
    model_frame = model_frame,
    margins = margins,
    vine = vine,
    stats = stats,
    order = var_nms[selected_vars],
    selected_vars = selected_vars
  )
  class(out) <- "vinereg"
  out
}

check_order <- function(order, var_nms) {
  stopifnot(length(order) > 0)
  if (!all(order %in% var_nms)) {
    stop(
      "unknown variable name in 'order'; ",
      "allowed values: '", paste(var_nms[-1], collapse = "', '"), "'."
    )
  }
  if (any(order == var_nms[1])) {
    stop(
      "response variable '", var_nms[1],
      "' must not appear in 'order'."
    )
  }
}
