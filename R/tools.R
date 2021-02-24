#' brings newdata in a format appropriate for applying rvinecopulib functions.
#' @noRd
prepare_newdata <- function(newdata, object, use_response = FALSE) {
  newdata <- as.data.frame(newdata)
  if (!use_response) {
    object$model_frame <- object$model_frame[-1]
  }
  check_newdata(newdata, object)

  # factors must be expanded to dummy numeric variables
  newdata <- expand_factors(newdata)
  newdata <- remove_unused(newdata, object, use_response)
}

#' checks if newdata has appropriate columns and sorts according to the order
#' used for fitting.
#' @noRd
check_newdata <- function(newdata, object) {
  check_var_availability(newdata, names(object$model_frame))

  # the check_x functions expect variables in newdata and model_frame in
  # the same order
  newdata <- newdata[names(object$model_frame)]
  check_types(newdata, object$model_frame)
  check_levels(newdata, object$model_frame)
}

#' checks if all *selected* covariates are in newdata.
#' @noRd
check_var_availability <- function(newdata, vars) {
  vars_avail <- match(vars, colnames(newdata))
  if (any(is.na(vars_avail))) {
    vars_missing <- paste(vars[is.na(vars_avail)], collapse = ", ")
    stop("'newdata' is missing variables ", vars_missing)
  }
}

#' checks if variable types are equal in original data and new data.
#' @importFrom utils capture.output
#' @noRd
check_types <- function(actual, expected) {
  different_type <- sapply(
    seq_along(actual),
    function(i) !identical(class(actual[[i]])[1], class(expected[[i]])[1])
  )
  if (any(different_type)) {
    errors <- data.frame(
      expected = sapply(actual[different_type], function(x) class(x)[1]),
      actual = sapply(expected[different_type], function(x) class(x)[1])
    )
    errors <- paste(capture.output(print(errors)), collapse = "\n")
    stop("some columns have incorrect type:\n", errors, call. = FALSE)
  }
}

#' checks if factor levels are equal in original data and new data.
#' @noRd
check_levels <- function(actual, expected) {
  # only check factors
  actual <- actual[sapply(actual, is.factor)]
  expected <- expected[sapply(expected, is.factor)]
  if (length(expected) == 0) {
    return(TRUE)
  }

  different_levels <- sapply(
    seq_along(actual),
    function(i) !identical(levels(actual[[i]]), levels(expected[[i]]))
  )
  if (any(different_levels)) {
    errors <- data.frame(
      expected = sapply(
        actual[different_levels],
        function(x) paste(levels(x), collapse = ",")
      ),
      actual = sapply(
        expected[different_levels],
        function(x) paste(levels(x), collapse = ",")
      )
    )
    errors <- paste(capture.output(print(errors)), collapse = "\n")
    stop("some factors have incorrect levels\n", errors, call. = FALSE)
  }
}

#' removes unused variables and returns newdata in the order used for fitting.
#' @noRd
remove_unused <- function(newdata, object, use_response = FALSE) {
  # x must be sorted in the order of the data used for fitting
  fit_order <- object$order[order(object$selected_vars)]
  if (use_response) {
    fit_order <- c(names(object$model_frame)[1], fit_order)
  }
  newdata[, fit_order, drop = FALSE]
}

#' transforms data to uniform scale with probability integral transform.
#' For discrete variables, the output has dimension 2*d
#' @noRd
to_uscale <- function(data, margins, add_response = FALSE) {
  u_sub <- list()
  u <- lapply(seq_along(margins), function(k) pkde1d(data[[k]], margins[[k]]))

  if (any(sapply(margins, function(m) nlevels(m$x) > 0))) {
    compute_u_sub <- function(k) {
      if (nlevels(margins[[k]]$x) > 0) {
        data[, k] <- ordered(data[, k], levels = levels(margins[[k]]$x))
        lv <- as.numeric(data[, k]) - 1
        lv0 <- which(lv == 0)
        lv[lv0] <- 1
        xlv <- ordered(levels(margins[[k]]$x)[lv],
                       levels = levels(margins[[k]]$x))
        u_sub <- pkde1d(xlv, margins[[k]])
        u_sub[lv0] <- 0
        return(u_sub)
      } else {
        return(u[[k]])
      }
    }
    u_sub <- lapply(seq_along(margins), compute_u_sub)
  }

  if (add_response) {
    u <- c(list(0.5), u)
    if (length(u_sub) > 0)
      u_sub <- c(list(0.5), u_sub)
  }
  u <- truncate_u(cbind(do.call(cbind, u), do.call(cbind, u_sub)))
  if ((length(u) == 1) & (NROW(data) > 1))
    u <- matrix(u, NROW(data))
  u
}

#' ensures that u-scale data does not contain zeros or ones.
#' @noRd
truncate_u <- function(u) {
  pmin(pmax(u, 1e-10), 1 - 1e-10)
}

#' transforms predicted response back to orginal variable scale.
#' @noRd
to_yscale <- function(u, object) {
  nms <- colnames(u)
  u <- lapply(u, qkde1d, obj = object$margins[[1]])
  u <- as.data.frame(u)
  names(u) <- nms
  u
}

#' @importFrom stats model.matrix
#' @noRd
expand_factors <- function(data) {
  if (is.data.frame(data)) {
    data <- lapply(data, function(x) {
      if (is.numeric(x) | is.ordered(x)) {
        return(x)
      }
      lvs <- levels(x)
      x <- model.matrix(~x)[, -1, drop = FALSE]
      x <- as.data.frame(x)
      x <- lapply(x, function(y) ordered(y, levels = 0:1))
      names(x) <- lvs[-1]
      x
    })
  }
  as.data.frame(data)
}


process_par_1d <- function(data, pars) {
  d <- ncol(data)
  if (!is.null(pars$xmin)) {
    if (length(pars$xmin) != d)
      stop("'xmin'  must be a vector with one value for each variable")
  } else {
    pars$xmin = rep(NaN, d)
  }
  if (!is.null(pars$xmax)) {
    if (length(pars$xmax) != d)
      stop("'xmax'  must be a vector with one value for each variable")
  } else {
    pars$xmax = rep(NaN, d)
  }


  if (is.null(pars$bw))
    pars$bw <- NA
  if (length(pars$bw) == 1)
    pars$bw <- rep(pars$bw, d)
  if (is.null(pars$mult))
    pars$mult <- 1
  if (length(pars$mult) == 1)
    pars$mult <- rep(pars$mult, d)

  if (is.null(pars$deg))
    pars$deg <- 2
  if (length(pars$deg) == 1)
    pars$deg <- rep(pars$deg, d)

  check_par_1d(data, pars)
  pars
}


#' @importFrom assertthat assert_that
check_par_1d <- function(data, ctrl) {
  nms <- colnames(data)
  if (is.null(nms)) {
    nms <- as.character(seq_len(ncol(data)))
  }
  lapply(seq_len(NCOL(data)), function(k) {
    msg_var <- paste0("Problem with par_1d for variable ", nms[k], ": ")
    tryCatch(
      assert_that(
        is.numeric(ctrl$mult[k]), ctrl$mult[k] > 0,
        is.numeric(ctrl$xmin[k]), is.numeric(ctrl$xmax[k]),
        is.na(ctrl$bw[k]) | (is.numeric(ctrl$bw[k]) & (ctrl$bw[k] > 0)),
        is.numeric(ctrl$deg[k])
      ),
      error = function(e) stop(msg_var, e$message)
    )

    if (is.ordered(data[, k]) & (!is.nan(ctrl$xmin[k]) | !is.nan(ctrl$xmax[k]))) {
      stop(msg_var, "xmin and xmax are not meaningful for x of type ordered.")
    }

    if (!is.nan(ctrl$xmax[k]) & !is.nan(ctrl$xmin[k])) {
      if (ctrl$xmin[k] > ctrl$xmax[k]) {
        stop(msg_var, "xmin is larger than xmax.")
      }
    }
    if (!is.nan(ctrl$xmin[k])) {
      if (any(data[, k] < ctrl$xmin[k])) {
        stop(msg_var, "not all data are larger than xmin.")
      }
    }
    if (!is.nan(ctrl$xmax[k])) {
      if (any(data[, k] > ctrl$xmax[k])) {
        stop(msg_var, "not all data are samller than xmax.")
      }
    }
    if (!(ctrl$deg[k] %in% 0:2)) {
      stop(msg_var, "deg must be either 0, 1, or 2.")
    }
  })
}

prep_for_kde1d <- function(data) {
  data <- lapply(data, function(x) if (is.ordered(x)) as.numeric(x) - 1 else x)
  do.call(cbind, data)
}

finalize_margins <- function(margins, data) {
  for (k in seq_along(margins)) {
    margins[[k]]$x <- data[[k]]
    margins[[k]]$nobs <- nrow(data)
    margins[[k]]$var_name <- names(margins)[k] <- colnames(data)[k]
  }
  margins
}

