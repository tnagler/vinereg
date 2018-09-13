#' brings newdata in a format appropriate for applying rvinecopulib functions.
#' @noRd
prepare_newdata <- function(newdata, object, use_response = FALSE) {
    newdata <- as.data.frame(newdata)

    if (!use_response)
        object$model_frame <- object$model_frame[-1]
    check_newdata(newdata, object)

    # factors must be expanded to dummy numeric variables
    newdata <- cctools::expand_as_numeric(newdata)
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
    actual   <- actual[sapply(actual, is.factor)]
    expected <- expected[sapply(expected, is.factor)]
    if (length(expected) == 0)
        return(TRUE)

    different_levels <- sapply(
        seq_along(actual),
        function(i) !identical(levels(actual[[i]]), levels(expected[[i]]))
    )
    if (any(different_levels)) {
        errors <- data.frame(
            expected = sapply(actual[different_levels],
                              function(x) paste(levels(x), collapse = ",")),
            actual = sapply(expected[different_levels],
                            function(x) paste(levels(x), collapse = ","))
        )
        errors <- paste(capture.output(print(errors)), collapse = "\n")
        stop("some factors have incorrect levels\n", errors, call. = FALSE)
    }

}

#' removes unused variables and returns newdata in the order used for fitting.
#' @noRd
remove_unused <- function(newdata, object, use_response = FALSE) {
    # x must be sorted in the order of the data used for fitting
    fit_order <- names(sort(object$selected_vars))
    if (use_response)
        fit_order <- c(names(object$model_frame)[1], fit_order)
    newdata[, fit_order, drop = FALSE]
}


#' checks if margins were estimated.
#' @noRd
adjust_uscale <- function(object, uscale) {
    if (!is.null(object$margins[[1]]$u) & (!uscale)) {
        warning("no margins have been estimated, setting uscale = TRUE")
        uscale <- TRUE
    }
    uscale
}

#' transforms data to uniform scale with probability integral transform.
#' @noRd
to_uscale <- function(x, object) {
    for (var in colnames(x))
        x[, var] <- truncate_u(pkde1d(x[, var], object$margins[[var]]))
    x
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
    if (inherits(object$model_frame[[1]], "ordered")) {
        # when response is discrete, we need to adjust the quantiles
        lvls <- levels(object$model_frame[[1]])
        u <- lapply(u, with_levels, lvls = lvls)
    }

    u <- as.data.frame(u)
    names(u) <- nms
    u
}

#' returns the quantile predictions as order variable with appropriate levels.
#' @noRd
with_levels <- function(q, lvls) {
    q <- ceiling(q)
    q <- pmax(q, 1)
    q <- pmin(q, length(lvls))
    ordered(lvls[q], levels = lvls)
}
