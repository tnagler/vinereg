#' brings newdata in a format appropriate for applying rvinecopulib functions.
#' @noRd
prepare_newdata <- function(newdata, object, use_response = FALSE) {
    newdata <- as.data.frame(newdata)
    if (!use_response)
        object$model_frame <- object$model_frame[-1]
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
    fit_order <- object$order[order(object$selected_vars)]
    if (use_response)
        fit_order <- c(names(object$model_frame)[1], fit_order)
    newdata[, fit_order, drop = FALSE]
}

#' transforms data to uniform scale with probability integral transform.
#' For discrete variables, the output has dimension 2*d
#' @noRd
to_uscale <- function(data, margins, add_response = FALSE) {
    compute_u <- function(k) {
        pkde1d(data[[k]], margins[[k]])
    }
    compute_u_sub <- function(k) {
        if (is.factor(data[, k])) {
            lv <- as.numeric(data[[k]]) - 1
            lv[lv == 0] <- NA
            us <- pkde1d(levels(data[[k]])[lv], margins[[k]])
            us[is.na(us) & !is.na(data[[k]])] <- 0
            return(us)
        } else {
            return(u[[k]])
        }
    }
    u_sub <- list()
    if (!is.null(margins[[1]]$u)) {
        # data are uniform, no need for PIT
        u <- lapply(margins, function(m) m$u)
    } else {
        u <- furrr::future_map(seq_along(margins), compute_u)
        if (any(sapply(margins, function(x) length(x$jitter_info$i_disc))))
            u_sub <- furrr::future_map(seq_along(margins), compute_u_sub)
    }
    if (add_response) {
        u <- c(list(0.5), u)
        if (length(u_sub) > 0)
            u_sub <- c(list(0.5), u_sub)
    }
    truncate_u(cbind(do.call(cbind, u), do.call(cbind, u_sub)))
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
            if (is.numeric(x) | is.ordered(x))
                return(x)
            x <- model.matrix(~ x)[, -1, drop = FALSE]
            x <- as.data.frame(x)
            x <- lapply(x, function(y) ordered(y, levels = 0:1))
            names(x) <- seq_along(x)
            x
        })
    }
    as.data.frame(data)
}
