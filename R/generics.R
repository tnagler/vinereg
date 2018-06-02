#' @export
print.vinereg <- function(x, ...) {
    cat("D-vine regression model: ")
    n_predictors <- length(x$order)
    if (n_predictors <= 10) {
        predictors <- paste(x$order, collapse = ", ")
    } else {
        predictors <- paste(x$order[1:10], collapse = ", ")
        predictors <- paste0(predictors, ", ... (", n_predictors - 10, " more)")
    }
    cat(names(x$model_frame)[1], "|", predictors, "\n")
    stats <- unlist(x$stats[1:5])
    stats <- paste(names(stats), round(stats, 2), sep = " = ")
    cat(paste(stats, collapse = ", "))
}

#' @export
summary.vinereg <- function(object, ...) {
    data.frame(
        var = c(names(object$model_frame)[1], object$order),
        edf = object$stats$var_edf,
        cll = object$stats$var_cll,
        caic = object$stats$var_caic,
        cbic = object$stats$var_cbic,
        p_value = object$stats$var_p_value
    )
}
