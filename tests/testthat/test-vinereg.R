context("vinereg()")

# simulate data
x <- matrix(rnorm(300), 100, 3)
y <- x %*% c(1, -1, 2)
dat <- data.frame(y = y, x = x, z = as.factor(rbinom(100, 3, 0.5)))


test_that("catches wrong arguments", {
    expect_error(vinereg(y ~ ., as.matrix(dat)))
    expect_error(vinereg(dat))
    expect_error(vinereg(asdf ~ ., dat))
    expect_error(vinereg(y ~ ., selcrit = "asdf"))
    expect_error(vinereg(y ~ ., order = c("a", "s", "d", "f")))
    expect_error(vinereg(y ~ ., order = c("y", "x.1")))
    expect_error(vinereg(z ~ .))
})

test_that("all selcrits work", {
    order <- c("x.3", "x.1")
    for (selcrit in c("loglik", "aic", "bic"))
        expect_equal(
            vinereg(y ~ x.1 + x.3, dat, fam = "gauss", selcrit = selcrit)$order,
            order
        )
})

test_that("works with discrete variables", {
    order <- c("x.1", "x.2", "x.3")
    expect_setequal(
        vinereg(y ~ ., dat, fam = "gauss", selcrit = "bic")$order,
        order
    )
    dat$y <- as.ordered(dat$y)
    expect_silent(vinereg(y ~ ., dat, fam = "tll")$order)
})

test_that("works with fixed order", {
    order <- c("x.3", "x.1", "x.2")
    expect_equal(
        vinereg(y ~ ., dat[-5], fam = "gauss", order = order)$order,
        order
    )

    fit_auto <- vinereg(y ~ ., dat[-5], selcrit = "bic")
    fit_ord <- vinereg(y ~ ., dat[-5], selcrit = "bic", order = fit_auto$order)
    expect_equal(fit_auto$vine, fit_ord$vine)
})

test_that("works on uscale", {
    fit <- vinereg(y ~ ., dat[-5])
    u <- vinereg:::get_pits(fit$margins, 1)
    fit_uscale <- vinereg(y ~ ., as.data.frame(u), uscale = TRUE)

    expect_equal(fit$vine, fit_uscale$vine)
})

test_that("works with threshold", {
    expect_silent(vinereg(y ~ ., dat[-5], threshold = 0.3))
})


test_that("works in parallel", {
    fit <- vinereg(y ~ ., dat[-5])
    fit_par <- vinereg(y ~ ., dat[-5], cores = 2)
    expect_equal(fit$vine, fit_par$vine)
})

## -------------------------------------------------------------
context("predict.vinereg()")

test_that("catches missing variables", {
    fit <- vinereg(y ~ ., dat[1:3])
    expect_error(predict(fit, dat[2]))
})

test_that("handles alpha correctly", {
    fit <- vinereg(y ~ ., dat[1:3])
    expect_equal(colnames(predict(fit, alpha = c(NA, 0.5))), c("mean", "0.5"))
    expect_equal(colnames(predict(fit, alpha = c(0.1, 0.5))), c("0.1", "0.5"))
    expect_error(predict(fit, alpha = NULL))
    expect_error(predict(fit, alpha = "0.2"))
    expect_error(predict(fit, alpha = 1.1))
})

test_that("works in bivariate case", {
    fit <- vinereg(y ~ ., dat[1:2])
    expect_silent(predict(fit, dat[2]))
})

test_that("works with continuous response", {
    fit <- vinereg(y ~ ., dat, selcrit = "loglik")
    expect_equal(fitted(fit), predict(fit, dat))
    expect_equal(
        cbind(fitted(fit, alpha = 0.2), fitted(fit, alpha = 0.8)),
        fitted(fit, alpha = c(0.2, 0.8))
    )
    expect_equal(
        unname(coef(lm(dat$y ~ fitted(fit)[[1]]))),
        c(0, 1),
        tol = 1e-1
    )
})

test_that("works with discrete response", {
    dat$y <- as.ordered(dat$y)
    fit <- vinereg(y ~ ., dat, fam = "tll")
    expect_equal(fitted(fit), predict(fit, dat))
    expect_error(fitted(fit, alpha = NA))
})

fit <- vinereg(y ~ ., dat[-5])
u <- vinereg:::get_pits(fit$margins, 1)
fit_uscale <- vinereg(y ~ ., as.data.frame(u), uscale = TRUE)

test_that("works on uscale", {
    expect_warning(
        expect_equal(predict(fit, u, uscale = TRUE), fitted(fit_uscale))
    )
})


## -------------------------------------------------------------
context("generics")

test_that("print() works", {
    expect_output(test <- print(fit))
    expect_equal(test, fit)
    expect_output(print(fit_uscale))
})

test_that("summary() works", {
    expect_silent(smr <- summary(fit))
    expect_silent(smr_uscale <- summary(fit_uscale))

    smr_vars <- c("var", "edf", "cll", "caic", "cbic", "p_value")
    expect_equal(colnames(smr), smr_vars)
    expect_equal(colnames(smr_uscale), smr_vars)
    expect_equal(nrow(smr), 4)
    expect_equal(nrow(smr_uscale),  4)
    expect_equal(unname(unlist(smr_uscale[1, -1])), c(rep(0, 4), NA))
})


test_that("plot_effects()", {
    expect_s3_class(plot_effects(fit, NA), "gg")
    expect_warning(plot_effects(fit_uscale))
    expect_error(plot_effects(fit, vars = "asdf"))
})

