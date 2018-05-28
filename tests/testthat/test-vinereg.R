context("vinereg() function")

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
    expect_warning(
        expect_setequal(
            vinereg(y ~ ., dat, fam = "gauss", selcrit = "bic")$order,
            order
        )
    )
})

test_that("works with fixed order", {
    order <- c("x.3", "x.1", "x.2")
    expect_equal(
        vinereg(y ~ ., dat[-5], fam = "gauss", order = order)$order,
        order
    )

    fit_auto <- vinereg(y ~ ., dat[-5], selcrit = "bic")
    fit_ord <- vinereg(y ~ ., dat[-5], selcrit = "bic", order = fit_auto$order)
    expect_equal(fit_auto, fit_ord)
})

context("methods")

test_that("predict method works", {
    expect_warning(fit <- vinereg(y ~ ., dat, selcrit = "bic"))
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

