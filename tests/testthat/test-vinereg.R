context("vinereg()")

# simulate data
set.seed(5)
x <- matrix(rnorm(120), 40, 3)
y <- x %*% c(1, -1, 2)
dat <- data.frame(y = y, x = x, z = as.factor(rbinom(40, 3, 0.5)))


test_that("catches wrong arguments", {
  expect_error(vinereg(y ~ ., as.matrix(dat)))
  expect_error(vinereg(dat))
  expect_error(vinereg(asdf ~ ., dat))
  expect_error(vinereg(y ~ ., dat, selcrit = "asdf"))
  expect_error(vinereg(y ~ ., dat, order = c("a", "s", "d", "f")))
  expect_error(vinereg(y ~ ., dat, order = c("y", "x.1")))
  expect_error(vinereg(y ~ ., dat, par_1d = list(deg = 10)))
  expect_error(vinereg(y ~ ., dat, par_1d = list(xmin = 1)))
  expect_error(vinereg(y ~ ., dat, par_1d = list(xmax = 1)))
  expect_error(vinereg(y ~ ., dat, par_1d = list(xmin = rep(100, 7))))
  expect_error(vinereg(y ~ ., dat, par_1d = list(xmax = rep(-100, 7))))
  expect_error(vinereg(y ~ ., dat, par_1d = list(xmin = rep(100, 7),
                                                 xmax = rep(100, 7))))
  expect_error(vinereg(y ~ ., dat, par_1d = list(xmin = rep(100, 5))))
  expect_error(vinereg(y ~ ., dat, par_1d = list(xmax = rep(100, 5))))
  expect_error(vinereg(z ~ .))
})

test_that("all selcrits work", {
  order <- c("x.3", "x.1")
  for (selcrit in c("loglik", "aic", "bic")) {
    expect_equal(
      vinereg(y ~ x.1 + x.3, dat, fam = "gauss", selcrit = selcrit)$order,
      order
    )
  }
})

test_that("works with discrete variables", {
  order <- c("x.1", "x.2", "x.3")
  expect_setequal(
    vinereg(y ~ ., dat, fam = "gauss", selcrit = "bic")$order[1:3],
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
  expect_equal(summary(fit_auto$vine), summary(fit_ord$vine))
})

test_that("works in parallel", {
  fit <- vinereg(y ~ ., dat[-5])
  fit_par <- vinereg(y ~ ., dat[-5], family = "par", cores = 2)
  expect_equal(fit$vine, fit_par$vine)
})

test_that("works on uscale", {
  dat[] <- runif(100)
  expect_silent(fit <- vinereg(y ~ ., dat, uscale = TRUE))
})
