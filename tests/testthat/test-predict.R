context("predict.vinereg()")

# simulate data
set.seed(3)
x <- matrix(rnorm(30), 10, 3)
y <- x %*% c(1, -1, 2)
dat <- data.frame(y = y, x = x, z = factor(rbinom(20, 3, 0.5), 0:3))
fit <- vinereg(y ~ ., dat, family = "gauss", selcrit = "loglik")

test_that("catches missing variables", {
  fit <- vinereg(y ~ ., dat[1:3])
  expect_error(predict(fit, dat[2]))
})

test_that("catches incorrect type", {
  fit <- vinereg(y ~ ., dat[1:2])
  dat[2] <- ordered(1:10)
  expect_error(predict(fit, dat[2]))
})

test_that("catches incorrect levels", {
  levels(dat[[5]]) <- 1:50
  expect_error(predict(fit, dat))
})

test_that("handles alpha correctly", {
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
  expect_equal(fitted(fit), predict(fit, dat))
  expect_equal(
    cbind(fitted(fit, alpha = 0.2), fitted(fit, alpha = 0.8)),
    fitted(fit, alpha = c(0.2, 0.8))
  )
  expect_equal(cor(dat$y, fitted(fit)[[1]]), 1, tol = 1e-1)
})

test_that("works with discrete response", {
  dat$y <- as.ordered(dat$y)
  fit <- vinereg(y ~ ., dat, fam = "tll")
  expect_equal(fitted(fit), predict(fit, dat))
  expect_error(fitted(fit, alpha = NA))
})


test_that("works on uscale", {
  dat[] <- runif(100)
  expect_silent(fit <- vinereg(y ~ ., dat, uscale = TRUE))
  q <- predict(fit, dat)
  expect_true(all(q >= 0 & q <= 1))
})

