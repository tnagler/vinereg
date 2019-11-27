context("cpit()")

# simulate data
set.seed(3)
x <- matrix(rnorm(30), 10, 3)
y <- x %*% c(1, -1, 2)
dat <- data.frame(y = y, x = x, z = as.factor(rbinom(10, 3, 0.5)))
fit <- vinereg(y ~ ., family = "gauss", dat)

test_that("catches missing variables", {
  expect_error(cpit(fit, dat[1:2]))
})

test_that("catches incorrect type", {
  dat[2] <- ordered(1:10)
  expect_error(cpit(fit, dat))
})

test_that("catches incorrect levels", {
  levels(dat[[5]]) <- 1:50
  expect_error(cpit(fit, dat))
})

test_that("works in bivariate case", {
  fit <- vinereg(y ~ ., dat[1:2])
  expect_silent(cpit(fit, dat[1:2]))
})

test_that("works with continuous response", {
  expect_gt(ks.test(cpit(fit, dat), "punif")$p, 0.01)
})

test_that("works with discrete response", {
  dat$y <- as.ordered(dat$y)
  fit <- vinereg(y ~ ., dat, fam = "tll")
  expect_silent(cpit(fit, dat))
})
