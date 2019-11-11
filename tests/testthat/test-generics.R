context("generics")

dat <- data.frame(y = runif(20), x = replicate(2, runif(20)))
fit <- vinereg(y ~ ., dat, order = paste0("x.", 1:2))

test_that("print() works", {
  expect_output(test <- print(fit))
  expect_equal(test, fit)
})

test_that("summary() works", {
  expect_silent(smr <- summary(fit))
  smr_vars <- c("var", "edf", "cll", "caic", "cbic", "p_value")
  expect_equal(colnames(smr), smr_vars)
  expect_equal(nrow(smr), 3)
})

test_that("plot_effects()", {
  expect_s3_class(plot_effects(fit, NA), "gg")
  expect_error(plot_effects(fit, vars = "asdf"))
})
