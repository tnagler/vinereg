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

test_that("model information methods work", {
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), fit$stats$cll)
  expect_equal(attr(ll, "df"), fit$stats$edf)
  expect_equal(attr(ll, "nobs"), fit$stats$nobs)
  expect_equal(nobs(fit), fit$stats$nobs)
  expect_equal(formula(fit), fit$formula)
  expect_equal(model.frame(fit), fit$model_frame)
  expect_equal(AIC(fit), fit$stats$caic)
  expect_equal(BIC(fit), fit$stats$cbic)
})

test_that("plot_effects()", {
  expect_s3_class(plot_effects(fit, NA), "gg")
  expect_error(plot_effects(fit, vars = "asdf"))
})
