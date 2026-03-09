test_that("proper infinite-bound route remains finite", {
  set.seed(42)
  n <- 100
  x <- runif(n)
  y <- x + rnorm(n)

  fit <- bkcde(
    x = x,
    y = y,
    proper = TRUE,
    y.ub = Inf,
    y.lb = -Inf,
    x.ub = Inf,
    x.lb = -Inf,
    degree.min = 0,
    degree.max = 0,
    nmulti = 1,
    optim.cores = "manual",
    optim.degree.cores = 1,
    optim.nmulti.cores = 1,
    optim.ksum.cores = 1,
    fitted.cores = 1,
    proper.cores = 1,
    display.warnings = FALSE
  )

  expect_identical(fit$degree, 0L)
  expect_equal(length(fit$h), 2L)
  expect_true(all(is.finite(fit$h)))
  expect_true(all(is.finite(fit$f)))
})

test_that("predict.bkcde accepts newdata data frames", {
  set.seed(7)
  n <- 80
  x <- runif(n, -0.25, 0.25)
  y <- rbeta(n, 1 + x, 1.5 + x)

  fit <- bkcde(
    h = c(0.15, 0.10),
    degree = 0,
    x = x,
    y = y,
    proper = TRUE,
    fitted.cores = 1,
    proper.cores = 1,
    display.warnings = FALSE
  )

  pred <- predict(
    fit,
    data.frame(x = c(0, 0.1), y = c(0.25, 0.5)),
    fitted.cores = 1,
    proper.cores = 1
  )

  expect_true(is.numeric(pred))
  expect_length(pred, 2L)
  expect_true(all(is.finite(pred)))
})
