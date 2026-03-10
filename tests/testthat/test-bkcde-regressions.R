make_objective_fixture <- function(seed = 101, n = 36) {
  set.seed(seed)
  x <- runif(n, min = 0, max = 1)
  mean_y <- plogis(1.5 * sin(2 * pi * x) + 0.5 * x)
  phi <- 18
  shape1 <- pmax(1.2, mean_y * phi)
  shape2 <- pmax(1.2, (1 - mean_y) * phi)
  y <- rbeta(n, shape1 = shape1, shape2 = shape2)
  list(
    x = x,
    y = y,
    x.lb = 0,
    x.ub = 1,
    y.lb = 0,
    y.ub = 1
  )
}

build_objective_args <- function(data,
                                 bwmethod,
                                 proper.cv,
                                 degree,
                                 h = NULL,
                                 n.integrate = 31,
                                 poly.raw = FALSE) {
  ns <- asNamespace("bkcde")
  essdee <- get("EssDee", envir = ns, inherits = FALSE)
  n <- length(data$y)
  if (is.null(h)) {
    h <- essdee(cbind(data$y, data$x)) * n^(-1 / 6)
  }
  y.seq <- if (isTRUE(proper.cv) || identical(bwmethod, "cv.ls")) {
    seq(data$y.lb, data$y.ub, length.out = n.integrate)
  } else {
    NULL
  }
  X <- if (degree > 0L) {
    cbind(1, poly(data$x, raw = poly.raw, degree = degree))
  } else {
    matrix(1, nrow = n, ncol = 1)
  }
  list(
    h = h,
    x = data$x,
    y = data$y,
    y.lb = data$y.lb,
    y.ub = data$y.ub,
    x.lb = data$x.lb,
    x.ub = data$x.ub,
    poly.raw = poly.raw,
    degree = degree,
    n.integrate = n.integrate,
    optim.ksum.cores = 1L,
    cv.penalty.method = "smooth",
    cv.penalty.cutoff = .Machine$double.xmin,
    verbose = FALSE,
    bwmethod = bwmethod,
    proper.cv = proper.cv,
    X = X,
    y.seq = y.seq,
    cv.binned = FALSE,
    n.binned = 100L,
    x.ub.finite = TRUE,
    x.lb.finite = TRUE,
    y.ub.finite = TRUE,
    y.lb.finite = TRUE
  )
}

evaluate_objective <- function(args) {
  objective <- get("bkcde_optim_fn", envir = asNamespace("bkcde"), inherits = FALSE)
  do.call(objective, args)
}

test_that("proper infinite-bound route remains finite", {
  set.seed(42)
  n <- 100
  x <- runif(n)
  y <- x + rnorm(n)

  fit <- bkcde::bkcde(
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

  fit <- bkcde::bkcde(
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

test_that("bkcde_optim_fn matches cv.only objective at fitted bandwidths", {
  fixture <- make_objective_fixture(seed = 202, n = 40)
  cases <- list(
    list(id = "cv_ml_proper_deg0", bwmethod = "cv.ml", proper.cv = TRUE, degree = 0L),
    list(id = "cv_ml_proper_deg1", bwmethod = "cv.ml", proper.cv = TRUE, degree = 1L),
    list(id = "cv_ls_deg0", bwmethod = "cv.ls", proper.cv = FALSE, degree = 0L)
  )

  for (case in cases) {
    fit <- bkcde::bkcde(
      x = fixture$x,
      y = fixture$y,
      x.lb = fixture$x.lb,
      x.ub = fixture$x.ub,
      y.lb = fixture$y.lb,
      y.ub = fixture$y.ub,
      bwmethod = case$bwmethod,
      proper.cv = case$proper.cv,
      cv.only = TRUE,
      degree.min = case$degree,
      degree.max = case$degree,
      nmulti = 1,
      n.integrate = 31,
      optim.cores = "manual",
      optim.degree.cores = 1,
      optim.nmulti.cores = 1,
      optim.ksum.cores = 1,
      fitted.cores = 1,
      proper.cores = 1,
      display.warnings = FALSE,
      progress = FALSE,
      seed = 77
    )

    objective_value <- evaluate_objective(
      build_objective_args(
        data = fixture,
        bwmethod = case$bwmethod,
        proper.cv = case$proper.cv,
        degree = case$degree,
        h = fit$h,
        n.integrate = 31
      )
    )

    expect_equal(
      unname(objective_value),
      unname(fit$value),
      tolerance = 1e-10,
      info = case$id
    )
  }
})

test_that("bkcde_optim_fn oracle values remain fixed on canonical fixtures", {
  fixture <- make_objective_fixture(seed = 101, n = 36)
  oracle_cases <- list(
    cv_ml_plain_deg0 = list(
      bwmethod = "cv.ml",
      proper.cv = FALSE,
      degree = 0L,
      expected = 20.24223990107549
    ),
    cv_ml_proper_deg0 = list(
      bwmethod = "cv.ml",
      proper.cv = TRUE,
      degree = 0L,
      expected = 19.38437099701448
    ),
    cv_ml_proper_deg1 = list(
      bwmethod = "cv.ml",
      proper.cv = TRUE,
      degree = 1L,
      expected = 21.2981516098956
    ),
    cv_ls_deg1 = list(
      bwmethod = "cv.ls",
      proper.cv = FALSE,
      degree = 1L,
      expected = 2.183656494478436
    )
  )

  for (case_id in names(oracle_cases)) {
    case <- oracle_cases[[case_id]]
    objective_value <- evaluate_objective(
      build_objective_args(
        data = fixture,
        bwmethod = case$bwmethod,
        proper.cv = case$proper.cv,
        degree = case$degree,
        n.integrate = 31
      )
    )

    expect_equal(
      unname(objective_value),
      case$expected,
      tolerance = 1e-10,
      info = case_id
    )
  }
})
