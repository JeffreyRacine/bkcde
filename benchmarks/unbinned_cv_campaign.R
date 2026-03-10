`%||%` <- function(x, y) {
  if (is.null(x) || is.na(x) || identical(x, "")) y else x
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- parts[1]
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
    out[[gsub("-", "_", key)]] <- value
  }
  out
}

require_arg <- function(opts, key) {
  value <- opts[[key]]
  if (is.null(value) || identical(value, "")) {
    stop(sprintf("Missing required option --%s", gsub("_", "-", key)), call. = FALSE)
  }
  value
}

as_int <- function(x, default = NULL) {
  if (is.null(x)) {
    return(default)
  }
  as.integer(x)
}

as_num <- function(x, default = NULL) {
  if (is.null(x)) {
    return(default)
  }
  as.numeric(x)
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("true", "1", "yes", "y")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

run_system <- function(cmd, args, stdout = "", stderr = "") {
  status <- system2(cmd, args = args, stdout = stdout, stderr = stderr)
  if (!identical(status, 0L)) {
    stop(sprintf("Command failed (%s %s) with exit status %s",
                 cmd,
                 paste(args, collapse = " "),
                 status),
         call. = FALSE)
  }
  invisible(status)
}

install_repo <- function(repo, lib, log_path) {
  ensure_dir(lib)
  run_system(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--preclean", sprintf("--library=%s", lib), repo),
    stdout = log_path,
    stderr = log_path
  )
}

make_dataset <- function(seed, n) {
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

make_scenarios <- function() {
  list(
    list(
      id = "cv_ml_plain_deg0_c1",
      bwmethod = "cv.ml",
      proper.cv = FALSE,
      degree = 0L,
      optim.ksum.cores = 1L
    ),
    list(
      id = "cv_ml_proper_deg0_c1",
      bwmethod = "cv.ml",
      proper.cv = TRUE,
      degree = 0L,
      optim.ksum.cores = 1L
    ),
    list(
      id = "cv_ls_deg0_c1",
      bwmethod = "cv.ls",
      proper.cv = FALSE,
      degree = 0L,
      optim.ksum.cores = 1L
    ),
    list(
      id = "cv_ml_plain_deg1_c1",
      bwmethod = "cv.ml",
      proper.cv = FALSE,
      degree = 1L,
      optim.ksum.cores = 1L
    ),
    list(
      id = "cv_ml_proper_deg1_c1",
      bwmethod = "cv.ml",
      proper.cv = TRUE,
      degree = 1L,
      optim.ksum.cores = 1L
    ),
    list(
      id = "cv_ls_deg1_c1",
      bwmethod = "cv.ls",
      proper.cv = FALSE,
      degree = 1L,
      optim.ksum.cores = 1L
    ),
    list(
      id = "cv_ml_proper_deg0_c2",
      bwmethod = "cv.ml",
      proper.cv = TRUE,
      degree = 0L,
      optim.ksum.cores = 2L
    )
  )
}

get_scenarios <- function(selected = NULL) {
  scenarios <- make_scenarios()
  if (is.null(selected) || !length(selected)) {
    return(scenarios)
  }
  ids <- vapply(scenarios, `[[`, character(1), "id")
  keep <- ids %in% selected
  if (!all(selected %in% ids)) {
    missing <- setdiff(selected, ids)
    stop(sprintf("Unknown scenario id(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  scenarios[keep]
}

get_ns_fun <- function(name) {
  get(name, envir = asNamespace("bkcde"), inherits = FALSE)
}

build_objective_args <- function(scenario, dataset, n.integrate, h = NULL, poly.raw = FALSE) {
  x <- dataset$x
  y <- dataset$y
  n <- length(y)
  essdee <- get_ns_fun("EssDee")
  if (is.null(h)) {
    h <- essdee(cbind(y, x)) * n^(-1 / 6)
  }
  y.seq <- if (isTRUE(scenario$proper.cv) || identical(scenario$bwmethod, "cv.ls")) {
    seq(dataset$y.lb, dataset$y.ub, length.out = n.integrate)
  } else {
    NULL
  }
  X <- if (scenario$degree > 0L) {
    cbind(1, poly(x, raw = poly.raw, degree = scenario$degree))
  } else {
    matrix(1, nrow = n, ncol = 1)
  }
  list(
    h = h,
    x = x,
    y = y,
    y.lb = dataset$y.lb,
    y.ub = dataset$y.ub,
    x.lb = dataset$x.lb,
    x.ub = dataset$x.ub,
    poly.raw = poly.raw,
    degree = scenario$degree,
    n.integrate = n.integrate,
    optim.ksum.cores = scenario$optim.ksum.cores,
    cv.penalty.method = "smooth",
    cv.penalty.cutoff = .Machine$double.xmin,
    verbose = FALSE,
    bwmethod = scenario$bwmethod,
    proper.cv = scenario$proper.cv,
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

evaluate_objective_scenario <- function(scenario, dataset, n.integrate, n, h = NULL) {
  objective <- get_ns_fun("bkcde_optim_fn")
  args <- build_objective_args(
    scenario = scenario,
    dataset = dataset,
    n.integrate = n.integrate,
    h = h
  )
  elapsed <- system.time({
    value <- do.call(objective, args)
  })[["elapsed"]]
  data.frame(
    scenario = scenario$id,
    bwmethod = scenario$bwmethod,
    proper_cv = scenario$proper.cv,
    degree_target = scenario$degree,
    optim_ksum_cores = scenario$optim.ksum.cores,
    n = n,
    n_integrate = n.integrate,
    nmulti = NA_integer_,
    elapsed_wall = as.numeric(elapsed),
    secs_optim_elapsed = NA_real_,
    value = as.numeric(value),
    degree = scenario$degree,
    h_y = as.numeric(args$h[1]),
    h_x = as.numeric(args$h[2]),
    convergence = NA_integer_,
    stringsAsFactors = FALSE
  )
}

fit_scenario <- function(scenario, dataset, optim_seed, n.integrate, nmulti, n) {
  set.seed(optim_seed)
  elapsed <- system.time({
    fit <- bkcde::bkcde(
      x = dataset$x,
      y = dataset$y,
      x.lb = dataset$x.lb,
      x.ub = dataset$x.ub,
      y.lb = dataset$y.lb,
      y.ub = dataset$y.ub,
      bwmethod = scenario$bwmethod,
      cv = "full",
      cv.only = TRUE,
      cv.binned = FALSE,
      proper.cv = scenario$proper.cv,
      degree.min = scenario$degree,
      degree.max = scenario$degree,
      nmulti = nmulti,
      n.integrate = n.integrate,
      optim.cores = "manual",
      optim.degree.cores = 1L,
      optim.nmulti.cores = 1L,
      optim.ksum.cores = scenario$optim.ksum.cores,
      fitted.cores = 1L,
      proper.cores = 1L,
      display.warnings = FALSE,
      progress = FALSE,
      seed = optim_seed
    )
  })[["elapsed"]]

  data.frame(
    scenario = scenario$id,
    bwmethod = scenario$bwmethod,
    proper_cv = scenario$proper.cv,
    degree_target = scenario$degree,
    optim_ksum_cores = scenario$optim.ksum.cores,
    n = n,
    n_integrate = n.integrate,
    nmulti = nmulti,
    elapsed_wall = as.numeric(elapsed),
    secs_optim_elapsed = as.numeric(fit$secs.optim.elapsed %||% NA_real_),
    value = as.numeric(fit$value %||% NA_real_),
    degree = as.integer(fit$degree %||% NA_integer_),
    h_y = as.numeric(fit$h[1] %||% NA_real_),
    h_x = as.numeric(fit$h[2] %||% NA_real_),
    convergence = as.integer(fit$convergence %||% NA_integer_),
    stringsAsFactors = FALSE
  )
}

write_lines <- function(lines, path) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

run_profile <- function(repo, out_dir, scenario_id, n, n.integrate, nmulti) {
  lib_dir <- file.path(out_dir, "lib")
  install_log <- file.path(out_dir, "install.log")
  install_repo(repo = repo, lib = lib_dir, log_path = install_log)
  .libPaths(c(lib_dir, .libPaths()))
  suppressPackageStartupMessages(library(bkcde))

  scenarios <- get_scenarios(scenario_id)
  scenario <- Filter(function(x) identical(x$id, scenario_id), scenarios)
  if (length(scenario) != 1L) {
    stop(sprintf("Unknown scenario id for profiling: %s", scenario_id), call. = FALSE)
  }
  scenario <- scenario[[1L]]

  dataset <- make_dataset(seed = 4101L, n = n)
  profile_path <- file.path(out_dir, "Rprof.out")
  summary_path <- file.path(out_dir, "profile_summary.txt")
  result_path <- file.path(out_dir, "profile_result.csv")

  Rprof(profile_path, interval = 0.001)
  on.exit(Rprof(NULL), add = TRUE)
  result <- fit_scenario(
    scenario = scenario,
    dataset = dataset,
    optim_seed = 9001L,
    n.integrate = n.integrate,
    nmulti = nmulti,
    n = n
  )
  Rprof(NULL)

  prof <- summaryRprof(profile_path)
  top_self <- capture.output(print(head(prof$by.self, 20)))
  top_total <- capture.output(print(head(prof$by.total, 20)))
  write_lines(
    c(
      sprintf("repo=%s", normalizePath(repo, winslash = "/", mustWork = TRUE)),
      sprintf("scenario=%s", scenario_id),
      sprintf("n=%d", n),
      sprintf("n_integrate=%d", n.integrate),
      sprintf("nmulti=%d", nmulti),
      "",
      "[by.self]",
      top_self,
      "",
      "[by.total]",
      top_total
    ),
    summary_path
  )
  write.csv(result, result_path, row.names = FALSE)
}

run_objective_profile <- function(repo, out_dir, scenario_id, n, n.integrate) {
  lib_dir <- file.path(out_dir, "lib")
  install_log <- file.path(out_dir, "install.log")
  install_repo(repo = repo, lib = lib_dir, log_path = install_log)
  .libPaths(c(lib_dir, .libPaths()))
  suppressPackageStartupMessages(library(bkcde))

  scenarios <- get_scenarios(scenario_id)
  scenario <- Filter(function(x) identical(x$id, scenario_id), scenarios)
  if (length(scenario) != 1L) {
    stop(sprintf("Unknown scenario id for objective profiling: %s", scenario_id), call. = FALSE)
  }
  scenario <- scenario[[1L]]

  dataset <- make_dataset(seed = 4101L, n = n)
  profile_path <- file.path(out_dir, "Rprof.out")
  summary_path <- file.path(out_dir, "profile_summary.txt")
  result_path <- file.path(out_dir, "profile_result.csv")

  Rprof(profile_path, interval = 0.001)
  on.exit(Rprof(NULL), add = TRUE)
  result <- evaluate_objective_scenario(
    scenario = scenario,
    dataset = dataset,
    n.integrate = n.integrate,
    n = n
  )
  Rprof(NULL)

  prof <- summaryRprof(profile_path)
  top_self <- capture.output(print(head(prof$by.self, 20)))
  top_total <- capture.output(print(head(prof$by.total, 20)))
  write_lines(
    c(
      sprintf("repo=%s", normalizePath(repo, winslash = "/", mustWork = TRUE)),
      sprintf("scenario=%s", scenario_id),
      sprintf("n=%d", n),
      sprintf("n_integrate=%d", n.integrate),
      "mode=objective",
      "",
      "[by.self]",
      top_self,
      "",
      "[by.total]",
      top_total
    ),
    summary_path
  )
  write.csv(result, result_path, row.names = FALSE)
}

run_benchmarks <- function(repo, out_dir, times, seed_mode, n, n.integrate, nmulti, selected = NULL) {
  lib_dir <- file.path(out_dir, "lib")
  install_log <- file.path(out_dir, "install.log")
  install_repo(repo = repo, lib = lib_dir, log_path = install_log)
  .libPaths(c(lib_dir, .libPaths()))
  suppressPackageStartupMessages(library(bkcde))

  scenarios <- get_scenarios(selected)
  results <- vector("list", length(scenarios) * times)
  idx <- 1L
  for (scenario in scenarios) {
    for (rep_idx in seq_len(times)) {
      data_seed <- if (identical(seed_mode, "vary")) 7000L + rep_idx else 7000L
      optim_seed <- if (identical(seed_mode, "vary")) 9000L + rep_idx else 9000L
      dataset <- make_dataset(seed = data_seed, n = n)
      row <- fit_scenario(
        scenario = scenario,
        dataset = dataset,
        optim_seed = optim_seed,
        n.integrate = n.integrate,
        nmulti = nmulti,
        n = n
      )
      row$replicate <- rep_idx
      row$seed_mode <- seed_mode
      row$data_seed <- data_seed
      row$optim_seed <- optim_seed
      results[[idx]] <- row
      idx <- idx + 1L
    }
  }

  results <- do.call(rbind, results)
  results <- results[, c(
    "scenario", "replicate", "seed_mode", "data_seed", "optim_seed",
    "bwmethod", "proper_cv", "degree_target", "optim_ksum_cores",
    "n", "n_integrate", "nmulti", "elapsed_wall", "secs_optim_elapsed",
    "value", "degree", "h_y", "h_x", "convergence"
  )]
  write.csv(results, file.path(out_dir, "results.csv"), row.names = FALSE)

  lines <- c(
    sprintf("repo=%s", normalizePath(repo, winslash = "/", mustWork = TRUE)),
    sprintf("times=%d", times),
    sprintf("seed_mode=%s", seed_mode),
    sprintf("n=%d", n),
    sprintf("n_integrate=%d", n.integrate),
    sprintf("nmulti=%d", nmulti),
    sprintf("cases=%d", length(scenarios))
  )
  write_lines(lines, file.path(out_dir, "run_meta.txt"))
}

run_objective_benchmarks <- function(repo, out_dir, times, seed_mode, n, n.integrate, selected = NULL) {
  lib_dir <- file.path(out_dir, "lib")
  install_log <- file.path(out_dir, "install.log")
  install_repo(repo = repo, lib = lib_dir, log_path = install_log)
  .libPaths(c(lib_dir, .libPaths()))
  suppressPackageStartupMessages(library(bkcde))

  scenarios <- get_scenarios(selected)
  results <- vector("list", length(scenarios) * times)
  idx <- 1L
  for (scenario in scenarios) {
    for (rep_idx in seq_len(times)) {
      data_seed <- if (identical(seed_mode, "vary")) 7000L + rep_idx else 7000L
      dataset <- make_dataset(seed = data_seed, n = n)
      row <- evaluate_objective_scenario(
        scenario = scenario,
        dataset = dataset,
        n.integrate = n.integrate,
        n = n
      )
      row$replicate <- rep_idx
      row$seed_mode <- seed_mode
      row$data_seed <- data_seed
      row$optim_seed <- NA_integer_
      results[[idx]] <- row
      idx <- idx + 1L
    }
  }

  results <- do.call(rbind, results)
  results <- results[, c(
    "scenario", "replicate", "seed_mode", "data_seed", "optim_seed",
    "bwmethod", "proper_cv", "degree_target", "optim_ksum_cores",
    "n", "n_integrate", "nmulti", "elapsed_wall", "secs_optim_elapsed",
    "value", "degree", "h_y", "h_x", "convergence"
  )]
  write.csv(results, file.path(out_dir, "results.csv"), row.names = FALSE)

  lines <- c(
    sprintf("repo=%s", normalizePath(repo, winslash = "/", mustWork = TRUE)),
    sprintf("times=%d", times),
    sprintf("seed_mode=%s", seed_mode),
    sprintf("n=%d", n),
    sprintf("n_integrate=%d", n.integrate),
    sprintf("cases=%d", length(scenarios)),
    "mode=objective"
  )
  write_lines(lines, file.path(out_dir, "run_meta.txt"))
}

run_single_fit <- function(lib_dir, out_csv, scenario_id, data_seed, optim_seed, n, n.integrate, nmulti) {
  .libPaths(c(lib_dir, .libPaths()))
  suppressPackageStartupMessages(library(bkcde))
  scenario <- get_scenarios(scenario_id)[[1L]]
  dataset <- make_dataset(seed = data_seed, n = n)
  result <- fit_scenario(
    scenario = scenario,
    dataset = dataset,
    optim_seed = optim_seed,
    n.integrate = n.integrate,
    nmulti = nmulti,
    n = n
  )
  write.csv(result, out_csv, row.names = FALSE)
}

run_single_objective <- function(lib_dir, out_csv, scenario_id, data_seed, n, n.integrate) {
  .libPaths(c(lib_dir, .libPaths()))
  suppressPackageStartupMessages(library(bkcde))
  scenario <- get_scenarios(scenario_id)[[1L]]
  dataset <- make_dataset(seed = data_seed, n = n)
  result <- evaluate_objective_scenario(
    scenario = scenario,
    dataset = dataset,
    n.integrate = n.integrate,
    n = n
  )
  write.csv(result, out_csv, row.names = FALSE)
}

run_system_stdout <- function(cmd, args) {
  status <- system2(cmd, args = args)
  if (!identical(status, 0L)) {
    stop(sprintf("Command failed (%s %s) with exit status %s",
                 cmd,
                 paste(args, collapse = " "),
                 status),
         call. = FALSE)
  }
  invisible(status)
}

run_interleaved_compare <- function(baseline_repo,
                                    candidate_repo,
                                    out_dir,
                                    times,
                                    seed_mode,
                                    n,
                                    n.integrate,
                                    nmulti,
                                    selected = NULL) {
  baseline_lib <- ensure_dir(file.path(out_dir, "lib-baseline"))
  candidate_lib <- ensure_dir(file.path(out_dir, "lib-candidate"))
  install_repo(baseline_repo, baseline_lib, file.path(out_dir, "install-baseline.log"))
  install_repo(candidate_repo, candidate_lib, file.path(out_dir, "install-candidate.log"))

  scenarios <- get_scenarios(selected)
  this_script <- normalizePath(sys.frame(1)$ofile %||% commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))],
                               winslash = "/",
                               mustWork = FALSE)
  if (!length(this_script) || identical(this_script, "")) {
    stop("Unable to resolve script path for interleaved compare mode", call. = FALSE)
  }
  this_script <- sub("^--file=", "", this_script[1L])

  rows <- vector("list", length(scenarios) * times * 2L)
  idx <- 1L
  for (scenario in scenarios) {
    for (rep_idx in seq_len(times)) {
      data_seed <- if (identical(seed_mode, "vary")) 7000L + rep_idx else 7000L
      optim_seed <- if (identical(seed_mode, "vary")) 9000L + rep_idx else 9000L
      order_tag <- if ((rep_idx %% 2L) == 1L) c("baseline", "candidate") else c("candidate", "baseline")
      row_cache <- list()

      for (tag in order_tag) {
        lib_dir <- if (identical(tag, "baseline")) baseline_lib else candidate_lib
        row_path <- file.path(out_dir, sprintf("%s_%s_rep%03d.csv", scenario$id, tag, rep_idx))
        run_system_stdout(
          file.path(R.home("bin"), "Rscript"),
          c(
            "--no-save",
            this_script,
            "--mode=single-fit",
            sprintf("--lib-dir=%s", lib_dir),
            sprintf("--out-csv=%s", row_path),
            sprintf("--scenario=%s", scenario$id),
            sprintf("--data-seed=%d", data_seed),
            sprintf("--optim-seed=%d", optim_seed),
            sprintf("--n=%d", n),
            sprintf("--n-integrate=%d", n.integrate),
            sprintf("--nmulti=%d", nmulti)
          )
        )
        row <- read.csv(row_path, stringsAsFactors = FALSE)
        row$variant <- tag
        row$replicate <- rep_idx
        row$seed_mode <- seed_mode
        row$data_seed <- data_seed
        row$optim_seed <- optim_seed
        row$order_slot <- match(tag, order_tag)
        row_cache[[tag]] <- row
      }

      rows[[idx]] <- row_cache[["baseline"]]
      rows[[idx + 1L]] <- row_cache[["candidate"]]
      idx <- idx + 2L
    }
  }

  rows <- do.call(rbind, rows)
  rows <- rows[, c(
    "variant", "scenario", "replicate", "seed_mode", "data_seed", "optim_seed", "order_slot",
    "bwmethod", "proper_cv", "degree_target", "optim_ksum_cores",
    "n", "n_integrate", "nmulti", "elapsed_wall", "secs_optim_elapsed",
    "value", "degree", "h_y", "h_x", "convergence"
  )]
  write.csv(rows, file.path(out_dir, "interleaved_rows.csv"), row.names = FALSE)

  baseline <- rows[rows$variant == "baseline", , drop = FALSE]
  candidate <- rows[rows$variant == "candidate", , drop = FALSE]
  write.csv(baseline, file.path(out_dir, "baseline_results.csv"), row.names = FALSE)
  write.csv(candidate, file.path(out_dir, "candidate_results.csv"), row.names = FALSE)
}

run_objective_interleaved_compare <- function(baseline_repo,
                                              candidate_repo,
                                              out_dir,
                                              times,
                                              seed_mode,
                                              n,
                                              n.integrate,
                                              selected = NULL) {
  baseline_lib <- ensure_dir(file.path(out_dir, "lib-baseline"))
  candidate_lib <- ensure_dir(file.path(out_dir, "lib-candidate"))
  install_repo(baseline_repo, baseline_lib, file.path(out_dir, "install-baseline.log"))
  install_repo(candidate_repo, candidate_lib, file.path(out_dir, "install-candidate.log"))

  scenarios <- get_scenarios(selected)
  this_script <- normalizePath(sys.frame(1)$ofile %||% commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))],
                               winslash = "/",
                               mustWork = FALSE)
  if (!length(this_script) || identical(this_script, "")) {
    stop("Unable to resolve script path for interleaved compare mode", call. = FALSE)
  }
  this_script <- sub("^--file=", "", this_script[1L])

  rows <- vector("list", length(scenarios) * times * 2L)
  idx <- 1L
  for (scenario in scenarios) {
    for (rep_idx in seq_len(times)) {
      data_seed <- if (identical(seed_mode, "vary")) 7000L + rep_idx else 7000L
      order_tag <- if ((rep_idx %% 2L) == 1L) c("baseline", "candidate") else c("candidate", "baseline")
      row_cache <- list()

      for (tag in order_tag) {
        lib_dir <- if (identical(tag, "baseline")) baseline_lib else candidate_lib
        row_path <- file.path(out_dir, sprintf("%s_%s_rep%03d.csv", scenario$id, tag, rep_idx))
        run_system_stdout(
          file.path(R.home("bin"), "Rscript"),
          c(
            "--no-save",
            this_script,
            "--mode=single-objective",
            sprintf("--lib-dir=%s", lib_dir),
            sprintf("--out-csv=%s", row_path),
            sprintf("--scenario=%s", scenario$id),
            sprintf("--data-seed=%d", data_seed),
            sprintf("--n=%d", n),
            sprintf("--n-integrate=%d", n.integrate)
          )
        )
        row <- read.csv(row_path, stringsAsFactors = FALSE)
        row$variant <- tag
        row$replicate <- rep_idx
        row$seed_mode <- seed_mode
        row$data_seed <- data_seed
        row$optim_seed <- NA_integer_
        row$order_slot <- match(tag, order_tag)
        row_cache[[tag]] <- row
      }

      rows[[idx]] <- row_cache[["baseline"]]
      rows[[idx + 1L]] <- row_cache[["candidate"]]
      idx <- idx + 2L
    }
  }

  rows <- do.call(rbind, rows)
  rows <- rows[, c(
    "variant", "scenario", "replicate", "seed_mode", "data_seed", "optim_seed", "order_slot",
    "bwmethod", "proper_cv", "degree_target", "optim_ksum_cores",
    "n", "n_integrate", "nmulti", "elapsed_wall", "secs_optim_elapsed",
    "value", "degree", "h_y", "h_x", "convergence"
  )]
  write.csv(rows, file.path(out_dir, "interleaved_rows.csv"), row.names = FALSE)

  baseline <- rows[rows$variant == "baseline", , drop = FALSE]
  candidate <- rows[rows$variant == "candidate", , drop = FALSE]
  write.csv(baseline, file.path(out_dir, "baseline_results.csv"), row.names = FALSE)
  write.csv(candidate, file.path(out_dir, "candidate_results.csv"), row.names = FALSE)
}

bootstrap_median_ci <- function(x, reps = 2000L) {
  if (!length(x)) {
    return(c(NA_real_, NA_real_))
  }
  boot <- replicate(reps, median(sample(x, replace = TRUE)))
  as.numeric(stats::quantile(boot, probs = c(0.025, 0.975), names = FALSE, type = 7))
}

compare_runs <- function(baseline_csv, candidate_csv, out_dir, mei_abs, mei_rel, mei_tail_abs, mei_tail_rel, spread_tol_rel) {
  baseline <- read.csv(baseline_csv, stringsAsFactors = FALSE)
  candidate <- read.csv(candidate_csv, stringsAsFactors = FALSE)

  key <- c("scenario", "replicate", "seed_mode", "data_seed", "optim_seed")
  merged <- merge(
    baseline,
    candidate,
    by = key,
    suffixes = c("_baseline", "_candidate"),
    all = FALSE,
    sort = TRUE
  )

  if (!nrow(merged)) {
    stop("No paired rows available after merging baseline and candidate results", call. = FALSE)
  }

  scenarios <- unique(merged$scenario)
  summaries <- vector("list", length(scenarios))
  pvals <- numeric(length(scenarios))

  for (i in seq_along(scenarios)) {
    scenario <- scenarios[i]
    rows <- merged[merged$scenario == scenario, , drop = FALSE]
    delta <- rows$elapsed_wall_candidate - rows$elapsed_wall_baseline
    base <- rows$elapsed_wall_baseline
    mean_delta <- mean(delta)
    median_delta <- median(delta)
    rel_mean <- 100 * mean_delta / mean(base)
    rel_median <- 100 * median_delta / median(base)
    if (length(delta) > 1L && is.finite(stats::sd(delta))) {
      mean_ci_half <- qt(0.975, df = length(delta) - 1L) * sd(delta) / sqrt(length(delta))
      mean_ci <- c(mean_delta - mean_ci_half, mean_delta + mean_ci_half)
    } else {
      mean_ci <- c(NA_real_, NA_real_)
    }
    median_ci <- bootstrap_median_ci(delta)
    pvals[i] <- if (length(delta) > 1L && sd(delta) > 0) {
      t.test(delta, alternative = "less", mu = 0)$p.value
    } else if (all(delta < 0)) {
      0
    } else {
      1
    }

    parity_degree <- identical(rows$degree_baseline, rows$degree_candidate)
    parity_conv <- identical(rows$convergence_baseline, rows$convergence_candidate)
    parity_hy <- max(abs(rows$h_y_candidate - rows$h_y_baseline))
    parity_hx <- max(abs(rows$h_x_candidate - rows$h_x_baseline))
    parity_value <- max(abs(rows$value_candidate - rows$value_baseline))

    summaries[[i]] <- data.frame(
      scenario = scenario,
      n_pairs = nrow(rows),
      baseline_mean = mean(base),
      candidate_mean = mean(rows$elapsed_wall_candidate),
      mean_delta = mean_delta,
      mean_delta_ci_lo = mean_ci[1],
      mean_delta_ci_hi = mean_ci[2],
      median_delta = median_delta,
      median_delta_ci_lo = median_ci[1],
      median_delta_ci_hi = median_ci[2],
      mean_pct = rel_mean,
      median_pct = rel_median,
      p90_delta = as.numeric(stats::quantile(delta, probs = 0.90, names = FALSE)),
      p95_delta = as.numeric(stats::quantile(delta, probs = 0.95, names = FALSE)),
      sd_delta = stats::sd(delta),
      iqr_delta = IQR(delta),
      mad_delta = mad(delta),
      parity_degree = parity_degree,
      parity_convergence = parity_conv,
      parity_max_abs_h_y = parity_hy,
      parity_max_abs_h_x = parity_hx,
      parity_max_abs_value = parity_value,
      p_value = pvals[i],
      keep_gate = (
        mean_delta <= -mei_abs &&
          median_delta <= -mei_abs &&
          rel_mean <= -mei_rel &&
          rel_median <= -mei_rel &&
          isTRUE(mean_ci[2] < 0) &&
          isTRUE(median_ci[2] < 0) &&
          as.numeric(stats::quantile(delta, probs = 0.95, names = FALSE)) <= mei_tail_abs &&
          (100 * as.numeric(stats::quantile(delta / base, probs = 0.95, names = FALSE))) <= mei_tail_rel &&
          parity_degree &&
          parity_conv &&
          parity_hy <= 1e-12 &&
          parity_hx <= 1e-12 &&
          parity_value <= 1e-10
      ),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summaries)
  summary_df$p_value_holm <- p.adjust(summary_df$p_value, method = "holm")
  summary_df$spread_ok <- with(
    summary_df,
    abs(mean_pct) >= mei_rel |
      (
        sd_delta <= spread_tol_rel * pmax(abs(mean_delta), mei_abs) &
          iqr_delta <= spread_tol_rel * pmax(abs(median_delta), mei_abs)
      )
  )
  summary_df$keep_gate <- summary_df$keep_gate & summary_df$spread_ok & summary_df$p_value_holm <= 0.05

  write.csv(merged, file.path(out_dir, "paired_rows.csv"), row.names = FALSE)
  write.csv(summary_df, file.path(out_dir, "summary.csv"), row.names = FALSE)

  lines <- c(
    sprintf("baseline_csv=%s", normalizePath(baseline_csv, winslash = "/", mustWork = TRUE)),
    sprintf("candidate_csv=%s", normalizePath(candidate_csv, winslash = "/", mustWork = TRUE)),
    "",
    "Declared MEI",
    sprintf("center_abs_seconds=%.3f", mei_abs),
    sprintf("center_rel_percent=%.1f", mei_rel),
    sprintf("tail_abs_seconds=%.3f", mei_tail_abs),
    sprintf("tail_rel_percent=%.1f", mei_tail_rel),
    sprintf("spread_tolerance_ratio=%.2f", spread_tol_rel),
    "",
    "Per-scenario summary"
  )
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    lines <- c(
      lines,
      sprintf(
        "%s mean_delta=%.4f median_delta=%.4f mean_pct=%.2f median_pct=%.2f p95_delta=%.4f p_holm=%.4g parity[h_y=%.3g h_x=%.3g value=%.3g degree=%s conv=%s] keep=%s",
        row$scenario,
        row$mean_delta,
        row$median_delta,
        row$mean_pct,
        row$median_pct,
        row$p95_delta,
        row$p_value_holm,
        row$parity_max_abs_h_y,
        row$parity_max_abs_h_x,
        row$parity_max_abs_value,
        row$parity_degree,
        row$parity_convergence,
        row$keep_gate
      )
    )
  }
  write_lines(lines, file.path(out_dir, "summary.txt"))
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  mode <- opts$mode %||% "run"

  if (identical(mode, "profile")) {
    out_dir <- ensure_dir(require_arg(opts, "out_dir"))
    run_profile(
      repo = require_arg(opts, "repo"),
      out_dir = out_dir,
      scenario_id = opts$scenario %||% "cv_ls_deg1_c1",
      n = as_int(opts$n, 180L),
      n.integrate = as_int(opts$n_integrate, 41L),
      nmulti = as_int(opts$nmulti, 2L)
    )
    return(invisible(NULL))
  }

  if (identical(mode, "objective-profile")) {
    out_dir <- ensure_dir(require_arg(opts, "out_dir"))
    run_objective_profile(
      repo = require_arg(opts, "repo"),
      out_dir = out_dir,
      scenario_id = opts$scenario %||% "cv_ml_proper_deg1_c1",
      n = as_int(opts$n, 180L),
      n.integrate = as_int(opts$n_integrate, 41L)
    )
    return(invisible(NULL))
  }

  if (identical(mode, "run")) {
    out_dir <- ensure_dir(require_arg(opts, "out_dir"))
    run_benchmarks(
      repo = require_arg(opts, "repo"),
      out_dir = out_dir,
      times = as_int(opts$times, 5L),
      seed_mode = opts$seed_mode %||% "fixed",
      n = as_int(opts$n, 180L),
      n.integrate = as_int(opts$n_integrate, 41L),
      nmulti = as_int(opts$nmulti, 2L),
      selected = if (is.null(opts$scenarios)) NULL else strsplit(opts$scenarios, ",", fixed = TRUE)[[1L]]
    )
    return(invisible(NULL))
  }

  if (identical(mode, "objective-run")) {
    out_dir <- ensure_dir(require_arg(opts, "out_dir"))
    run_objective_benchmarks(
      repo = require_arg(opts, "repo"),
      out_dir = out_dir,
      times = as_int(opts$times, 5L),
      seed_mode = opts$seed_mode %||% "fixed",
      n = as_int(opts$n, 180L),
      n.integrate = as_int(opts$n_integrate, 41L),
      selected = if (is.null(opts$scenarios)) NULL else strsplit(opts$scenarios, ",", fixed = TRUE)[[1L]]
    )
    return(invisible(NULL))
  }

  if (identical(mode, "single-fit")) {
    run_single_fit(
      lib_dir = require_arg(opts, "lib_dir"),
      out_csv = require_arg(opts, "out_csv"),
      scenario_id = require_arg(opts, "scenario"),
      data_seed = as_int(require_arg(opts, "data_seed")),
      optim_seed = as_int(require_arg(opts, "optim_seed")),
      n = as_int(opts$n, 180L),
      n.integrate = as_int(opts$n_integrate, 41L),
      nmulti = as_int(opts$nmulti, 2L)
    )
    return(invisible(NULL))
  }

  if (identical(mode, "single-objective")) {
    run_single_objective(
      lib_dir = require_arg(opts, "lib_dir"),
      out_csv = require_arg(opts, "out_csv"),
      scenario_id = require_arg(opts, "scenario"),
      data_seed = as_int(require_arg(opts, "data_seed")),
      n = as_int(opts$n, 180L),
      n.integrate = as_int(opts$n_integrate, 41L)
    )
    return(invisible(NULL))
  }

  if (identical(mode, "interleaved-compare")) {
    out_dir <- ensure_dir(require_arg(opts, "out_dir"))
    run_interleaved_compare(
      baseline_repo = require_arg(opts, "baseline_repo"),
      candidate_repo = require_arg(opts, "candidate_repo"),
      out_dir = out_dir,
      times = as_int(opts$times, 5L),
      seed_mode = opts$seed_mode %||% "fixed",
      n = as_int(opts$n, 180L),
      n.integrate = as_int(opts$n_integrate, 41L),
      nmulti = as_int(opts$nmulti, 2L),
      selected = if (is.null(opts$scenarios)) NULL else strsplit(opts$scenarios, ",", fixed = TRUE)[[1L]]
    )
    return(invisible(NULL))
  }

  if (identical(mode, "objective-interleaved-compare")) {
    out_dir <- ensure_dir(require_arg(opts, "out_dir"))
    run_objective_interleaved_compare(
      baseline_repo = require_arg(opts, "baseline_repo"),
      candidate_repo = require_arg(opts, "candidate_repo"),
      out_dir = out_dir,
      times = as_int(opts$times, 5L),
      seed_mode = opts$seed_mode %||% "fixed",
      n = as_int(opts$n, 180L),
      n.integrate = as_int(opts$n_integrate, 41L),
      selected = if (is.null(opts$scenarios)) NULL else strsplit(opts$scenarios, ",", fixed = TRUE)[[1L]]
    )
    return(invisible(NULL))
  }

  if (identical(mode, "compare")) {
    out_dir <- ensure_dir(require_arg(opts, "out_dir"))
    compare_runs(
      baseline_csv = require_arg(opts, "baseline_csv"),
      candidate_csv = require_arg(opts, "candidate_csv"),
      out_dir = out_dir,
      mei_abs = as_num(opts$mei_abs, 0.05),
      mei_rel = as_num(opts$mei_rel, 5.0),
      mei_tail_abs = as_num(opts$mei_tail_abs, 0.03),
      mei_tail_rel = as_num(opts$mei_tail_rel, 3.0),
      spread_tol_rel = as_num(opts$spread_tol_rel, 4.0)
    )
    return(invisible(NULL))
  }

  stop(sprintf("Unsupported --mode value: %s", mode), call. = FALSE)
}

main()
