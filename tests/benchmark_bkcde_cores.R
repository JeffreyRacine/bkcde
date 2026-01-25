suppressPackageStartupMessages({
  library(bkcde)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(parallel)
})

## ============================================================================
## Configuration: adjust these to suit your machine and exploration needs
## ============================================================================

# Data generation
set.seed(42)
n <- 250

# Core sweep parameters for optimization
# Note: Going beyond detectCores() can help if tasks finish at different rates
# (e.g., lower-degree polynomials finish faster, allowing idle cores to pick up work)
optim.degree.cores.min <- 1
optim.degree.cores.max <- 2
optim.degree.cores.by  <- 1

optim.nmulti.cores.min <- 1
optim.nmulti.cores.max <- 2
optim.nmulti.cores.by  <- 1

optim.ksum.cores.min <- 1
optim.ksum.cores.max <- 2
optim.ksum.cores.by  <- 1

# Core sweep parameters for fitting
fitted.cores.min <- 1
fitted.cores.max <- min(6, detectCores())  # Fitting often benefits less from parallelism
fitted.cores.by  <- 1

# Core sweep parameters for proper (ensuring proper density)
proper.cores.min <- 1
proper.cores.max <- min(6, detectCores())  # Proper enforcement also benefits less from parallelism
proper.cores.by  <- 1

## ============================================================================
## Data generation
## ============================================================================

x <- runif(n, 0, 1)
y <- rnorm(n, mean = 2 * sin(4 * pi * x), sd = 1 + abs(x))

## ============================================================================
## Build optimization grid
## ============================================================================

grid <- expand_grid(
  optim.degree.cores = seq(optim.degree.cores.min, optim.degree.cores.max, by = optim.degree.cores.by),
  optim.nmulti.cores = seq(optim.nmulti.cores.min, optim.nmulti.cores.max, by = optim.nmulti.cores.by),
  optim.ksum.cores   = seq(optim.ksum.cores.min, optim.ksum.cores.max, by = optim.ksum.cores.by)
)

cat("Detected cores:", detectCores(), "\n")
cat("Optimization grid size:", nrow(grid), "configurations\n")
cat("Core ranges: degree=[", min(grid$optim.degree.cores), ",", max(grid$optim.degree.cores),
    "], nmulti=[", min(grid$optim.nmulti.cores), ",", max(grid$optim.nmulti.cores),
    "], ksum=[", min(grid$optim.ksum.cores), ",", max(grid$optim.ksum.cores), "]\n\n", sep = "")

## ============================================================================
## Benchmark optimization (cv.only)
## ============================================================================

run_one <- function(optim.degree.cores, optim.nmulti.cores, optim.ksum.cores) {
  cat(sprintf("Running deg=%d nmulti=%d ksum=%d ...\n",
              optim.degree.cores, optim.nmulti.cores, optim.ksum.cores))
  tryCatch({
    fit <- bkcde(
      x = x, y = y,
      bwmethod = "cv.ml",
      cv.only = TRUE,                # benchmark optimization only
      optim.cores = "manual",       # honor manual settings
      optim.degree.cores = optim.degree.cores,
      optim.nmulti.cores = optim.nmulti.cores,
      optim.ksum.cores = optim.ksum.cores,
      fitted.cores = 1,
      proper.cores = 1
    )

    # Compute basic metrics from the returned object (summary.bkcde only prints).
    optim_cpu <- sum(fit$secs.optim.mat, na.rm = TRUE)
    optim_elapsed <- fit$secs.optim.elapsed
    opt_cores <- max(1, fit$optim.ksum.cores * fit$optim.degree.cores * fit$optim.nmulti.cores)
    speedup_opt <- optim_cpu / max(1e-9, optim_elapsed)
    eff_opt <- speedup_opt / opt_cores

    tibble(
      optim.degree.cores = optim.degree.cores,
      optim.nmulti.cores = optim.nmulti.cores,
      optim.ksum.cores = optim.ksum.cores,
      elapsed_total = fit$secs.elapsed,
      optim_cpu = optim_cpu,
      optim_elapsed = optim_elapsed,
      speedup_opt = speedup_opt,
      eff_opt = eff_opt,
      convergence = fit$convergence,
      error = NA_character_
    )
  }, error = function(e) {
    tibble(
      optim.degree.cores = optim.degree.cores,
      optim.nmulti.cores = optim.nmulti.cores,
      optim.ksum.cores = optim.ksum.cores,
      elapsed_total = NA_real_,
      optim_cpu = NA_real_,
      optim_elapsed = NA_real_,
      speedup_opt = NA_real_,
      eff_opt = NA_real_,
      convergence = NA_integer_,
      error = conditionMessage(e)
    )
  })
}

results <- pmap_dfr(grid, run_one)

print(results)

write.csv(results, "benchmark_results.csv", row.names = FALSE)

# Identify the fastest valid configuration (lowest optim_elapsed) with good convergence
best <- results %>%
  filter(is.na(error) | error == "",
         !is.na(optim_elapsed),
         !is.na(convergence),
         convergence == 0) %>%
  arrange(optim_elapsed) %>%
  slice_head(n = 1)

cat("\nBest (lowest optimization elapsed) configuration:\n")
print(best)

if (nrow(best) == 1) {
  cat("\nSuggested settings: optim.cores='manual', optim.degree.cores=",
      best$optim.degree.cores, ", optim.nmulti.cores=", best$optim.nmulti.cores,
      ", optim.ksum.cores=", best$optim.ksum.cores, "\n", sep = "")
  
  # Now benchmark fitting with optimal h and degree from cv.only run
  cat("\n=== Running fitting benchmark with optimal h and degree ===\n")
  
  # Get optimal h and degree from a cv.only run with best config
  fit_opt <- bkcde(
    x = x, y = y,
    bwmethod = "cv.ml",
    cv.only = TRUE,
    optim.cores = "manual",
    optim.degree.cores = best$optim.degree.cores,
    optim.nmulti.cores = best$optim.nmulti.cores,
    optim.ksum.cores = best$optim.ksum.cores,
    fitted.cores = 1,
    proper.cores = 1
  )
  
  cat("Optimal h.y =", fit_opt$h[1], ", h.x =", fit_opt$h[2], ", degree =", fit_opt$degree, "\n")
  
  fitted_grid <- tibble(
    fitted.cores = seq(fitted.cores.min, fitted.cores.max, by = fitted.cores.by)
  )
  
  run_fitted <- function(fitted.cores) {
    cat(sprintf("  Testing fitted.cores=%d ...\n", fitted.cores))
    tryCatch({
      fit <- bkcde(
        x = x, y = y,
        h = fit_opt$h,
        degree = fit_opt$degree,
        fitted.cores = fitted.cores,
        proper.cores = 1
      )
      
      tibble(
        fitted.cores = fitted.cores,
        elapsed_total = fit$secs.elapsed,
        elapsed_fit = fit$secs.estimate,
        error = NA_character_
      )
    }, error = function(e) {
      tibble(
        fitted.cores = fitted.cores,
        elapsed_total = NA_real_,
        elapsed_fit = NA_real_,
        error = conditionMessage(e)
      )
    })
  }
  
  fitted_results <- pmap_dfr(fitted_grid, run_fitted)
  
  print(fitted_results)
  
  # Pick best: lowest elapsed_fit; in ties prefer fewest cores
  best_fitted <- fitted_results %>%
    filter(is.na(error) | error == "", !is.na(elapsed_fit)) %>%
    arrange(elapsed_fit, fitted.cores) %>%
    slice_head(n = 1)
  
  cat("\nBest fitted.cores (lowest fitting elapsed, fewest cores in ties):\n")
  print(best_fitted)
  
  if (nrow(best_fitted) == 1) {
    cat("\nSuggested fitted.cores =", best_fitted$fitted.cores, "\n")
    
    # Now benchmark proper.cores with proper=TRUE
    cat("\n=== Running proper.cores benchmark with proper=TRUE ===\n")
    
    proper_grid <- tibble(
      proper.cores = seq(proper.cores.min, proper.cores.max, by = proper.cores.by)
    )
    
    run_proper <- function(proper.cores) {
      cat(sprintf("  Testing proper.cores=%d ...\n", proper.cores))
      tryCatch({
        fit <- bkcde(
          x = x, y = y,
          h = fit_opt$h,
          degree = fit_opt$degree,
          proper = TRUE,
          fitted.cores = best_fitted$fitted.cores,
          proper.cores = proper.cores
        )
        
        tibble(
          proper.cores = proper.cores,
          elapsed_total = fit$secs.elapsed,
          elapsed_proper = fit$secs.elapsed - fit$secs.estimate,  # Approximate proper time
          error = NA_character_
        )
      }, error = function(e) {
        tibble(
          proper.cores = proper.cores,
          elapsed_total = NA_real_,
          elapsed_proper = NA_real_,
          error = conditionMessage(e)
        )
      })
    }
    
    proper_results <- pmap_dfr(proper_grid, run_proper)
    
    print(proper_results)
    
    # Pick best: lowest elapsed_total; in ties prefer fewest cores
    best_proper <- proper_results %>%
      filter(is.na(error) | error == "", !is.na(elapsed_total)) %>%
      arrange(elapsed_total, proper.cores) %>%
      slice_head(n = 1)
    
    cat("\nBest proper.cores (lowest total elapsed, fewest cores in ties):\n")
    print(best_proper)
    
    if (nrow(best_proper) == 1) {
      cat("\nSuggested proper.cores =", best_proper$proper.cores, "\n")
      
      # Final summary: all five parameters
      cat("\n", strrep("=", 70), "\n", sep = "")
      cat("FINAL RECOMMENDATIONS (TUNED)\n")
      cat(strrep("=", 70), "\n", sep = "")
      cat("\nOptimization (cv.only) parameters:\n")
      cat("  optim.degree.cores =", best$optim.degree.cores, "\n")
      cat("  optim.nmulti.cores =", best$optim.nmulti.cores, "\n")
      cat("  optim.ksum.cores   =", best$optim.ksum.cores, "\n")
      cat("\nFitting parameters:\n")
      cat("  fitted.cores =", best_fitted$fitted.cores, "\n")
      cat("  proper.cores =", best_proper$proper.cores, "\n")
      cat("\nSuggested bkcde() call (tuned):\n")
      cat("  f.yx <- bkcde(x=x, y=y, bwmethod='cv.ml',\n")
      cat("               optim.cores='manual',\n")
      cat("               optim.degree.cores=", best$optim.degree.cores, ",\n", sep = "")
      cat("               optim.nmulti.cores=", best$optim.nmulti.cores, ",\n", sep = "")
      cat("               optim.ksum.cores=", best$optim.ksum.cores, ",\n", sep = "")
      cat("               fitted.cores=", best_fitted$fitted.cores, ",\n", sep = "")
      cat("               proper=TRUE, proper.cores=", best_proper$proper.cores, ")\n", sep = "")
      cat("\n", strrep("=", 70), "\n", sep = "")
      
      # Now run with default options (optim.cores='auto')
      cat("\n=== Running with DEFAULT options (optim.cores='auto') ===\n")
      
      fit_default <- tryCatch({
        bkcde(
          x = x, y = y,
          bwmethod = "cv.ml",
          proper = TRUE
          # All other params use defaults, optim.cores='auto'
        )
      }, error = function(e) {
        cat("Error running default configuration:", conditionMessage(e), "\n")
        NULL
      })
      
      if (!is.null(fit_default)) {
        cat("\nDefault configuration completed.\n")
        cat("Optimal h.y =", fit_default$h[1], ", h.x =", fit_default$h[2], 
            ", degree =", fit_default$degree, "\n")
        
        # Compare timings: tuned should include optimization + fitting + proper
        # Best optimization time (cv.only stage)
        best_optim_elapsed <- best$optim_elapsed
        # Best fitting time (already included in proper stage)
        best_fit_elapsed <- best_fitted$elapsed_fit
        # Best proper time (includes fitting)
        best_proper_elapsed <- best_proper$elapsed_total
        
        # Total tuned: optimization stage + proper stage (proper includes fitting)
        tuned_total <- best_optim_elapsed + best_proper_elapsed
        default_total <- fit_default$secs.elapsed
        speedup_factor <- default_total / tuned_total
        time_saved <- default_total - tuned_total
        
        cat("\n", strrep("=", 70), "\n", sep = "")
        cat("PERFORMANCE COMPARISON (FULL WORKFLOW)\n")
        cat(strrep("=", 70), "\n", sep = "")
        cat("\nTuned configuration breakdown:\n")
        cat("  Optimization (cv.only):  ", formatC(best_optim_elapsed, format="f", digits=4), " sec\n", sep = "")
        cat("  Fitting + Proper:        ", formatC(best_proper_elapsed, format="f", digits=4), " sec\n", sep = "")
        cat("  ────────────────────────────────\n")
        cat("  TOTAL (tuned):           ", formatC(tuned_total, format="f", digits=4), " sec\n\n", sep = "")
        cat("Default (auto) configuration: ", formatC(default_total, format="f", digits=4), " sec\n", sep = "")
        
        if (speedup_factor > 1) {
          cat("\nResult: Tuned is FASTER\n")
          cat("  ", formatC(speedup_factor, format="f", digits=2), "x faster\n", sep = "")
          cat("  Time saved: ", formatC(time_saved, format="f", digits=4), " sec\n", sep = "")
        } else if (speedup_factor < 1) {
          cat("\nResult: Tuned is SLOWER\n")
          cat("  ", formatC(1/speedup_factor, format="f", digits=2), "x slower\n", sep = "")
          cat("  Time cost: ", formatC(-time_saved, format="f", digits=4), " sec\n", sep = "")
        } else {
          cat("\nResult: Equivalent performance\n")
        }
        cat("\n", strrep("=", 70), "\n", sep = "")
        
        cat("\n", strrep("=", 70), "\n", sep = "")
        cat("PARAMETER COMPARISON\n")
        cat(strrep("=", 70), "\n", sep = "")
        cat("\nTuned configuration:\n")
        cat("  degree =", fit_opt$degree, "\n")
        cat("  h.y    =", formatC(fit_opt$h[1], format="f", digits=8), "\n")
        cat("  h.x    =", formatC(fit_opt$h[2], format="f", digits=8), "\n")
        cat("\nDefault (auto) configuration:\n")
        cat("  degree =", fit_default$degree, "\n")
        cat("  h.y    =", formatC(fit_default$h[1], format="f", digits=8), "\n")
        cat("  h.x    =", formatC(fit_default$h[2], format="f", digits=8), "\n")
        
        # Check if parameters match
        params_match <- (fit_opt$degree == fit_default$degree &&
                         abs(fit_opt$h[1] - fit_default$h[1]) < 1e-6 &&
                         abs(fit_opt$h[2] - fit_default$h[2]) < 1e-6)
        
        if (params_match) {
          cat("\nNote: Both configurations converged to THE SAME parameters.\n")
        } else {
          cat("\nNote: Configurations converged to DIFFERENT parameters.\n")
        }
        cat("\n", strrep("=", 70), "\n", sep = "")
      }
    }
    
    write.csv(fitted_results, "benchmark_fitted_results.csv", row.names = FALSE)
    write.csv(proper_results, "benchmark_proper_results.csv", row.names = FALSE)
  }
}
