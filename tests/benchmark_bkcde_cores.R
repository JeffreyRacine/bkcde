suppressPackageStartupMessages({
  library(bkcde)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

## ============================================================================
## Configuration: adjust these to suit your machine and exploration needs
## ============================================================================

# Data generation
set.seed(42)
n <- 1000

# Core sweep parameters for optimization
# Note: Going beyond detectCores() can help if tasks finish at different rates
# (e.g., lower-degree polynomials finish faster, allowing idle cores to pick up work)
optim.degree.cores.min <- 1
optim.degree.cores.max <- 1.5*detectCores()
optim.degree.cores.by  <- 2

optim.nmulti.cores.min <- 1
optim.nmulti.cores.max <- 1.5*detectCores()
optim.nmulti.cores.by  <- 2

optim.ksum.cores.min <- 1
optim.ksum.cores.max <- 2  # Keep modest; ksum parallelism often has overhead
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
      cat("FINAL RECOMMENDATIONS\n")
      cat(strrep("=", 70), "\n", sep = "")
      cat("\nOptimization (cv.only) parameters:\n")
      cat("  optim.degree.cores =", best$optim.degree.cores, "\n")
      cat("  optim.nmulti.cores =", best$optim.nmulti.cores, "\n")
      cat("  optim.ksum.cores   =", best$optim.ksum.cores, "\n")
      cat("\nFitting parameters:\n")
      cat("  fitted.cores =", best_fitted$fitted.cores, "\n")
      cat("  proper.cores =", best_proper$proper.cores, "\n")
      cat("\nSuggested bkcde() call:\n")
      cat("  f.yx <- bkcde(x=x, y=y, bwmethod='cv.ml',\n")
      cat("               optim.cores='manual',\n")
      cat("               optim.degree.cores=", best$optim.degree.cores, ",\n", sep = "")
      cat("               optim.nmulti.cores=", best$optim.nmulti.cores, ",\n", sep = "")
      cat("               optim.ksum.cores=", best$optim.ksum.cores, ",\n", sep = "")
      cat("               fitted.cores=", best_fitted$fitted.cores, ",\n", sep = "")
      cat("               proper=TRUE, proper.cores=", best_proper$proper.cores, ")\n", sep = "")
      cat("\n", strrep("=", 70), "\n", sep = "")
    }
    
    write.csv(fitted_results, "benchmark_fitted_results.csv", row.names = FALSE)
    write.csv(proper_results, "benchmark_proper_results.csv", row.names = FALSE)
  }
}

# Optional: simple heat map of optimization elapsed time by degree/nmulti, faceted by ksum
p <- ggplot(results, aes(x = factor(optim.degree.cores),
                        y = factor(optim.nmulti.cores),
                        fill = optim_elapsed)) +
  geom_tile(color = "white") +
  facet_wrap(~ optim.ksum.cores, labeller = label_both) +
  scale_fill_viridis_c(name = "Elapsed (s)") +
  labs(x = "optim.degree.cores", y = "optim.nmulti.cores",
       title = "bkcde cv-only elapsed time by core allocation")

ggsave("benchmark_heatmap.png", p, width = 8, height = 4, dpi = 150)
cat("\nSaved results to benchmark_results.csv and plot to benchmark_heatmap.png\n")
