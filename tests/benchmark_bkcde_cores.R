suppressPackageStartupMessages({
  library(bkcde)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
})

set.seed(42)

n <- 500
x <- runif(n, 0, 1)
y <- rnorm(n, mean = 2 * sin(4 * pi * x), sd = 1 + abs(x))

# Adjust the grids to suit your machine.
grid <- expand_grid(
  optim.degree.cores = c(1, 2, 4, 8),
  optim.nmulti.cores = c(1, 2, 4, 8),
  optim.ksum.cores   = c(1, 2)
)

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
  
  fitted_grid <- tibble(fitted.cores = c(1, 2, 4, 8))
  
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
  }
  
  write.csv(fitted_results, "benchmark_fitted_results.csv", row.names = FALSE)
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
