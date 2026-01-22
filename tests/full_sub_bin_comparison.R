################################################################################
## Bivariate Boundary Binning Performance Test
## Jeff Racine - Feasible Resampling vs. Binned vs. Unbinned
################################################################################

## User-defined parameters
M <- 100                     # Number of replications
n_vec <- c(200, 400, 800, 1600) 
num_bins <- 32              
n_sub_default <- 100        # Fixed subsample size for feasible resampling
eval_grid_size <- 20        
plot_enabled <- TRUE        

library(bkcde)

# DGP: Beta distribution with parameters depending on X to create complex boundary
dgp_beta <- function(n) {
  x <- runif(n)
  shape1 <- 2 + 5 * x
  shape2 <- 5 - 2 * x
  y <- rbeta(n, shape1, shape2)
  return(data.frame(x = x, y = y))
}

true_density <- function(x_grid, y_grid) {
  s1 <- 2 + 5 * x_grid
  s2 <- 5 - 2 * x_grid
  dbeta(y_grid, s1, s2)
}

methods_list <- c("Unbinned", "Binned", "Subsample")

# Data frames to store simulation results
results_time <- expand.grid(M = 1:M, n = n_vec, method = methods_list, stringsAsFactors = FALSE)
results_time$seconds <- NA 

results_mse <- expand.grid(M = 1:M, n = n_vec, method = methods_list, stringsAsFactors = FALSE)
results_mse$mse <- NA

bw_storage <- expand.grid(M = 1:M, n = n_vec, method = methods_list, stringsAsFactors = FALSE)
bw_storage$hx <- NA
bw_storage$hy <- NA

for (m in 1:M) {
  for (n_curr in n_vec) {
    cat(sprintf("\rReplication %d of %d | n = %-4d", m, M, n_curr))
    flush.console()
    
    # 1. Generate Data and Evaluation Grid
    data <- dgp_beta(n_curr)
    x.lb <- min(data$x); x.ub <- max(data$x)
    y.lb <- min(data$y); y.ub <- max(data$y)
    
    x_seq <- seq(x.lb, x.ub, length.out = eval_grid_size)
    y_seq <- seq(y.lb, y.ub, length.out = eval_grid_size)
    eval_grid <- expand.grid(x = x_seq, y = y_seq)
    
    x.eval <- as.numeric(eval_grid$x)
    y.eval <- as.numeric(eval_grid$y)
    f_true <- true_density(x.eval, y.eval)
    
    # --- METHOD 1: Unbinned Full CV ---
    t1 <- proc.time()
    fit_unb <- bkcde(x = data$x, y = data$y, x.eval = x.eval, y.eval = y.eval, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "full", cv.binned = FALSE)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Unbinned", "seconds"] <- (proc.time() - t1)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Unbinned", "mse"] <- mean((fit_unb$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Unbinned", "hx"] <- fit_unb$h[1]
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Unbinned", "hy"] <- fit_unb$h[2]
    
    # --- METHOD 2: Binned Full CV ---
    t2 <- proc.time()
    fit_bin <- bkcde(x = data$x, y = data$y, x.eval = x.eval, y.eval = y.eval, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "full", cv.binned = TRUE, nbins = num_bins)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Binned", "seconds"] <- (proc.time() - t2)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Binned", "mse"] <- mean((fit_bin$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Binned", "hx"] <- fit_bin$h[1]
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Binned", "hy"] <- fit_bin$h[2]
    
    # --- METHOD 3: Feasible Resampling (Subsample CV) ---
    # [cite_start]Resamples rows from the data to avoid O(n^2k) recomputation[cite: 1, 14, 94].
    t3 <- proc.time()
    fit_sub <- bkcde(x = data$x, y = data$y, x.eval = x.eval, y.eval = y.eval, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "sub", n.sub = n_sub_default, cv.binned = FALSE)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Subsample", "seconds"] <- (proc.time() - t3)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Subsample", "mse"] <- mean((fit_sub$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Subsample", "hx"] <- fit_sub$h[1]
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Subsample", "hy"] <- fit_sub$h[2]
    
    # Live Plotting
    if (plot_enabled && (m > 1 || n_curr >= max(n_vec))) {
      dev.hold()
      df_mse_sub <- results_mse[!is.na(results_mse$mse), ]
      df_mse_sub$method <- factor(df_mse_sub$method, levels = methods_list)
      df_time_sub <- results_time[!is.na(results_time$seconds), ]
      df_time_sub$method <- factor(df_time_sub$method, levels = methods_list)
      
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
      
      # 1. SCALING PLOT: Median Lines + Jittered Raw Points
      cols <- c("black", "blue", "red")
      medians_plot <- sapply(methods_list, function(meth) {
        sapply(n_vec, function(nv) median(df_time_sub$seconds[df_time_sub$method == meth & df_time_sub$n == nv], na.rm=TRUE))
      })
      
      matplot(n_vec, medians_plot, type = "l", col = cols, lty = 1, lwd = 2,
              xlab = "n", ylab = "Time (s)", main = paste("Scaling (M =", m, ")"))
      
      for(i in seq_along(methods_list)) {
        sub_data <- df_time_sub[df_time_sub$method == methods_list[i], ]
        points(jitter(sub_data$n, amount = 15), sub_data$seconds, col = cols[i], pch = 20, cex = 0.5)
      }
      legend("topleft", legend = methods_list, col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.8)
      
      # 2. MSE at Smallest n (Added outline=FALSE, notch=TRUE)
      d_small <- df_mse_sub[abs(df_mse_sub$n - min(n_vec)) < 1e-7, ]
      if(nrow(d_small) > 0) {
        boxplot(mse ~ method, data = d_small, 
                main = paste("MSE (n =", min(n_vec), ")"), 
                col = "lightgray", outline = FALSE, notch = TRUE)
      } else plot.new()
      
      # 3. MSE at Largest n (Added outline=FALSE, notch=TRUE)
      d_large <- df_mse_sub[abs(df_mse_sub$n - max(n_vec)) < 1e-7, ]
      if(nrow(d_large) > 0) {
        boxplot(mse ~ method, data = d_large, 
                main = paste("MSE (n =", max(n_vec), ")"), 
                col = "lightblue", outline = FALSE, notch = TRUE)
      } else plot.new()
      
      # 4. BANDWIDTHS: Compare h.x and h.y for ALL methods (Added outline=FALSE, notch=TRUE)
      bw_subset <- bw_storage[abs(bw_storage$n - max(n_vec)) < 1e-7 & !is.na(bw_storage$hx), ]
      if(nrow(bw_subset) > 0) {
        bw_list <- list(
          "U.hx" = bw_subset$hx[bw_subset$method == "Unbinned"],
          "B.hx" = bw_subset$hx[bw_subset$method == "Binned"],
          "S.hx" = bw_subset$hx[bw_subset$method == "Subsample"],
          "U.hy" = bw_subset$hy[bw_subset$method == "Unbinned"],
          "B.hy" = bw_subset$hy[bw_subset$method == "Binned"],
          "S.hy" = bw_subset$hy[bw_subset$method == "Subsample"]
        )
        bw_cols <- rep(c("black", "blue", "red"), 2)
        boxplot(bw_list, main = "Bandwidths (Max n)", 
                col = bw_cols, las = 2, cex.axis = 0.7, 
                outline = FALSE, notch = TRUE)
        abline(v = 3.5, lty = 2, col = "gray") 
      } else {
        plot.new()
      }
      
      dev.flush()
    }
  }
}

cat("\n\nSimulation Complete.\n")

## Final Summary and Export
final_stats <- do.call(rbind, lapply(n_vec, function(n_val) {
  idx_n <- abs(results_time$n - n_val) < 1e-7
  t_unb <- results_time$seconds[idx_n & results_time$method == "Unbinned"]
  t_bin <- results_time$seconds[idx_n & results_time$method == "Binned"]
  speedup <- t_unb / t_bin
  mse_bin <- results_mse$mse[abs(results_mse$n - n_val) < 1e-7 & results_mse$method == "Binned"]
  data.frame(n = n_val, 
             Avg_Speedup = mean(speedup, na.rm=TRUE), 
             Med_Speedup = median(speedup, na.rm=TRUE),
             Avg_MSE_Binned = mean(mse_bin, na.rm=TRUE))
}))

print(final_stats, row.names = FALSE)
write.csv(final_stats, "benchmark_summary_table.csv", row.names = FALSE)