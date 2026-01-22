################################################################################
## Bivariate Boundary Binning Performance Test
## Jeff Racine
################################################################################

## User-defined parameters
M <- 10                     
n_vec <- c(200, 400, 800, 1600) 
num_bins <- 32              
n_sub_default <- 100        
eval_grid_size <- 20        
plot_enabled <- TRUE        

library(bkcde)

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
    
    data <- dgp_beta(n_curr)
    x.lb <- min(data$x); x.ub <- max(data$x)
    y.lb <- min(data$y); y.ub <- max(data$y)
    
    x_seq <- seq(x.lb, x.ub, length.out = eval_grid_size)
    y_seq <- seq(y.lb, y.ub, length.out = eval_grid_size)
    eval_grid <- expand.grid(x = x_seq, y = y_seq)
    
    x.eval <- as.numeric(eval_grid$x)
    y.eval <- as.numeric(eval_grid$y)
    f_true <- true_density(x.eval, y.eval)
    
    # 1. Unbinned Full
    t1 <- proc.time()
    fit_unb <- bkcde(x = data$x, y = data$y, x.eval = x.eval, y.eval = y.eval, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "full", cv.binned = FALSE)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Unbinned", "seconds"] <- (proc.time() - t1)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Unbinned", "mse"] <- mean((fit_unb$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Unbinned", "hx"] <- fit_unb$h[1]
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Unbinned", "hy"] <- fit_unb$h[2]
    
    # 2. Binned Full
    t2 <- proc.time()
    fit_bin <- bkcde(x = data$x, y = data$y, x.eval = x.eval, y.eval = y.eval, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "full", cv.binned = TRUE, nbins = num_bins)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Binned", "seconds"] <- (proc.time() - t2)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Binned", "mse"] <- mean((fit_bin$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Binned", "hx"] <- fit_bin$h[1]
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Binned", "hy"] <- fit_bin$h[2]
    
    # 3. Subsample Estimator
    t3 <- proc.time()
    fit_sub <- bkcde(x = data$x, y = data$y, x.eval = x.eval, y.eval = y.eval, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "sub", n.sub = n_sub_default, cv.binned = FALSE)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Subsample", "seconds"] <- (proc.time() - t3)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Subsample", "mse"] <- mean((fit_sub$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Subsample", "hx"] <- fit_sub$h[1]
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Subsample", "hy"] <- fit_sub$h[2]
    
    if (plot_enabled && (m > 1 || n_curr >= max(n_vec))) {
      dev.hold()
      df_mse_sub <- results_mse[!is.na(results_mse$mse), ]
      df_mse_sub$method <- factor(df_mse_sub$method, levels = methods_list)
      df_time_sub <- results_time[!is.na(results_time$seconds), ]
      df_time_sub$method <- factor(df_time_sub$method, levels = methods_list)
      
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
      
      # 1. SCALING PLOT with Raw Points and Median Lines
      cols <- c("black", "blue", "red")
      medians_plot <- sapply(methods_list, function(meth) {
        sapply(n_vec, function(nv) median(df_time_sub$seconds[df_time_sub$method == meth & df_time_sub$n == nv], na.rm=TRUE))
      })
      
      # Initialize plot with median trends
      matplot(n_vec, medians_plot, type = "l", col = cols, lty = 1, lwd = 2,
              xlab = "n", ylab = "Time (s)", main = paste("Scaling (M =", m, ")"))
      
      # Add raw data points for each method
      for(i in seq_along(methods_list)) {
        sub_data <- df_time_sub[df_time_sub$method == methods_list[i], ]
        # Add slight jitter to n for visibility if points overlap
        points(jitter(sub_data$n, amount = 10), sub_data$seconds, col = cols[i], pch = 20, cex = 0.6)
      }
      legend("topleft", legend = methods_list, col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.8)
      
      # 2. MSE Smallest n
      d_small <- df_mse_sub[abs(df_mse_sub$n - min(n_vec)) < 1e-7, ]
      if(nrow(d_small) > 0) boxplot(mse ~ method, data = d_small, main = paste("MSE (n =", min(n_vec), ")"), col = "lightgray") else plot.new()
      
      # 3. MSE Largest n
      d_large <- df_mse_sub[abs(df_mse_sub$n - max(n_vec)) < 1e-7, ]
      if(nrow(d_large) > 0) boxplot(mse ~ method, data = d_large, main = paste("MSE (n =", max(n_vec), ")"), col = "lightblue") else plot.new()
      
      # 4. Bandwidths
      bw_subset_plot <- bw_storage[abs(bw_storage$n - max(n_vec)) < 1e-7 & bw_storage$method == "Binned" & !is.na(bw_storage$hx), ]
      if(nrow(bw_subset_plot) > 0) boxplot(list(hx = bw_subset_plot$hx, hy = bw_subset_plot$hy), main = "Binned h (Max n)", col = c("orange", "green"), names = c("h.x", "h.y")) else plot.new()
      
      dev.flush()
    }
  }
}
cat("\n\nSimulation Complete.\n")