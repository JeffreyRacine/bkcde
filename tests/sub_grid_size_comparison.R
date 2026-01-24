################################################################################
## Bivariate Boundary Binning: Subsample Size Comparison
## Comparison of n.sub = 50, 100, and 200
################################################################################

## User-defined parameters
M <- 100                     
n_vec <- c(500, 1000, 1500) 
n_sub_vec <- c(100, 150, 200) # The three subsample sizes to test
eval_grid_size <- 25
plot_enabled <- TRUE        

library(bkcde)

# DGP: Beta distribution with parameters depending on X
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

# Define names based on the subsample sizes
methods_list <- paste0("Sub", n_sub_vec)

# Results Storage
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
    f_true <- true_density(eval_grid$x, eval_grid$y)
    
    # Loop through the three subsample sizes
    for (i in seq_along(n_sub_vec)) {
      sub_size <- n_sub_vec[i]
      m_name <- methods_list[i]
      
      t_start <- proc.time()
      # Use cv="sub" and varying n.sub values
      fit <- bkcde(x = data$x, y = data$y, 
                   x.eval = eval_grid$x, y.eval = eval_grid$y, 
                   x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                   cv = "sub", n.sub = sub_size, display.warnings = FALSE)
      
      results_time[results_time$n == n_curr & results_time$M == m & results_time$method == m_name, "seconds"] <- (proc.time() - t_start)[3]
      results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == m_name, "mse"] <- mean((fit$f - f_true)^2)
      bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == m_name, c("hy","hx")] <- fit$h[1:2]
    }
    
    if (plot_enabled && (m > 2)) {
      dev.hold()
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
      
      # 1. Scaling Plot (Time vs Subsample Size)
      cols <- c("darkgreen", "blue", "red")
      df_t <- results_time[!is.na(results_time$seconds), ]
      medians <- sapply(methods_list, function(mt) sapply(n_vec, function(nv) median(df_t$seconds[df_t$method==mt & df_t$n==nv], na.rm=T)))
      
      matplot(n_vec, medians, type="l", col=cols, lty=1, lwd=2, 
              xlab="n", ylab="Seconds", main="Scaling by Subsample Size")
      
      # Add Points vertically (no jitter)
      for(i in seq_along(methods_list)) {
        sub_data <- df_t[df_t$method == methods_list[i], ]
        points(sub_data$n, sub_data$seconds, col = cols[i], pch = 20, cex = 0.4)
      }
      legend("topleft", methods_list, col=cols, lty=1, bty="n", cex=0.8)
      
      # 2. MSE Smallest n
      d_s <- results_mse[results_mse$n == min(n_vec) & !is.na(results_mse$mse), ]
      boxplot(mse ~ method, data=d_s, main=paste("MSE (n =", min(n_vec), ")"), 
              outline=FALSE, notch=TRUE, col="lightgray")
      
      # 3. MSE Largest n
      d_l <- results_mse[results_mse$n == max(n_vec) & !is.na(results_mse$mse), ]
      boxplot(mse ~ method, data=d_l, main=paste("MSE (n =", max(n_vec), ")"), 
              outline=FALSE, notch=TRUE, col="lightblue")
      
      # 4. Bandwidths Comparison (h.x and h.y)
      bw_sub <- bw_storage[bw_storage$n == max(n_vec) & !is.na(bw_storage$hx), ]
      bw_list <- list(
        "S50.hx" = bw_sub$hx[bw_sub$method==methods_list[1]], 
        "S100.hx"= bw_sub$hx[bw_sub$method==methods_list[2]], 
        "S200.hx"= bw_sub$hx[bw_sub$method==methods_list[3]],
        "S50.hy" = bw_sub$hy[bw_sub$method==methods_list[1]], 
        "S100.hy"= bw_sub$hy[bw_sub$method==methods_list[2]], 
        "S200.hy"= bw_sub$hy[bw_sub$method==methods_list[3]]
      )
      boxplot(bw_list, main="Bandwidths (Max n)", col=rep(cols, 2), 
              las=2, cex.axis=0.7, outline=FALSE, notch=TRUE)
      abline(v=3.5, lty=2, col="gray")
      
      dev.flush()
    }
  }
}

cat("\n\nSimulation Complete.\n")

# Final summary table
final_stats <- do.call(rbind, lapply(n_vec, function(n_val) {
  mse_n <- results_mse[results_mse$n == n_val, ]
  res_n <- results_time[results_time$n == n_val, ]
  data.frame(n = n_val, 
             Avg_Time_50  = mean(res_n$seconds[res_n$method==methods_list[1]]),
             Avg_Time_200 = mean(res_n$seconds[res_n$method==methods_list[3]]),
             Avg_MSE_50   = mean(mse_n$mse[mse_n$method==methods_list[1]]),
             Avg_MSE_100  = mean(mse_n$mse[mse_n$method==methods_list[2]]),
             Avg_MSE_200  = mean(mse_n$mse[mse_n$method==methods_list[3]]))
}))
print(final_stats, row.names = FALSE)
