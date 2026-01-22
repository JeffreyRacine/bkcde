################################################################################
## Bivariate Boundary Binning Performance Test
## Jeff Racine - Feasible Resampling vs. Binned vs. Unbinned
################################################################################

## User-defined parameters
M <- 100                    # Number of replications
n_vec <- c(500, 1000, 1500) 
num_bins <- 50              
n_sub_val <- 250            # n.sub parameter for feasible resampling
eval_grid_size <- 50        # Used to create evaluation data
plot_enabled <- TRUE        

library(bkcde)
library(parallel)

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
    
    # 1. Generate Data
    data <- dgp_beta(n_curr)
    x.lb <- min(data$x); x.ub <- max(data$x)
    y.lb <- min(data$y); y.ub <- max(data$y)
    
    # 2. Define Evaluation Grid
    x_seq <- seq(x.lb, x.ub, length.out = eval_grid_size)
    y_seq <- seq(y.lb, y.ub, length.out = eval_grid_size)
    eval_grid <- expand.grid(x = x_seq, y = y_seq)
    f_true <- true_density(eval_grid$x, eval_grid$y)
    
    # --- METHOD 1: Unbinned (Full CV) ---
    # Verified args: cv="full", cv.binned=FALSE
    t1 <- proc.time()
    fit_unb <- bkcde(x = data$x, y = data$y, 
                     x.eval = eval_grid$x, y.eval = eval_grid$y, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "full", cv.binned = FALSE)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Unbinned", "seconds"] <- (proc.time() - t1)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Unbinned", "mse"] <- mean((fit_unb$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Unbinned", c("hx","hy")] <- fit_unb$h[1:2]
    
    # --- METHOD 2: Binned (Full CV with Binning) ---
    # Verified args: cv="full", cv.binned=TRUE, n.binned=num_bins
    t2 <- proc.time()
    fit_bin <- bkcde(x = data$x, y = data$y, 
                     x.eval = eval_grid$x, y.eval = eval_grid$y, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "full", cv.binned = TRUE, n.binned = num_bins)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Binned", "seconds"] <- (proc.time() - t2)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Binned", "mse"] <- mean((fit_bin$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Binned", c("hx","hy")] <- fit_bin$h[1:2]
    
    # --- METHOD 3: Subsample (Feasible Resampling CV) ---
    # Verified args: cv="sub", n.sub=n_sub_val
    # Implements dimension-free weight-based resampling[cite: 14, 119].
    t3 <- proc.time()
    fit_sub <- bkcde(x = data$x, y = data$y, 
                     x.eval = eval_grid$x, y.eval = eval_grid$y, 
                     x.lb = x.lb, x.ub = x.ub, y.lb = y.lb, y.ub = y.ub,
                     cv = "sub", n.sub = n_sub_val)
    results_time[results_time$n == n_curr & results_time$M == m & results_time$method == "Subsample", "seconds"] <- (proc.time() - t3)[3]
    results_mse[results_mse$n == n_curr & results_mse$M == m & results_mse$method == "Subsample", "mse"] <- mean((fit_sub$f - f_true)^2)
    bw_storage[bw_storage$n == n_curr & bw_storage$M == m & bw_storage$method == "Subsample", c("hx","hy")] <- fit_sub$h[1:2]
    
    # Live Plotting
    if (plot_enabled && (m > 2)) {
      dev.hold()
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
      
      cols <- c("black", "blue", "red")
      df_t <- results_time[!is.na(results_time$seconds), ]
      medians <- sapply(methods_list, function(mt) sapply(n_vec, function(nv) median(df_t$seconds[df_t$method==mt & df_t$n==nv], na.rm=T)))
      
      matplot(n_vec, medians, type="l", col=cols, lty=1, lwd=2, xlab="n", ylab="Seconds", main="Scaling Performance")
      legend("topleft", methods_list, col=cols, lty=1, bty="n", cex=0.8)
      
      d_s <- results_mse[results_mse$n == min(n_vec) & !is.na(results_mse$mse), ]
      boxplot(mse ~ method, data=d_s, main=paste("MSE (n =", min(n_vec), ")"), 
              outline=FALSE, notch=TRUE, col="lightgray")
      
      d_l <- results_mse[results_mse$n == max(n_vec) & !is.na(results_mse$mse), ]
      boxplot(mse ~ method, data=d_l, main=paste("MSE (n =", max(n_vec), ")"), 
              outline=FALSE, notch=TRUE, col="lightblue")
      
      bw_sub <- bw_storage[bw_storage$n == max(n_vec) & !is.na(bw_storage$hx), ]
      bw_list <- list("U.hx"=bw_sub$hx[bw_sub$method=="Unbinned"], "B.hx"=bw_sub$hx[bw_sub$method=="Binned"], "S.hx"=bw_sub$hx[bw_sub$method=="Subsample"],
                      "U.hy"=bw_sub$hy[bw_sub$method=="Unbinned"], "B.hy"=bw_sub$hy[bw_sub$method=="Binned"], "S.hy"=bw_sub$hy[bw_sub$method=="Subsample"])
      boxplot(bw_list, main="Bandwidths (Max n)", col=rep(cols, 2), las=2, cex.axis=0.7, outline=FALSE, notch=TRUE)
      abline(v=3.5, lty=2, col="gray")
      
      dev.flush()
    }
  }
}

cat("\n\nSimulation Complete.\n")

# Summary Table
final_stats <- do.call(rbind, lapply(n_vec, function(n_val) {
  res_n <- results_time[results_time$n == n_val, ]
  mse_n <- results_mse[results_mse$n == n_val, ]
  
  data.frame(
    n = n_val, 
    Avg_Speedup_Bin = mean(res_n$seconds[res_n$method=="Unbinned"]) / mean(res_n$seconds[res_n$method=="Binned"]),
    Avg_MSE_Unbinned = mean(mse_n$mse[mse_n$method=="Unbinned"]),
    Avg_MSE_Binned   = mean(mse_n$mse[mse_n$method=="Binned"]),
    Avg_MSE_Subsample = mean(mse_n$mse[mse_n$method=="Subsample"])
  )
}))

print(final_stats, row.names = FALSE)