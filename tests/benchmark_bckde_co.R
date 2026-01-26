library(bkcde)
set.seed(42)
n <- 100
x <- runif(n, 0, 1)
y <- rnorm(n, mean = 2 * sin(4 * pi * x), sd = 1 + abs(x))

# Benchmark core allocations with a small grid
res <- bkcdeco(x, y,
                optim.degree.cores.min = 1,
                optim.degree.cores.max = 6,
                optim.nmulti.cores.min = 1,
                optim.nmulti.cores.max = 6,
                optim.ksum.cores.min = 1,
                optim.ksum.cores.max = 2,
                fitted.cores.max = 6,
                proper.cores.max = 6,
                progress = TRUE)

summary(res)

# Load microbenchmark
if (!require(microbenchmark)) install.packages("microbenchmark")
library(microbenchmark)

# Benchmark default and tuned bkcde models
timing_results <- microbenchmark(
  default = {
    f.yx.default <- bkcde(x = x, y = y,
                          bwmethod = "cv.ml",
                          proper = TRUE)
  },
  tuned = {
    f.yx.tuned <- bkcde(x = x, y = y,
                        bwmethod = "cv.ml",
                        optim.cores = "manual",
                        optim.degree.cores = res$best.optim$optim.degree.cores,
                        optim.nmulti.cores = res$best.optim$optim.nmulti.cores,
                        optim.ksum.cores = res$best.optim$optim.ksum.cores,
                        fitted.cores = res$best.fitted$fitted.cores,
                        proper = TRUE,
                        proper.cores = res$best.proper$proper.cores)
  },
  times = 10 # Adjust as needed
)

# Print timing comparison
print(timing_results)
# Extract mean and median times (in seconds)
mean_default <- mean(timing_results$time[timing_results$expr == "default"]) / 1e9
mean_tuned   <- mean(timing_results$time[timing_results$expr == "tuned"]) / 1e9
median_default <- median(timing_results$time[timing_results$expr == "default"]) / 1e9
median_tuned   <- median(timing_results$time[timing_results$expr == "tuned"]) / 1e9

# Calculate speedup
mean_speedup <- mean_default / mean_tuned
median_speedup <- median_default / median_tuned

# Print summary
cat(sprintf("Mean Tuned total (sec): %.4f\n", mean_tuned))
cat(sprintf("Mean Default total (sec): %.4f\n", mean_default))
cat(sprintf("Mean speedup vs default: %.2fx faster (saved %.4f sec)\n\n", mean_speedup, mean_default - mean_tuned))

cat(sprintf("Median Tuned total (sec): %.4f\n", median_tuned))
cat(sprintf("Median Default total (sec): %.4f\n", median_default))
cat(sprintf("Median speedup vs default: %.2fx faster (saved %.4f sec)\n", median_speedup, median_default - median_tuned))
boxplot(timing_results, main = "Timing: Default vs Tuned bkcde", outline = FALSE, notch = TRUE)
