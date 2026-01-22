## Benchmark script for bkcde optimization
## Based on the DGP from demo/normal.R
library(bkcde)
library(microbenchmark)

## 1. DATA GENERATING PROCESS (DGP)
set.seed(42)
n <- 5000  # Increased n to highlight the optimization gap
x <- runif(n, 0, 1)
y <- rnorm(n, mean = x, sd = 0.2)

## Define bounds
x.lb <- 0; x.ub <- 1
y.lb <- min(y); y.ub <- max(y)

## 2. BENCHMARKING
## We compare full optimization (cross-validation) for both methods.
## Note: This assumes the 'linear-binning' path now uses the 
## pre-computed X.act matrix internally.

bm <- microbenchmark(
  unbinned = {
    bkcde(x=x, y=y, 
          x.lb=x.lb, x.ub=x.ub, 
          y.lb=y.lb, y.ub=y.ub,
          cv.binned=FALSE,
          degree.min=1, degree.max=1, # Fixed degree to isolate speed
          nmulti=1,
          cv.only=TRUE)
  },
  binned_optimized = {
    bkcde(x=x, y=y, 
          x.lb=x.lb, x.ub=x.ub, 
          y.lb=y.lb, y.ub=y.ub,
          cv.binned=TRUE,
          n.binned=100,
          degree.min=1, degree.max=1,
          nmulti=1,
          cv.only=TRUE)
  },
  times = 5
)

print(bm)

## 3. VERIFICATION OF ACCURACY
## Ensure the bandwidths (h) found by both methods are similar
fit.unbinned <- bkcde(x=x, y=y, x.lb=x.lb, x.ub=x.ub, y.lb=y.lb, y.ub=y.ub,
                      cv.binned=FALSE, degree=1, nmulti=1)

fit.binned   <- bkcde(x=x, y=y, x.lb=x.lb, x.ub=x.ub, y.lb=y.lb, y.ub=y.ub,
                      cv.binned=TRUE, degree=1, nmulti=1)

cat("\nBandwidths (Unbinned):", fit.unbinned$h)
cat("\nBandwidths (Binned):  ", fit.binned$h)
