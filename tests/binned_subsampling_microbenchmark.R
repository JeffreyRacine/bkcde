## Benchmark script for bkcde optimization: Comparing Unbinned, Binned, and Sub-sampled
library(bkcde)
library(microbenchmark)

## 1. DATA GENERATING PROCESS (DGP)
set.seed(42)
n <- 2500  # Increased n to 5000 to clearly see the O(n^2) explosion
x <- runif(n, 0, 1)
y <- rnorm(n, mean = x, sd = 0.2)

## Define bounds based on data range to avoid "out of bounds" errors
x.lb <- 0; x.ub <- 1
y.lb <- min(y); y.ub <- max(y)

## 2. BENCHMARKING
## Sub-sampling at 10% (n=500) vs Binned (M=100) vs Full Unbinned (n=5000)
bm <- microbenchmark(
  unbinned_full = {
    bkcde(x=x, y=y, 
          x.lb=x.lb, x.ub=x.ub, 
          y.lb=y.lb, y.ub=y.ub,
          cv.binned=FALSE,
          degree.min=1, degree.max=1,
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
  sub_sampled = {
    bkcde(x=x, y=y, 
          x.lb=x.lb, x.ub=x.ub, 
          y.lb=y.lb, y.ub=y.ub,
          cv="sub",             # Enable sub-sampling
          cv.sub.size=500,      # Use 500 observations for CV
          degree.min=1, degree.max=1,
          nmulti=1,
          cv.only=TRUE)
  },
  times = 3 # Reduced times because 'unbinned_full' will be slow at n=5000
)

print(bm)

## 3. SUMMARY OF EXPECTATIONS
## - Unbinned: Should be the slowest by far (O(n^2)).
## - Sub-sampled: Speed depends on 'cv.sub.size'. At 500, it should be much faster than full.
## - Binned: Should be the most stable. Its speed is tied to 'n.binned', not 'n'.
