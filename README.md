# bkcde: Boundary Corrected Polynomial Adaptive Kernel Conditional Density Estimation

This repository contains the R package **bkcde** for fast, flexible, and robust nonparametric estimation of conditional densities $f(y|x)$, with automatic boundary correction and adaptive polynomial order selection. Written and maintained by Jeffrey S. Racine (<racinej@mcmaster.ca>).

## Features

- Estimates conditional densities $f(y|x)$ using boundary-corrected, polynomial-adaptive kernel methods
- Automatically selects bandwidths and polynomial order via likelihood or least-squares cross-validation
- Handles bounded responses, with or without known bounds
- Efficient for large data via sub-sampled cross-validation and parallel processing
- S3 methods: `plot`, `predict`, and `summary` for `bkcde` objects
- Includes benchmarking for optimal core allocation (`bkcdeco`)

## Installation

Install dependencies and the development version from GitHub:

```r
install.packages(c("MCPAN", "pbmcapply", "robustbase"), repos="https://cran.wu.ac.at/")
# Optionally: install.packages(c("lpcde", "progress"))
library(devtools)
install_github('JeffreyRacine/bkcde')
```

**Note:**
- Windows users may need [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- macOS users may need [Xcode](https://apps.apple.com/us/app/xcode/id497799835) and command line tools (`xcode-select --install`)
- Parallelism uses forking, which is not available on Windows

## Main Functions

- `bkcde()`: Main estimator. Many arguments for bandwidth, degree, bounds, cross-validation, parallelism, and more.
- `plot.bkcde()`: 2D/3D plotting with bootstrapped confidence intervals.
- `predict.bkcde()`: Predicts conditional density at new data.
- `sub.cv()`: Efficient bandwidth and degree selection for large data.
- `bkcdeco()`: Benchmarks core allocation for parallelism.

## Example Usage

```r
library(bkcde)
set.seed(42)
n <- 250
x <- runif(n, -.25, .25)
s1 <- 1
s2 <- 1.25
y <- rbeta(n, s1 + x, s2 + x)
# Estimate conditional density
f.yx <- bkcde(x = x, y = y)
summary(f.yx)
# Plot with Bonferroni confidence intervals
plot(f.yx, ci = TRUE, ci.method = "Bonferroni")
# Predict at new points
predict(f.yx, data.frame(x = c(0, .1), y = c(.25, 0.5)))
```

### Large Data Example

```r
# Efficient sub-sampled cross-validation for large n
set.seed(42)
n <- 1e6
x <- runif(n, 1, 4)
y <- rbeta(n, 1.25 + x, 1.25 + x)
f.yx <- bkcde(x = x, y = y, n.grid = 25, cv="sub")
plot(f.yx, n.grid = 25, theta = 120, phi = 45, main = paste("Estimate (n =", n, ")"))
```

## Parallelism

Extensive use of parallel processing; users can control core allocation for different steps (see `?bkcde`, `?bkcdeco`).

## References

- Racine, J.S. (1993), "An Efficient Cross-Validation Algorithm For Window Width Selection for Nonparametric Kernel Regression," Communications in Statistics, 22(4), 1107-1114.
- Delaigle, A. and I. Gijbels (2004), "Bootstrap bandwidth selection in kernel density estimation from a contaminated sample," Annals of the Institute of Statistical Mathematics, 56, 19-47.

## Notes

- For large data, use `cv="sub"` for efficiency.
- See the package documentation and demos for advanced usage and simulation studies.


