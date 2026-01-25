\name{bkcde.co}
\alias{bkcde.co}
\alias{summary.bkcde.co}
\alias{print.bkcde.co}

\title{
Core Allocation Benchmarking for bkcde
}
\description{
\code{bkcde.co} is a helper function that benchmarks core allocation strategies for the \code{\link{bkcde}} function by sweeping over optimization, fitting, and proper cores parameters. It identifies the fastest core configuration and compares its performance against the default (\code{optim.cores="auto"}) configuration. The function returns an object containing detailed timing results, optimal core settings, and parameter comparisons.
}
\usage{
bkcde.co(x, y,
         optim.degree.cores.min = 1,
         optim.degree.cores.max = parallel::detectCores(),
         optim.degree.cores.by  = 1,
         optim.nmulti.cores.min = 1,
         optim.nmulti.cores.max = parallel::detectCores(),
         optim.nmulti.cores.by  = 1,
         optim.ksum.cores.min   = 1,
         optim.ksum.cores.max   = 2,
         optim.ksum.cores.by    = 1,
         fitted.cores.min = 1,
         fitted.cores.max = parallel::detectCores(),
         fitted.cores.by  = 1,
         proper.cores.min = 1,
         proper.cores.max = parallel::detectCores(),
         proper.cores.by  = 1,
         progress = TRUE,
         display.warnings = FALSE,
         ...)

\method{summary}{bkcde.co}(object, ...)

\method{print}{bkcde.co}(x, ...)
}

\arguments{
  \item{x}{
  a numeric vector of the predictor variable.
}
  \item{y}{
  a numeric vector of the response variable.
}
  \item{optim.degree.cores.min}{
  minimum number of cores to test for degree optimization. Default is 1.
}
  \item{optim.degree.cores.max}{
  maximum number of cores to test for degree optimization. Default is \code{\link[parallel]{detectCores}}.
}
  \item{optim.degree.cores.by}{
  increment for degree cores sweep. Default is 1.
}
  \item{optim.nmulti.cores.min}{
  minimum number of cores to test for multi-start optimization. Default is 1.
}
  \item{optim.nmulti.cores.max}{
  maximum number of cores to test for multi-start optimization. Default is \code{\link[parallel]{detectCores}}.
}
  \item{optim.nmulti.cores.by}{
  increment for nmulti cores sweep. Default is 1.
}
  \item{optim.ksum.cores.min}{
  minimum number of cores to test for kernel sum optimization. Default is 1.
}
  \item{optim.ksum.cores.max}{
  maximum number of cores to test for kernel sum optimization. Default is 2.
}
  \item{optim.ksum.cores.by}{
  increment for ksum cores sweep. Default is 1.
}
  \item{fitted.cores.min}{
  minimum number of cores to test for fitting. Default is 1.
}
  \item{fitted.cores.max}{
  maximum number of cores to test for fitting. Default is \code{\link[parallel]{detectCores}}.
}
  \item{fitted.cores.by}{
  increment for fitted cores sweep. Default is 1.
}
  \item{proper.cores.min}{
  minimum number of cores to test for proper normalization. Default is 1.
}
  \item{proper.cores.max}{
  maximum number of cores to test for proper normalization. Default is \code{\link[parallel]{detectCores}}.
}
  \item{proper.cores.by}{
  increment for proper cores sweep. Default is 1.
}
  \item{progress}{
  logical; if TRUE (default), displays progress bars during benchmarking. Requires the \pkg{progress} package.
}
  \item{display.warnings}{
  logical; if FALSE (default), suppresses warnings from internal \code{bkcde} calls during benchmarking.
}
  \item{\dots}{
  additional arguments passed to \code{\link{bkcde}} (e.g., \code{bwmethod}, \code{degree.max}, etc.).
}
  \item{object}{
  an object of class \code{"bkcde.co"}.
}
}
\details{
\code{bkcde.co} performs a systematic grid search over core allocation parameters for the \code{\link{bkcde}} function. The benchmarking process consists of four stages:

\enumerate{
  \item \strong{Optimization sweep}: Tests all combinations of \code{optim.degree.cores}, \code{optim.nmulti.cores}, and \code{optim.ksum.cores} using \code{cv.only=TRUE} to identify the fastest optimization configuration.
  
  \item \strong{Fitting sweep}: Using the optimal bandwidth and degree from stage 1, tests different \code{fitted.cores} values to identify the fastest fitting configuration.
  
  \item \strong{Proper sweep}: Tests different \code{proper.cores} values with \code{proper=TRUE} using the optimal fitting cores from stage 2.
  
  \item \strong{Comparison}: Runs a complete \code{bkcde} call with the tuned core settings and compares total elapsed time against a baseline run using \code{optim.cores="auto"}.
}

The function automatically handles convergence failures by selecting the fastest non-error run when no configurations achieve \code{convergence==0}.

Progress bars (when \code{progress=TRUE}) provide visibility into each stage of the benchmarking process. Setting \code{display.warnings=FALSE} suppresses warnings about negative density values during intermediate sweeps.

The comparison between tuned and default configurations uses identical \code{...} arguments to ensure fair comparison. Both runs are complete end-to-end \code{bkcde} calls, capturing all overhead.
}
\value{
An object of class \code{"bkcde.co"}, which is a list containing:

\item{call}{the matched call}
\item{optim.grid}{data frame of all optimization core combinations tested}
\item{optim.results}{data frame of timing results for each optimization configuration}
\item{best.optim}{single-row data frame with the best optimization configuration}
\item{fit.opt}{the \code{bkcde} object from the optimization stage with best cores}
\item{fitted.grid}{data frame of all fitted cores values tested}
\item{fitted.results}{data frame of timing results for each fitted configuration}
\item{best.fitted}{single-row data frame with the best fitted configuration}
\item{proper.grid}{data frame of all proper cores values tested}
\item{proper.results}{data frame of timing results for each proper configuration}
\item{best.proper}{single-row data frame with the best proper configuration}
\item{tuned_total}{total elapsed time (seconds) for the tuned configuration}
\item{tuned_breakdown}{list with \code{optim} and \code{proper} timing breakdowns}
\item{tuned.fit}{the complete \code{bkcde} object from the tuned run}
\item{default.fit}{the complete \code{bkcde} object from the default run}
\item{default_total}{total elapsed time (seconds) for the default configuration}
\item{speedup_factor}{ratio of default_total / tuned_total}
\item{time_saved}{difference default_total - tuned_total (seconds)}
\item{params_match}{logical indicating if tuned and default converged to same parameters}
}
\references{
Racine, J.S. (1993), "An Efficient Cross-Validation Algorithm for Window Width Selection for Nonparametric Kernel Regression," Communications in Statistics - Simulation and Computation, 22, 1107-1120.
}
\author{
Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\seealso{
\code{\link{bkcde}}, \code{\link{sub.cv}}
}
\examples{
\dontrun{
# Generate synthetic data
set.seed(42)
n <- 250
x <- runif(n, 0, 1)
y <- rnorm(n, mean = 2 * sin(4 * pi * x), sd = 1 + abs(x))

# Benchmark core allocations with a small grid
res <- bkcde.co(x, y,
                optim.degree.cores.min = 1,
                optim.degree.cores.max = 2,
                optim.nmulti.cores.min = 1,
                optim.nmulti.cores.max = 2,
                optim.ksum.cores.min = 1,
                optim.ksum.cores.max = 2,
                fitted.cores.max = 4,
                proper.cores.max = 4,
                progress = TRUE)

# View summary
summary(res)

# Extract best core settings
cat("Best optim.degree.cores:", res$best.optim$optim.degree.cores, "\n")
cat("Best optim.nmulti.cores:", res$best.optim$optim.nmulti.cores, "\n")
cat("Best optim.ksum.cores:", res$best.optim$optim.ksum.cores, "\n")
cat("Best fitted.cores:", res$best.fitted$fitted.cores, "\n")
cat("Best proper.cores:", res$best.proper$proper.cores, "\n")

# Use tuned settings for production
f.yx <- bkcde(x = x, y = y,
              bwmethod = "cv.ml",
              optim.cores = "manual",
              optim.degree.cores = res$best.optim$optim.degree.cores,
              optim.nmulti.cores = res$best.optim$optim.nmulti.cores,
              optim.ksum.cores = res$best.optim$optim.ksum.cores,
              fitted.cores = res$best.fitted$fitted.cores,
              proper = TRUE,
              proper.cores = res$best.proper$proper.cores)
}
}
\keyword{nonparametric}
\keyword{smooth}
