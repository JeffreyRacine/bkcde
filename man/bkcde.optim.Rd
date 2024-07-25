\name{bkcde.optim}
\alias{bkcde.optim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimize the Bandwidth and Polynomial Degree of a Bivariate Kernel Conditional Density Estimator
}
\description{
This function optimizes the bandwidth and polynomial degree of a bivariate kernel conditional density estimator using leave-one-out likelihood cross-validation.
}
\usage{
bkcde.optim(x = x, 
            y = y, 
            y.lb = y.lb, 
            y.ub = y.ub, 
            x.lb = x.lb, 
            x.ub = x.ub, 
            poly.raw = poly.raw, 
            degree.min = degree.min, 
            degree.max = degree.max, 
            nmulti = nmulti, 
            ksum.cores = ksum.cores, 
            degree.cores = degree.cores, 
            nmulti.cores = nmulti.cores, 
            max.pen.neg.loo = max.pen.neg.loo,
            ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a vector of length $n$ containing the $x$-coordinates of the data points
}
  \item{y}{
  a vector of length $n$ containing the $y$-coordinates of the data points
}
  \item{y.lb}{
  the lower bound of the $y$-variable
}
  \item{y.ub}{
  the upper bound of the $y$-variable
}
  \item{x.lb}{
  the lower bound of the $x$-variable
}
  \item{x.ub}{
  the upper bound of the $x$-variable
}
  \item{poly.raw}{
  a logical indicating whether the raw polynomial basis should be used
}
  \item{degree.min}{
  the minimum polynomial degree to consider
}
  \item{degree.max}{
  the maximum polynomial degree to consider
}
  \item{nmulti}{
  the number of random initializations to use
}
  \item{ksum.cores}{
  the number of cores to use for the kernel sum
}
  \item{degree.cores}{
  the number of cores to use for the polynomial degree optimization
}
  \item{nmulti.cores}{
  the number of cores to use for the random initialization optimization
}
  \item{max.pen.neg.loo}{
  a logical value that indicates whether the maximum penalty for negative leave-one-out likelihoods should be used \code{.Machine\$double.xmin}}
}
  \item{\dots}{
  additional arguments to be passed to \code{optim}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x = x, y = y, y.lb = y.lb, y.ub = y.ub, x.lb = x.lb, 
    x.ub = x.ub, poly.raw = poly.raw, degree.min = degree.min, 
    degree.max = degree.max, nmulti = nmulti, ksum.cores = ksum.cores, 
    degree.cores = degree.cores, nmulti.cores = nmulti.cores, 
    max.pen.neg.loo = max.pen.neg.loo, ...) 
{
    if (degree.min < 0 | degree.max >= length(y)) 
        stop("degree.min must lie in [0,1,...,", length(y) - 
            1, "] (i.e., [0,1,dots, n-1]) in bkcde.optim()")
    if (degree.max < 0 | degree.max >= length(y)) 
        stop("degree.max must lie in [0,1,...,", length(y) - 
            1, "] (i.e., [0,1,dots, n-1]) in bkcde.optim()")
    if (degree.max < degree.min) 
        stop("degree.max must be >= degree.min in bkcde.optim()")
    if (missing(x)) 
        stop("must provide x in bkcde.optim()")
    if (missing(y)) 
        stop("must provide y in bkcde.optim()")
    if (missing(y.lb)) 
        stop("must provide y.lb in bkcde.optim()")
    if (missing(y.ub)) 
        stop("must provide y.ub in bkcde.optim()")
    if (missing(x.lb)) 
        stop("must provide x.lb in bkcde.optim()")
    if (missing(x.ub)) 
        stop("must provide x.ub in bkcde.optim()")
    if (!is.logical(poly.raw)) 
        stop("poly.raw must be logical in bkcde.optim()")
    n <- length(y)
    lower <- 0.1 * c(EssDee(y), EssDee(x)) * n^{
        -1/6
    }
    upper <- 1000 * c(EssDee(y), EssDee(x))
    degree.return <- mclapply(degree.min:degree.max, function(p) {
        nmulti.return <- mclapply(1:nmulti, function(i) {
            if (i == 1) {
                init <- c(EssDee(y), EssDee(x)) * n^{
                  -1/6
                }
            }
            else {
                init <- runif(2, 0.5, 5) * c(EssDee(y), EssDee(x)) * 
                  n^{
                    -1/6
                  }
            }
            st <- system.time(optim.return <- optim(par = init, 
                fn = bkcde.loo, x = x, y = y, y.lb = y.lb, y.ub = y.ub, 
                x.lb = x.lb, x.ub = x.ub, poly.raw = poly.raw, 
                degree = p, ksum.cores = ksum.cores, max.pen.neg.loo = max.pen.neg.loo, 
                lower = lower, upper = upper, method = "L-BFGS-B", 
                control = list(fnscale = -1)))
            optim.return$secs.optim <- st["elapsed"]
            optim.return$degree <- p
            optim.return
        }, mc.cores = nmulti.cores)
        optim.out <- nmulti.return[[which.max(sapply(nmulti.return, 
            function(x) x$value))]]
        optim.out$value.vec <- sapply(nmulti.return, function(x) x$value)
        optim.out$degree.vec <- sapply(nmulti.return, function(x) x$degree)
        optim.out$convergence.vec <- sapply(nmulti.return, function(x) x$convergence)
        optim.out$secs.optim.vec <- sapply(nmulti.return, function(x) x$secs.optim)
        optim.out
    }, mc.cores = degree.cores)
    output.return <- degree.return[[which.max(sapply(degree.return, 
        function(x) x$value))]]
    output.return$par.mat <- t(sapply(degree.return, function(x) x$par))
    output.return$value.vec <- sapply(degree.return, function(x) x$value)
    output.return$value.mat <- t(sapply(degree.return, function(x) x$value.vec))
    output.return$convergence.vec <- sapply(degree.return, function(x) x$convergence)
    output.return$convergence.mat <- t(sapply(degree.return, 
        function(x) x$convergence.vec))
    output.return$degree.mat <- t(sapply(degree.return, function(x) x$degree.vec))
    output.return$secs.optim <- sapply(degree.return, function(x) x$secs.optim)
    output.return$secs.optim.mat <- t(sapply(degree.return, function(x) x$secs.optim.vec))
    return(output.return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
