\name{bkcde.default}
\alias{bkcde.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default bkcde Function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bkcde.default(h = NULL, 
              x = NULL, 
              y = NULL, 
              x.eval = NULL, 
              y.eval = NULL, 
              x.lb = NULL, 
              y.lb = NULL, 
              x.ub = NULL, 
              y.ub = NULL, 
              degree.cores = NULL, 
              degree.max = 5, 
              degree.min = 0, 
              degree = 0, 
              ksum.cores = 1, 
              max.pen.neg.loo = TRUE, 
              n.integrate = 100, 
              nmulti.cores = NULL, 
              nmulti = 5, 
              poly.raw = TRUE, 
              proper = TRUE, 
              ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
  a bandwidth vector of length 2, where \code{h[1]} is the bandwidth for the response variable and \code{h[2]} is the bandwidth for the predictor variable
}
  \item{x}{
  a vector of predictor variable values
}
  \item{y}{
  a vector of response variable values
}
  \item{x.eval}{
  a vector of predictor variable values at which to evaluate the density estimate
}
  \item{y.eval}{
  a vector of response variable values at which to evaluate the density estimate
}
  \item{x.lb}{
  the lower bound of the predictor variable
}
  \item{y.lb}{
  the lower bound of the response variable
}
  \item{x.ub}{
  the upper bound of the predictor variable
}
  \item{y.ub}{
  the upper bound of the response variable
}
  \item{degree.cores}{
  the number of cores to use in parallel processing for the degree optimization
}
  \item{degree.max}{
  the maximum degree of the polynomial to use in the density estimate
}
  \item{degree.min}{
  the minimum degree of the polynomial to use in the density estimate
}
  \item{degree}{
  the degree of the polynomial to use in the density estimate
}
  \item{ksum.cores}{
  the number of cores to use in parallel processing for the kernel sum
}
  \item{max.pen.neg.loo}{
  a logical indicating whether to use the maximum penalty for negative leave-one-out cross-validation estimates
}
  \item{n.integrate}{
  the number of points to use in the numerical integration of the density estimate
}
  \item{nmulti.cores}{
  the number of cores to use in parallel processing for the multi-dimensional optimization
}
  \item{nmulti}{
  the number of multi-dimensional optimization runs to use in the density estimate
}
  \item{poly.raw}{
  a logical indicating whether to use raw or orthogonal polynomials in the density estimate
}
  \item{proper}{
  a logical indicating whether to normalize the density estimate to integrate to 1
}
  \item{\dots}{
  additional arguments to be passed to the optimization function
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
Jeffrey S. Racine <racinej@mcmaster.ca>
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
function (h = NULL, x = NULL, y = NULL, x.eval = NULL, y.eval = NULL, 
    x.lb = NULL, y.lb = NULL, x.ub = NULL, y.ub = NULL, degree.cores = NULL, 
    degree.max = 5, degree.min = 0, degree = 0, ksum.cores = 1, 
    max.pen.neg.loo = TRUE, n.integrate = 100, nmulti.cores = NULL, 
    nmulti = 5, poly.raw = TRUE, proper = TRUE, ...) 
{
    if (is.null(x)) 
        stop("must provide x in bkcde()")
    if (is.null(y)) 
        stop("must provide y in bkcde()")
    if (is.null(x.eval)) 
        stop("must provide x.eval in bkcde()")
    if (is.null(y.eval)) 
        stop("must provide y.eval in bkcde()")
    if (is.null(y.lb)) 
        y.lb <- min(y)
    if (is.null(y.ub)) 
        y.ub <- max(y)
    if (is.null(x.lb)) 
        x.lb <- min(x)
    if (is.null(x.ub)) 
        x.ub <- max(x)
    if (any(y < y.lb) | any(y > y.ub)) 
        stop("y must lie in [y.lb,y.ub] in bkcde()")
    if (any(y.eval < y.lb) | any(y.eval > y.ub)) 
        stop("y.eval must lie in [y.lb,y.ub] in bkcde()")
    if (any(x < x.lb) | any(x > x.ub)) 
        stop("x must lie in [x.lb,x.ub] in bkcde()")
    if (any(x.eval < x.lb) | any(x.eval > x.ub)) 
        stop("x.eval must lie in [x.lb,x.ub] in bkcde()")
    if (y.lb >= y.ub) 
        stop("y.lb must be less than y.ub in bkcde()")
    if (x.lb >= x.ub) 
        stop("x.lb must be less than x.ub in bkcde()")
    if (!is.logical(poly.raw)) 
        stop("poly.raw must be logical in bkcde()")
    if (!is.logical(proper)) 
        stop("proper must be logical in bkcde()")
    if (!is.logical(max.pen.neg.loo)) 
        stop("max.pen.neg.loo must be logical in bkcde()")
    if (nmulti < 1) 
        stop("nmulti must be at least 1 in bkcde()")
    if (n.integrate < 1) 
        stop("n.integrate must be at least 1 in bkcde()")
    if (degree < 0 | degree >= length(y)) 
        stop("degree must lie in [0,1,...,", length(y) - 1, "] (i.e., [0,1,dots, n-1]) in bkcde()")
    if (degree.min < 0 | degree.min >= length(y)) 
        stop("degree.min must lie in [0,1,...,", length(y) - 
            1, "] (i.e., [0,1,dots, n-1]) in bkcde()")
    if (degree.max < 0 | degree.max >= length(y)) 
        stop("degree.max must lie in [0,1,...,", length(y) - 
            1, "] (i.e., [0,1,dots, n-1]) in bkcde()")
    if (degree.min > degree.max) 
        stop("degree.min must be <= degree.max in bkcde()")
    if (ksum.cores < 1) 
        stop("ksum.cores must be at least 1 in bkcde()")
    if (is.null(degree.cores)) 
        degree.cores <- degree.max - degree.min + 1
    if (is.null(nmulti.cores)) 
        nmulti.cores <- nmulti
    secs.start.total <- Sys.time()
    if (is.null(h)) {
        optim.out <- bkcde.optim(x = x, y = y, y.lb = y.lb, y.ub = y.ub, 
            x.lb = x.lb, x.ub = x.ub, poly.raw = poly.raw, degree.min = degree.min, 
            degree.max = degree.max, nmulti = nmulti, ksum.cores = ksum.cores, 
            degree.cores = degree.cores, nmulti.cores = nmulti.cores, 
            max.pen.neg.loo = max.pen.neg.loo, ...)
        h <- optim.out$par
        h.mat <- optim.out$par.mat
        degree <- optim.out$degree
        degree.mat <- optim.out$degree.mat
        value <- optim.out$value
        value.vec <- optim.out$value.vec
        value.mat <- optim.out$value.mat
        convergence <- optim.out$convergence
        convergence.vec <- optim.out$convergence.vec
        convergence.mat <- optim.out$convergence.mat
        secs.optim <- optim.out$secs.optim
        secs.optim.mat <- optim.out$secs.optim.mat
    }
    else {
        h.mat <- NULL
        degree.mat <- NULL
        value <- NULL
        value.vec <- NULL
        value.mat <- NULL
        convergence <- NULL
        convergence.vec <- NULL
        convergence.mat <- NULL
        secs.optim <- NULL
        secs.optim.mat <- NULL
    }
    secs.start.estimate <- Sys.time()
    if (degree == 0) {
        f.yx <- as.numeric(mcmapply(function(i) {
            kernel.bk.x <- kernel.bk(x.eval[i], x, h[2], x.lb, 
                x.ub)
            mean(kernel.bk(y.eval[i], y, h[1], y.lb, y.ub) * 
                kernel.bk.x)/NZD(mean(kernel.bk.x))
        }, 1:length(y.eval), mc.cores = ksum.cores))
    }
    else {
        X.poly <- poly(x, raw = poly.raw, degree = degree)
        X <- cbind(1, X.poly)
        f.yx <- as.numeric(mcmapply(function(i) {
            coef(lm.wfit(x = X, y = kernel.bk(y.eval[i], y, h[1], 
                y.lb, y.ub), w = NZD(kernel.bk(x.eval[i], x, 
                h[2], x.lb, x.ub)))) \%*\% t(cbind(1, predict(X.poly, 
                x.eval[i])))
        }, 1:length(y.eval), mc.cores = ksum.cores))
    }
    if (proper) {
        if (is.finite(y.lb) && is.finite(y.ub)) 
            y.seq <- seq(y.lb, y.ub, length = n.integrate)
        if (is.finite(y.lb) && !is.finite(y.ub)) 
            y.seq <- seq(y.lb, extendrange(y, f = 10)[2], length = n.integrate)
        if (!is.finite(y.lb) && is.finite(y.ub)) 
            y.seq <- seq(extendrange(y, f = 10)[1], y.ub, length = n.integrate)
        if (!is.finite(y.lb) && !is.finite(y.ub)) 
            y.seq <- seq(extendrange(y, f = 10)[1], extendrange(y, 
                f = 10)[2], length = n.integrate)
        K <- kernel.bk(x.eval[1], x, h[2], x.lb, x.ub)
        if (degree == 0) {
            f.seq <- as.numeric(mcmapply(function(i) {
                mean(kernel.bk(y.seq[i], y, h[1], y.lb, y.ub) * 
                  K)/NZD(mean(K))
            }, 1:n.integrate, mc.cores = ksum.cores))
        }
        else {
            X.poly <- poly(x, raw = poly.raw, degree = degree)
            X <- cbind(1, X.poly)
            X.eval <- cbind(1, predict(X.poly, x.eval[1]))
            f.seq <- as.numeric(mcmapply(function(i) {
                coef(lm.wfit(x = X, y = kernel.bk(y.seq[i], y, 
                  h[1], y.lb, y.ub), w = NZD(K))) \%*\% t(X.eval)
            }, 1:n.integrate, mc.cores = ksum.cores))
        }
        if (any(!is.finite(f.yx))) 
            warning("non-finite density estimate reset to 0 via option proper=TRUE in bkcde()")
        f.yx[!is.finite(f.yx)] <- 0
        if (any(f.yx < 0)) 
            warning("negative density estimate reset to 0 via option proper=TRUE in bkcde()")
        f.yx[f.yx < 0] <- 0
        f.seq[!is.finite(f.seq) | f.seq < 0] <- 0
        int.f.seq <- integrate.trapezoidal(y.seq, f.seq)[length(y.seq)]
        f.yx <- f.yx/int.f.seq
    }
    else {
        int.f.seq <- NULL
        if (any(!is.finite(f.yx))) 
            warning("non-finite density estimate encountered with option proper=FALSE, reset to 0 in bkcde()")
        f.yx[!is.finite(f.yx)] <- 0
        if (any(f.yx < 0)) 
            warning("negative density estimate encountered, consider option proper=TRUE in bkcde()")
    }
    return.list <- list(convergence.mat = convergence.mat, convergence.vec = convergence.vec, 
        convergence = convergence, degree.cores = degree.cores, 
        degree.mat = degree.mat, degree.max = degree.max, degree.min = degree.min, 
        degree = degree, f.yx.integral = int.f.seq, f = f.yx, 
        h.mat = h.mat, h = h, ksum.cores = ksum.cores, nmulti.cores = nmulti.cores, 
        max.pen.neg.loo = max.pen.neg.loo, secs.elapsed = as.numeric(difftime(Sys.time(), 
            secs.start.total, units = "secs")), secs.estimate = as.numeric(difftime(Sys.time(), 
            secs.start.estimate, units = "secs")), secs.optim.mat = secs.optim.mat, 
        value.mat = value.mat, value.vec = value.vec, value = value, 
        x.eval = x.eval, x.lb = x.lb, x.ub = x.ub, x = x, y.eval = y.eval, 
        y.lb = y.lb, y.ub = y.ub, y = y)
    class(return.list) <- "bkcde"
    return(return.list)
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
