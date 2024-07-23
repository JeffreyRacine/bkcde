\name{bkcde.default}
\alias{bkcde.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bkcde.default(h = NULL, x = NULL, y = NULL, x.eval = NULL, y.eval = NULL, y.lb = NULL, y.ub = NULL, x.lb = NULL, x.ub = NULL, poly.raw = TRUE, degree = 0, degree.min = 0, degree.max = 5, proper = TRUE, nmulti = 5, n.int = 100, ksum.cores = 1, degree.cores = NULL, nmulti.cores = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x.eval}{
%%     ~~Describe \code{x.eval} here~~
}
  \item{y.eval}{
%%     ~~Describe \code{y.eval} here~~
}
  \item{y.lb}{
%%     ~~Describe \code{y.lb} here~~
}
  \item{y.ub}{
%%     ~~Describe \code{y.ub} here~~
}
  \item{x.lb}{
%%     ~~Describe \code{x.lb} here~~
}
  \item{x.ub}{
%%     ~~Describe \code{x.ub} here~~
}
  \item{poly.raw}{
%%     ~~Describe \code{poly.raw} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{degree.min}{
%%     ~~Describe \code{degree.min} here~~
}
  \item{degree.max}{
%%     ~~Describe \code{degree.max} here~~
}
  \item{proper}{
%%     ~~Describe \code{proper} here~~
}
  \item{nmulti}{
%%     ~~Describe \code{nmulti} here~~
}
  \item{n.int}{
%%     ~~Describe \code{n.int} here~~
}
  \item{ksum.cores}{
%%     ~~Describe \code{ksum.cores} here~~
}
  \item{degree.cores}{
%%     ~~Describe \code{degree.cores} here~~
}
  \item{nmulti.cores}{
%%     ~~Describe \code{nmulti.cores} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (h = NULL, x = NULL, y = NULL, x.eval = NULL, y.eval = NULL, 
    y.lb = NULL, y.ub = NULL, x.lb = NULL, x.ub = NULL, poly.raw = TRUE, 
    degree = 0, degree.min = 0, degree.max = 5, proper = TRUE, 
    nmulti = 5, n.int = 100, ksum.cores = 1, degree.cores = NULL, 
    nmulti.cores = NULL, ...) 
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
    if (nmulti < 1) 
        stop("nmulti must be at least 1 in bkcde()")
    if (n.int < 1) 
        stop("n.int must be at least 1 in bkcde()")
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
            ...)
        h <- optim.out$par
        h.mat <- optim.out$par.mat
        degree <- optim.out$degree
        value <- optim.out$value
        value.vec <- optim.out$value.vec
        convergence <- optim.out$convergence
        convergence.vec <- optim.out$convergence.vec
        secs.optim <- optim.out$secs.optim
    }
    else {
        h.mat <- NULL
        value <- NULL
        value.vec <- NULL
        convergence <- NULL
        convergence.vec <- NULL
        secs.optim <- NULL
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
    f.yx[!is.finite(f.yx) | f.yx <= 0] <- .Machine$double.xmin
    if (proper) {
        if (is.finite(y.lb) && is.finite(y.ub)) 
            y.seq <- seq(y.lb, y.ub, length = n.int)
        if (is.finite(y.lb) && !is.finite(y.ub)) 
            y.seq <- seq(y.lb, extendrange(y, f = 10)[2], length = n.int)
        if (!is.finite(y.lb) && is.finite(y.ub)) 
            y.seq <- seq(extendrange(y, f = 10)[1], y.ub, length = n.int)
        if (!is.finite(y.lb) && !is.finite(y.ub)) 
            y.seq <- seq(extendrange(y, f = 10)[1], extendrange(y, 
                f = 10)[2], length = n.int)
        K <- kernel.bk(x.eval[1], x, h[2], x.lb, x.ub)
        if (degree == 0) {
            f.seq <- as.numeric(mcmapply(function(i) {
                mean(kernel.bk(y.seq[i], y, h[1], y.lb, y.ub) * 
                  K)/NZD(mean(K))
            }, 1:n.int, mc.cores = ksum.cores))
        }
        else {
            X.poly <- poly(x, raw = poly.raw, degree = degree)
            X <- cbind(1, X.poly)
            X.eval <- cbind(1, predict(X.poly, x.eval[1]))
            f.seq <- as.numeric(mcmapply(function(i) {
                coef(lm.wfit(x = X, y = kernel.bk(y.seq[i], y, 
                  h[1], y.lb, y.ub), w = NZD(K))) \%*\% t(X.eval)
            }, 1:n.int, mc.cores = ksum.cores))
        }
        f.seq[!is.finite(f.seq) | f.seq <= 0] <- .Machine$double.xmin
        f.yx[!is.finite(f.yx) | f.yx <= 0] <- .Machine$double.xmin
        int.f.seq <- integrate.trapezoidal(y.seq, f.seq)[length(y.seq)]
        f.yx <- f.yx/int.f.seq
    }
    return.list <- list(f = f.yx, h = h, h.mat = h.mat, degree = degree, 
        value = value, value.vec = value.vec, convergence = convergence, 
        convergence.vec = convergence.vec, secs.elapsed = as.numeric(difftime(Sys.time(), 
            secs.start.total, units = "secs")), secs.optim.vec = as.numeric(secs.optim), 
        secs.estimate = as.numeric(difftime(Sys.time(), secs.start.estimate, 
            units = "secs")), f.yx.integral = int.f.seq, y.lb = y.lb, 
        y.ub = y.ub, x.lb = x.lb, x.ub = x.ub, y = y, x = x, 
        y.eval = y.eval, x.eval = x.eval)
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
