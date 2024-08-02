\name{bkcde}
\alias{bkcde}
\alias{bkcde.call}
\alias{bkcde.default}

\title{
Boundary Corrected Polynomial Adaptive Conditional Kernel Density Estimation
}
\description{
bkcde() is a function that estimates the conditional density of a response variable given a predictor variable using a boundary corrected polynomial adaptive kernel density estimator. The function is designed to be fast and memory efficient by using parallel processing and vectorized operations. The function is also designed to be flexible by allowing the user to specify the degree of the polynomial and the bandwidths of the kernel functions manually, if so desired.
}
\usage{
bkcde(...)

\method{bkcde}{default}(h = NULL, 
      x = NULL, 
      y = NULL, 
      x.eval = NULL, 
      y.eval = NULL, 
      x.lb = NULL, 
      y.lb = NULL, 
      x.ub = NULL, 
      y.ub = NULL, 
      degree.max = 3, 
      degree.min = 0, 
      degree = 0, 
      ksum.cores = 1, 
      n.integrate = 1000, 
      nmulti = 3, 
      optim.degree.cores = NULL, 
      optim.nmulti.cores = NULL, 
      penalty.cutoff = .Machine$double.xmin, 
      penalty.method = c("smooth", "constant", "trim"), 
      poly.raw = FALSE, 
      proper.cores = 12, 
      proper = TRUE, 
      verbose = FALSE, 
      ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
  a bandwidth vector of length 2, where h[1] is the bandwidth of the kernel function for the response variable and h[2] is the bandwidth of the kernel function for the predictor variable. If h is NULL, the function will estimate the optimal bandwidths using the bkcde.optim() function.
}
  \item{x}{
  a numeric vector of the predictor variable.
}
  \item{y}{
  a numeric vector of the response variable.
}
  \item{x.eval}{
  a numeric vector of the predictor variable at which to evaluate the conditional density. If x.eval is NULL, the function will evaluate the conditional density at the sample x values
}
  \item{y.eval}{
  a numeric vector of the response variable at which to evaluate the conditional density. If y.eval is NULL, the function will evaluate the conditional density at the sample y values
}
  \item{x.lb}{
  the lower bound of the predictor variable. If x.lb is NULL, the function will use the minimum of x as the lower bound
}
  \item{y.lb}{
  the lower bound of the response variable. If y.lb is NULL, the function will use the minimum of y as the lower bound
}
  \item{x.ub}{
  the upper bound of the predictor variable. If x.ub is NULL, the function will use the maximum of x as the upper bound
}
  \item{y.ub}{
  the upper bound of the response variable. If y.ub is NULL, the function will use the maximum of y as the upper bound
}
  \item{degree.max}{
  the maximum degree of the polynomial searched over. The default is 3
}
  \item{degree.min}{
  the minimum degree of the polynomial searched. The default is 0
}
  \item{degree}{
  the degree of the polynomial if specified manually. The default is 0
}
  \item{ksum.cores}{
  the number of cores to use for the kernel sum. The default is 1
}
  \item{n.integrate}{
  the number of points to use for numerical integration. The default is 1000
}
  \item{nmulti}{
  the number of multi-starts to use for optimization. The default is 3
}
  \item{optim.degree.cores}{
  the number of cores to use for the optimization of the degree. The default is degree.max - degree.min + 1
}
  \item{optim.nmulti.cores}{
  the number of cores to use for the optimization of the multi-starts. The default is nmulti
}
  \item{penalty.cutoff}{
  the cutoff value for the penalty used in the log likelihood function for negative density values. The default is .Machine\$double.xmin
}
  \item{penalty.method}{
  the method used to penalize negative density values. The default is "smooth"
}
  \item{poly.raw}{
  a logical value indicating whether to use raw or orthogonal polynomials. The default is FALSE (i.e., orthogonal polynomials are the default)
}
  \item{proper.cores}{
  the number of cores to use for the proper normalization. The default is 12
}
  \item{proper}{
  a logical value indicating whether to use proper normalization. The default is TRUE
}
  \item{verbose}{
  a logical value indicating whether to print warnings. The default is FALSE
}
  \item{\dots}{
  additional arguments to be passed to the bkcde.optim() function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
return.list <- list(convergence.mat = convergence.mat, convergence.vec = convergence.vec, 
        convergence = convergence, degree.mat = degree.mat, degree.max = degree.max, 
        degree.min = degree.min, degree = degree, f.yx.integral.post = int.f.seq.post, 
        f.yx.integral.pre.neg = int.f.seq.pre.neg, f.yx.integral = int.f.seq, 
        f = f.yx, h.mat = h.mat, h = h, ksum.cores = ksum.cores, 
        optim.degree.cores = optim.degree.cores, optim.nmulti.cores = optim.nmulti.cores, 
        proper.cores = proper.cores, proper = proper, secs.elapsed = as.numeric(difftime(Sys.time(), 
            secs.start.total, units = "secs")), secs.estimate = as.numeric(difftime(Sys.time(), 
            secs.start.estimate, units = "secs")), secs.optim.mat = secs.optim.mat, 
        value.mat = value.mat, value.vec = value.vec, value = value, 
        x.eval = x.eval, x.lb = x.lb, x.ub = x.ub, x = x, y.eval = y.eval, 
        y.lb = y.lb, y.ub = y.ub, y = y)
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{convergence.mat}{
a matrix of convergence values for each multi-start and degree
}
\item{convergence.vec}{
a vector of convergence values for each multi-start
}
\item{convergence}{
a scalar convergence values reported by optim()
}
\item{degree.mat}{
a matrix of degrees for each multi-start
}
\item{degree.max}{
the maximum degree of the polynomial searched over
}
\item{degree.min}{
the minimum degree of the polynomial searched
}
\item{degree}{
the degree of the polynomial if specified manually
}
\item{f.yx.integral.post}{
the integral of the conditional density after proper normalization
}
\item{f.yx.integral.pre.neg}{
the integral of the conditional density before proper normalization
}
\item{f.yx.integral}{
the integral of the conditional density
}
\item{f}{
the conditional density
}
\item{h.mat}{
a matrix of bandwidths for each multi-start
}
\item{h}{
a bandwidth vector of length 2, where h[1] is the bandwidth of the kernel function for the response variable and h[2] is the bandwidth of the kernel function for the predictor variable
}
\item{ksum.cores}{
the number of cores to use for the kernel sum
}
\item{optim.degree.cores}{
the number of cores to use for the optimization of the degree
}
\item{optim.nmulti.cores}{
the number of cores to use for the optimization of the multi-starts
}
\item{proper.cores}{
the number of cores to use for the proper normalization
}
\item{proper}{
a logical value indicating whether to use proper normalization
}
\item{secs.elapsed}{
the total time elapsed in seconds
}
\item{secs.estimate}{
the time elapsed for the estimation in seconds
}
\item{secs.optim.mat}{
a matrix of times elapsed for each multi-start
}
\item{value.mat}{
a matrix of values for each multi-start
}
\item{value.vec}{
a vector of values for each multi-start
}
\item{value}{
a scalar value reported by optim()
}
\item{x.eval}{
a numeric vector of the predictor variable at which to evaluate the conditional density
}
\item{x.lb}{
the lower bound of the predictor variable
}
\item{x.ub}{
the upper bound of the predictor variable
}
\item{x}{
a numeric vector of the predictor variable
}
\item{y.eval}{
a numeric vector of the response variable at which to evaluate the conditional density
}
\item{y.lb}{
the lower bound of the response variable
}
\item{y.ub}{
the upper bound of the response variable
}
\item{y}{
a numeric vector of the response variable
}
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
    x.lb = NULL, y.lb = NULL, x.ub = NULL, y.ub = NULL, degree.max = 3, 
    degree.min = 0, degree = 0, ksum.cores = 1, n.integrate = 1000, 
    nmulti = 3, optim.degree.cores = NULL, optim.nmulti.cores = NULL, 
    penalty.cutoff = .Machine$double.xmin, penalty.method = c("smooth", 
        "constant", "trim"), poly.raw = FALSE, proper.cores = 12, 
    proper = TRUE, verbose = FALSE, ...) 
{
    if (is.null(x)) 
        stop("must provide x in bkcde()")
    if (is.null(y)) 
        stop("must provide y in bkcde()")
    if (!is.null(x.eval) & is.null(y.eval) & length(x.eval) != 
        length(y)) 
        stop("length of x.eval must be equal to length of y in bkcde() when y.eval is NULL")
    if (!is.null(x.eval) & !is.null(y.eval) & length(x.eval) != 
        length(y.eval)) 
        stop("length of x.eval must be equal to length of y.eval in bkcde() when x.eval and y.eval are not NULL")
    if (!is.null(y.eval) & is.null(x.eval)) 
        stop("must provide x.eval in bkcde() when y.eval is not NULL")
    if (is.null(x.eval)) 
        x.eval <- x
    if (is.null(y.eval)) 
        y.eval <- y
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
    if (!is.logical(verbose)) 
        stop("verbose must be logical in bkcde()")
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
    if (proper.cores < 1) 
        stop("proper.cores must be at least 1 in bkcde()")
    if (is.null(optim.degree.cores)) 
        optim.degree.cores <- degree.max - degree.min + 1
    if (is.null(optim.nmulti.cores)) 
        optim.nmulti.cores <- nmulti
    penalty.method <- match.arg(penalty.method)
    if (penalty.cutoff <= 0) 
        stop("penalty.cutoff must be positive in bkcde()")
    secs.start.total <- Sys.time()
    if (is.null(h)) {
        optim.out <- bkcde.optim(x = x, y = y, y.lb = y.lb, y.ub = y.ub, 
            x.lb = x.lb, x.ub = x.ub, degree.max = degree.max, 
            degree.min = degree.min, ksum.cores = ksum.cores, 
            nmulti = nmulti, optim.degree.cores = optim.degree.cores, 
            optim.nmulti.cores = optim.nmulti.cores, penalty.cutoff = penalty.cutoff, 
            penalty.method = penalty.method, poly.raw = poly.raw, 
            verbose = verbose, ...)
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
            beta.hat <- coef(lm.wfit(x = X, y = kernel.bk(y.eval[i], 
                y, h[1], y.lb, y.ub), w = NZD(kernel.bk(x.eval[i], 
                x, h[2], x.lb, x.ub))))
            beta.hat[!is.na(beta.hat)] \%*\% t(cbind(1, predict(X.poly, 
                x.eval[i]))[, !is.na(beta.hat), drop = FALSE])
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
        int.f.seq.pre.neg <- numeric()
        int.f.seq <- numeric()
        int.f.seq.post <- numeric()
        x.eval.unique <- unique(x.eval)
        proper.out <- mclapply(1:length(x.eval.unique), function(j) {
            K <- kernel.bk(x.eval.unique[j], x, h[2], x.lb, x.ub)
            if (degree == 0) {
                f.seq <- as.numeric(mcmapply(function(i) {
                  mean(kernel.bk(y.seq[i], y, h[1], y.lb, y.ub) * 
                    K)/NZD(mean(K))
                }, 1:n.integrate, mc.cores = ksum.cores))
            }
            else {
                X.poly <- poly(x, raw = poly.raw, degree = degree)
                X <- cbind(1, X.poly)
                X.eval <- cbind(1, predict(X.poly, x.eval.unique[j]))
                f.seq <- as.numeric(mcmapply(function(i) {
                  beta.hat <- coef(lm.wfit(x = X, y = kernel.bk(y.seq[i], 
                    y, h[1], y.lb, y.ub), w = NZD(K)))
                  beta.hat[!is.na(beta.hat)] \%*\% t(X.eval[, !is.na(beta.hat), 
                    drop = FALSE])
                }, 1:n.integrate, mc.cores = ksum.cores))
            }
            if (verbose & any(f.yx[x.eval == x.eval.unique[j]] < 
                0)) 
                warning("negative density estimate reset to 0 via option proper=TRUE in bkcde() [degree = ", 
                  degree, ", j = ", length(x.eval == x.eval.unique[j]), 
                  " element(s), h.y = ", round(h[1], 5), ", h.x = ", 
                  round(h[2], 5), "]", immediate. = TRUE)
            int.f.seq.pre.neg[j] <- integrate.trapezoidal(y.seq, 
                f.seq)[length(y.seq)]
            f.seq[f.seq < 0] <- 0
            int.f.seq[j] <- integrate.trapezoidal(y.seq, f.seq)[length(y.seq)]
            int.f.seq.post[j] <- integrate.trapezoidal(y.seq, 
                f.seq/int.f.seq[j])[length(y.seq)]
            foo <- f.yx[x.eval == x.eval.unique[j]]
            foo[foo < 0] <- 0
            foo <- foo/int.f.seq[j]
            f.yx[x.eval == x.eval.unique[j]] <- foo
        }, mc.cores = proper.cores)
        int.f.seq.pre.neg <- mean(int.f.seq.pre.neg)
        int.f.seq <- mean(int.f.seq)
        int.f.seq.post <- mean(int.f.seq.post)
    }
    else {
        int.f.seq.pre.neg <- NA
        int.f.seq <- NA
        int.f.seq.post <- NA
        if (verbose & any(f.yx < 0)) 
            warning("negative density estimate encountered, consider option proper=TRUE in bkcde() [degree = ", 
                degree, ", ", length(f.yx[f.yx < 0]), " element(s), h.y = ", 
                round(h[1], 5), ", h.x = ", round(h[2], 5), "]", 
                immediate. = TRUE)
    }
    return.list <- list(convergence.mat = convergence.mat, convergence.vec = convergence.vec, 
        convergence = convergence, degree.mat = degree.mat, degree.max = degree.max, 
        degree.min = degree.min, degree = degree, f.yx.integral.post = int.f.seq.post, 
        f.yx.integral.pre.neg = int.f.seq.pre.neg, f.yx.integral = int.f.seq, 
        f = f.yx, h.mat = h.mat, h = h, ksum.cores = ksum.cores, 
        optim.degree.cores = optim.degree.cores, optim.nmulti.cores = optim.nmulti.cores, 
        proper.cores = proper.cores, proper = proper, secs.elapsed = as.numeric(difftime(Sys.time(), 
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
