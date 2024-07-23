\name{bkcde.optim}
\alias{bkcde.optim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bkcde.optim(x = x, y = y, y.lb = y.lb, y.ub = y.ub, x.lb = x.lb, x.ub = x.ub, poly.raw = poly.raw, degree.min = degree.min, degree.max = degree.max, nmulti = nmulti, ksum.cores = ksum.cores, degree.cores = degree.cores, nmulti.cores = nmulti.cores, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
  \item{degree.min}{
%%     ~~Describe \code{degree.min} here~~
}
  \item{degree.max}{
%%     ~~Describe \code{degree.max} here~~
}
  \item{nmulti}{
%%     ~~Describe \code{nmulti} here~~
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
function (x = x, y = y, y.lb = y.lb, y.ub = y.ub, x.lb = x.lb, 
    x.ub = x.ub, poly.raw = poly.raw, degree.min = degree.min, 
    degree.max = degree.max, nmulti = nmulti, ksum.cores = ksum.cores, 
    degree.cores = degree.cores, nmulti.cores = nmulti.cores, 
    ...) 
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
    output <- mclapply(degree.min:degree.max, function(p) {
        optim.return <- mclapply(1:nmulti, function(i) {
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
                degree = p, ksum.cores = ksum.cores, lower = lower, 
                upper = upper, method = "L-BFGS-B", control = list(fnscale = -1)))
            optim.return$secs.optim <- st["elapsed"]
            optim.return$degree <- p
            optim.return
        }, mc.cores = nmulti.cores)
        optim.out <- optim.return[[which.max(sapply(optim.return, 
            function(x) x$value))]]
        optim.out$value.vec <- sapply(optim.return, function(x) x$value)
        optim.out$convergence.vec <- sapply(optim.return, function(x) x$convergence)
        optim.out$par.mat <- t(sapply(optim.return, function(x) x$par))
        optim.out$secs.optim <- sapply(optim.return, function(x) x$secs.optim)
        optim.out
    }, mc.cores = degree.cores)
    output[[which.max(sapply(output, function(x) x$value))]]
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
