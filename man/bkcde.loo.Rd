\name{bkcde.loo}
\alias{bkcde.loo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bkcde.loo(h = NULL, x = NULL, y = NULL, y.lb = NULL, y.ub = NULL, x.lb = NULL, x.ub = NULL, poly.raw = TRUE, degree = 0, ksum.cores = 1, penalize.neg.loo = FALSE)
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
  \item{ksum.cores}{
%%     ~~Describe \code{ksum.cores} here~~
}
  \item{penalize.neg.loo}{
%%     ~~Describe \code{penalize.neg.loo} here~~
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
function (h = NULL, x = NULL, y = NULL, y.lb = NULL, y.ub = NULL, 
    x.lb = NULL, x.ub = NULL, poly.raw = TRUE, degree = 0, ksum.cores = 1, 
    penalize.neg.loo = FALSE) 
{
    if (y.lb >= y.ub) 
        stop("y.lb must be less than y.ub in bkcde.loo()")
    if (x.lb >= x.ub) 
        stop("x.lb must be less than x.ub in bkcde.loo()")
    if (is.null(x)) 
        stop("must provide x in bkcde.loo()")
    if (is.null(y)) 
        stop("must provide y in bkcde.loo()")
    if (!is.logical(poly.raw)) 
        stop("poly.raw must be logical in bkcde.loo()")
    if (degree < 0 | degree >= length(y)) 
        stop("degree must lie in [0,1,...,", length(y) - 1, "] (i.e., [0,1,dots, n-1]) in bkcde.loo()")
    if (degree == 0) {
        f.loo <- as.numeric(mcmapply(function(i) {
            kernel.bk.x <- kernel.bk(x[i], x[-i], h[2], x.lb, 
                x.ub)
            mean(kernel.bk(y[i], y[-i], h[1], y.lb, y.ub) * kernel.bk.x)/NZD(mean(kernel.bk.x))
        }, 1:length(y), mc.cores = ksum.cores))
    }
    else {
        X.poly <- poly(x, raw = poly.raw, degree = degree)
        X <- cbind(1, X.poly)
        f.loo <- as.numeric(mcmapply(function(i) {
            coef(lm.wfit(x = X[-i, , drop = FALSE], y = kernel.bk(y[i], 
                y[-i], h[1], y.lb, y.ub), w = NZD(kernel.bk(x[i], 
                x[-i], h[2], x.lb, x.ub)))) \%*\% t(X[i, , drop = FALSE])
        }, 1:length(y), mc.cores = ksum.cores))
    }
    if (penalize.neg.loo) {
        f.loo[!is.finite(f.loo) | f.loo <= 0] <- .Machine$double.xmin
        return(sum(log(f.loo)))
    }
    else {
        return(mean(log(f.loo[is.finite(f.loo) & f.loo > 0])))
    }
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
