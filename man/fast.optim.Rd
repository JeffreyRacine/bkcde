\name{fast.optim}
\alias{fast.optim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Optimization of Bandwidth Vector and Polynomial Degree for bkcde using Large Data Sets
}
\description{
This function optimizes the bandwidth vector and polynomial degree for the bkcde function. It is designed to be used with large data sets.
}
\usage{
fast.optim(x, y, n.sub = 1000, resamples = 10, proper = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector of data values.
}
  \item{y}{
A numeric vector of data values.
}
  \item{n.sub}{
An integer specifying the number of sub-samples to use in the optimization.
}
  \item{resamples}{
An integer specifying the number of resamples to use in the optimization.
}
  \item{proper}{
A logical value indicating whether the density estimate should be proper (not used hence set to FALSE to avoid unnecessary computation).
}
  \item{\dots}{
Additional arguments to be passed to the bkcde function.
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
Racine, J.S. (1993), \dQuote{An Efficient Cross-Validation Algorithm For Window Width Selection for Nonparametric Kernel Regression,} Communications in Statistics, October, Volume 22, Issue 4, pages 1107-1114.
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
library(bkcde)
## Generate some data from a beta distribution with the same parameters as the
## true density, and estimate the density then plot the results
set.seed(42)
s1 <- s2 <- 1
n <- 100000
n.grid <- 25
x <- runif(n,0,4)
y <- rbeta(n,s1+x,s2+x)
## Compute the optimal h and degree
optimal <- fast.optim(x=x,y=y,n.sub=250,resamples=10,proper=FALSE)
bkcde.fast <- bkcde(h=optimal$h,degree=optimal$degree,x=x,y=y,proper=FALSE)
summary(bkcde.fast)
plot(bkcde.fast)
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
