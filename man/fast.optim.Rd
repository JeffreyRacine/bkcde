\name{fast.optim}
\alias{fast.optim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Optimization of Bandwidth Vector and Polynomial Degree for bkcde When Using Large Data Sets
}
\description{
This function optimizes the bandwidth vector and polynomial degree for the bkcde function. It is designed to be used with large data sets.
}
\usage{
fast.optim(x, y, n.sub = 500, resamples = 11, progress = TRUE, ...)
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
\item{progress}{
A logical value indicating whether a progress bar should be displayed.
}
  \item{\dots}{
Additional arguments to be passed to the bkcde function.
}
}
\details{
Cross-validated selection of the bandwidth vector and polynomial degree is computationally intensive for large data sets. This function uses a computationally attractive algorithm to select an appropriate bandwidth vector and polynomial degree for the \code{\link{bkcde}} function. The algorithm is based on the work of Racine (1993) and is designed to be used with large data sets. The basic idea is that each (unknown) optimal bandwidth can be written as a (unknown) constant rescaled by the scale each variable (\eqn{\sigma_y} or \eqn{\sigma_y}) and by a power of the sample size (\eqn{n^{-1/6}}). The scale for each variable is based on a robust measure of spread (\code{\link{IQR}} rescaled to be 1 at the standard Gaussian). Resampling based on sub-samples of size \eqn{\ll n} is used to approximate the rescaling constants which should be appropriate regardless of the sample size where, for each resample, likelihood cross-validation is used to select the bandwidth and polynomial order. Importantly, bandwidth properties (rates and optimal values) differ with the degree of polynomial so it is not sensible to unconditionally (unconditionally) return median values across all degrees and bandwidths. Therefore, after computing the appropriate bandwidth vector and degree for each resample, we first select the median polynomial order across all resamples then take a robust measure of the \dQuote{typical} vector of bandwidths corresponding to all resamples of the same degree using \code{\link{covMcd()$center}} from the \pkg{robustbase} package. 
}
\value{
\item{h}{
A numeric vector of length 2 containing the optimal bandwidths.
}
\item{degree}{
An integer specifying the optimal polynomial degree.
}
\item{h.degree.mat}{
A matrix containing the rescaled bandwidths and polynomial degrees.
}
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
\code{\link{bkcde}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}
}
\examples{
library(bkcde)
## Generate a million observations from a beta distribution, estimate the density then plot
## the results along with the true density function from which the data was drawn.
set.seed(42)
s1 <- s2 <- 1.25
n <- 1000000
n.grid <- 25
x <- runif(n,0,4)
y <- rbeta(n,s1+x,s2+x)
## Compute the optimal bandwidth vector and polynomial degree then use them in 
## the call to bkcde()
optimal <- fast.optim(x=x,y=y,n.sub=250)
bkcde.fast <- bkcde(h=optimal$h,degree=optimal$degree,x=x,y=y,proper=FALSE,n.grid=n.grid)
summary(bkcde.fast)
par(mfrow=c(1,2))
plot(bkcde.fast,plot.3D.n.grid=n.grid,theta=120,phi=45,shade=0.5,expand=1.25,main="Estimate")
x.seq <- sort(unique(bkcde.fast$x.eval))
y.seq <- sort(unique(bkcde.fast$y.eval))
persp(x.seq,y.seq,
      matrix(dbeta(bkcde.fast$y.eval,s1+bkcde.fast$x.eval,s2+bkcde.fast$x.eval),n.grid,n.grid),
      theta=120,
      phi=45,
      ticktype="detailed",
      xlab="x",
      ylab="y",
      zlab="f(y|x)",
      shade=0.5,
      expand=1.25,
      main="True Density")
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
