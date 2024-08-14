\name{sub.cv}
\alias{sub.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sub-Sampled Optimization of Bandwidth Vector and Polynomial Degree for bkcde When Using Large Data Sets
}
\description{
This function optimizes the bandwidth vector and polynomial degree for the bkcde function. It is designed to be used with large data sets.
}
\usage{
sub.cv(x, y, 
       n.sub = 300, 
       progress = FALSE,
       replace = FALSE,
       resamples = 10, 
       ...)
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
An integer specifying the sub-sample size to use (\eqn{n_{sub}} must be \eqn{\le n}, the length of x and y).
}
\item{progress}{
A logical value indicating whether a progress bar should be displayed.
}
\item{replace}{
A logical value indicating whether the sub-samples should be drawn with replacement from x and y.
}
  \item{resamples}{
An integer specifying the number of resamples to use in the optimization.
}

  \item{\dots}{
Additional arguments to be passed to the bkcde function.
}
}
\details{
Cross-validated selection of the bandwidth vector and polynomial degree is computationally intensive for large data sets. This function uses a computationally attractive algorithm to select an appropriate bandwidth vector and polynomial degree for the \code{\link{bkcde}} function. The algorithm is based on the work of Racine (1993) and is designed to be used with large data sets. The basic idea is that each (unknown) optimal bandwidth can be written as a (unknown) constant rescaled by the scale of each variable (\eqn{\sigma_y} or \eqn{\sigma_y}) and by a power of the sample size (\eqn{n^{-1/6}}). The scale estimate for each variable is based on a robust measure of spread (\code{\link{IQR}} rescaled to be 1 at the standard Gaussian). Resampling based on sub-samples of size \eqn{\ll n} is used to approximate the rescaling constants which should be appropriate regardless of the sample size where, for each resample, likelihood cross-validation is used to select the bandwidth and polynomial order. Importantly, bandwidth properties (rates and optimal values) may differ with the polynomial degree so it is not sensible to return the \dQuote{typical} value taken across all polynomial degrees selected across the resamples. Therefore, after computing the appropriate bandwidth vector and degree for each resample, we first select the modal polynomial order across all resamples and then construct a robust measure of the \dQuote{typical} vector of bandwidths corresponding to all resamples of the same degree if possible (for \code{\link[robustbase:covMcd]{covMcd()$center}} it is the number of resamples having median degree \eqn{>3} using \code{\link[robustbase:covMcd]{covMcd()$center}} from the \pkg{robustbase} package, otherwise the mean bandwidth vector is returned). 
}
\value{
\item{h.median}{
A numeric vector of length 2 containing the optimal bandwidths for the x and y variables via the median criterion (column 1 for y, column 2 for x).
}
\item{h.mean}{
A numeric vector of length 2 containing the optimal bandwidths for the x and y variables via the mean criterion (column 1 for y, column 2 for x).
}
\item{h.covMcd}{
A numeric vector of length 2 containing the optimal bandwidths for the x and y variables via the robust mean criterion (column 1 for y, column 2 for x).
}
\item{h.ml}{
A numeric vector of length 2 containing the optimal bandwidths for the x and y variables via the maximum likelihood criterion (column 1 for y, column 2 for x).
}
\item{degree}{
An integer specifying the optimal polynomial degree.
}
\item{degree.modal.length}{
An integer specifying the number of resamples containing the modal polynomial degree.
}
\item{h.mat}{
A matrix containing the optimal bandwidths for the x and y variables for each resample (column 1 for y, column 2 for x).
}
\item{cv.vec}{
A numeric vector containing the cross-validated likelihood values for each resample.
}
\item{sf.mat}{
A matrix containing the scale factors of the bandwidths for the y and x variables (column 1 for y, column 2 for x).
}
}
\references{
Racine, J.S. (1993), \dQuote{An Efficient Cross-Validation Algorithm For Window Width Selection for Nonparametric Kernel Regression,} Communications in Statistics, October, Volume 22, Issue 4, pages 1107-1114.

Delaigle, A. and I. Gijbels (2004), \dQuote{Bootstrap bandwidth selection in kernel density estimation from a contaminated sample}, Annals of the Institute of Statistical Mathematics, 56, 19-47.
}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>
}
\note{
Note that sample with replacement is ill-advised, but the switch \code{resample=TRUE} can be used to allow for this. The default is \code{resample=FALSE}.
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
n <- 1e+06
n.grid <- 25
x <- runif(n,1,4)
y <- rbeta(n,s1+x,s2+x)
## Compute the sub-sampled bandwidth vector and polynomial degree then use them in 
## the call to bkcde()
scv.out <- sub.cv(x=x,y=y)
## We set proper=FALSE as the numerical integration involved when proper=TRUE can add 
## unnecessarily to the computation time for this DGP
f.yx <- bkcde(h=scv.out$h.median,degree=scv.out$degree,x=x,y=y,proper=FALSE,n.grid=n.grid)
## Create a 3D plot of the estimate and the true density
par(mfrow=c(1,2),cex=.8)
options(scipen=9)
plot(f.yx,
     plot.3D.n.grid=n.grid,
     theta=120,
     phi=45,
     shade=0.5,
     expand=1.25,
     main=paste("Estimate (n = ",n,")",sep=""))
x.seq <- sort(unique(f.yx$x.eval))
y.seq <- sort(unique(f.yx$y.eval))
persp(x.seq,y.seq,
      matrix(dbeta(f.yx$y.eval,s1+f.yx$x.eval,s2+f.yx$x.eval),n.grid,n.grid),
      theta=120,
      phi=45,
      ticktype="detailed",
      xlab="x",
      ylab="y",
      zlab="f(y|x)",
      shade=0.5,
      expand=1.25,
      main=paste("True Density (n = ",n,")",sep=""))
## summary() for the bkcde() results
summary(f.yx)
\dontrun{
## Example with a million observations from a bivariate normal distribution with 
## a nonlinear mean and heteroskedasticity, calling sub.cv() from within bkcde() 
## via cv="sub"
set.seed(42)
n <- 1e+06
n.grid <- 50
x <- runif(n,0,1)
y <- rnorm(n,mean=2*sin(4*pi*x),sd=1+abs(x))
f.yx <- bkcde(x=x,y=y,n.grid=n.grid,cv="sub",proper=FALSE,progress=TRUE)
par(mfrow=c(1,2),cex=.75)
x.seq <- sort(unique(f.yx$x.eval))
y.seq <- sort(unique(f.yx$y.eval))
dgp.mat <- matrix(dnorm(f.yx$y.eval,mean=2*(4*pi*f.yx$x.eval),sd=1+abs(f.yx$x.eval)),n.grid,n.grid)
persp(y=y.seq,
      x=x.seq,
      z=dgp.mat,
      xlab="X",
      ylab="Y",
      zlab="f(y|x)",
      zlim=range(dgp.mat,f.yx$f),
      theta=120,
      phi=45,
      main="True",
      ticktype="detailed",
      expand=0.75)
plot(f.yx,
     zlim=range(dgp.mat,f.yx$f),
     theta=120,
     phi=45,
     plot.3D.n.grid=n.grid,
     main="Estimate",
     progress=TRUE,
     proper=FALSE,
     expand=0.75)
summary(f.yx)
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
