\name{bkcde}
\alias{bkcde}
\alias{bkcde.call}
\alias{bkcde.default}

\title{
Boundary Corrected Polynomial Adaptive Conditional Kernel Density Estimation
}
\description{
bkcde() is a function that estimates the conditional density of a response variable given a predictor variable using a boundary corrected polynomial adaptive kernel density estimator. The function is designed to be fast and memory efficient by using parallel processing and vectorized operations. Key parameter choices such as the order of the polynomial and bandwidths tuned to the data being analyzed are obtained by a data-driven method (likelihood cross-validation). The function is also designed to be flexible by allowing the user to specify the degree of the polynomial and the bandwidths of the kernel functions manually, if so desired. The function is \dQuote{overloaded} with arguments that allow the user to control the number of cores used for parallel processing, the number of grid points used for evaluation, the number of points used for numerical integration, the number of multi-starts used for optimization, the number of resamples used for sub-sampled cross-validation, the method used to penalize negative density values, and whether to use raw or orthogonal polynomials. The function returns an object of class \code{bkcde} that can be used to plot the conditional density, predict the conditional density at new values, and summarize the conditional density.
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
      bwscaling = FALSE,
      cv = c("full","sub"),
      cv.only = FALSE,
      degree.max = 3, 
      degree.min = 0, 
      degree = NULL, 
      fitted.cores = 12,
      ksum.cores = 1, 
      n.grid = 25,
      n.integrate = 1000, 
      n.sub = 300,
      nmulti = 3, 
      optim.degree.cores = NULL, 
      optim.nmulti.cores = NULL, 
      penalty.cutoff = .Machine$double.xmin, 
      penalty.method = c("smooth", "constant", "trim"), 
      poly.raw = FALSE, 
      progress = FALSE,
      proper.cores = 12, 
      proper = FALSE, 
      resamples = 10,
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
\item{bwscaling}{
  a logical value that when set to TRUE the supplied bandwidths are interpreted as \sQuote{scale factors} 
}
  \item{cv}{
  the type of cross-validation to use. The default is "full"
  }
  \item{cv.only}{
  a logical value indicating whether to only use cross-validation to estimate the bandwidths and NOT also estimate the conditional density nor render it proper if specified (used if the ONLY thing to be extracted from the estimated object is the bandwidth vector; plot cannot be called when TRUE). The default is FALSE
  }
  \item{degree.max}{
  the maximum degree of the polynomial searched over. The default is 3
}
  \item{degree.min}{
  the minimum degree of the polynomial searched. The default is 0
}
  \item{degree}{
  the degree of the polynomial if specified manually. The default is NULL
}
  \item{fitted.cores}{
  the number of cores to use for the fitted values. The default is 12
}
  \item{ksum.cores}{
  the number of cores to use for the kernel sum. The default is 1
}
  \item{n.grid}{
  the number of grid points to use for evaluation data grid when x.eval and y.eval are not provided. The default is 25
  }
  \item{n.integrate}{
  the number of points to use for numerical integration. The default is 1000
}
\item{n.sub}{
  the number of observations to use when \code{cv="sub"}. The default is 300
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
  the cutoff value for the penalty used in the log likelihood function for negative density values. 
}
  \item{penalty.method}{
  the method used to penalize negative density values. The default is "smooth"
}
  \item{poly.raw}{
  a logical value indicating whether to use raw or orthogonal polynomials. The default is FALSE (i.e., orthogonal polynomials are the default)
}
\item{progress}{
  a logical value indicating whether to print progress. The default is FALSE
}
  \item{proper.cores}{
  the number of cores to use for the proper normalization. The default is 12
}
  \item{proper}{
  a logical value indicating whether to use proper normalization. The default is FALSE and a warning is issued if any negative values are detected in the estimate
}
\item{resamples}{
  the number of resamples to use for cross-validation when \code{cv="sub"}. The default is 10
  }
  \item{verbose}{
  a logical value indicating whether to print warnings and issue other low-level messages. The default is FALSE
}
  \item{\dots}{
  additional arguments to be passed to the bkcde.optim() function
}
}
\details{
A few words on computational cost and rendering the resultant estimate proper are in order.

When \code{proper=TRUE}, for each unique value of \code{x} in \code{x.eval} the conditional density is normalized so that a) the integral of the conditional density returned is equal to 1 and b) the conditional density returned is non-negative (it may be so regardless of whether \code{proper=TRUE} or \code{proper=FALSE}, and if polynomial \code{degree=0} the estimate is guaranteed to be non-negative, while if \code{y.lb=-Inf} and \code{y.ub=Inf} the estimate is guaranteed to integrate to 1 over its support regardless of the polynomial degree). 

To render an estimate proper, a sequence of uniformly spaced points for \code{y} of length \code{n.integrate=1000} is created (the length of this sequence can be changed by the user, its range is discussed shortly). When finite bounds \eqn{[a,b]} for \code{y}  exist  the sequence extends across this range (i.e., \code{seq(a,b,length=1000)}). Otherwise, the function \code{\link{extendrange}} is invoked with argument \code{f=2} to generate the sequence and \code{ceiling(n.integrate/(2*f+1))} evaluation points (200 by default) will lie inside the sample range, the remaining (800 by default) lying outside of this range (\code{f} denotes the fraction by which the range is extended; see \code{\link{extendrange}} for further details). For each unique value of \code{x} in \code{x.eval} the conditional density is estimated on the uniform grid, negative values set to zero, and the integral of the conditional density is re-computed and used to normalize the density returned. Summaries of the integrals prior to and after normalization are returned in the output object and can be accessed using the \code{summary} function. Note that when a large number of unique values of \code{x} are used for \code{x.eval} (e.g., all sample realizations in a large sample) this procedure can involve significant compute overhead. The setting for \code{proper} can be toggled on or off at any time, and when constructing confidence intervals the user can choose to use the raw estimate or the proper estimate (which may be indistinguishable). The \code{proper} argument is set to \code{FALSE} by default, though a warning is issued should any negative density values be produced (though integration to one is of course not guaranteed and the user ought to confirm that results with and without \code{proper=FALSE} are comparable). The argument \code{progress} can be set to obtain more granular details of progress (set to \code{TRUE} to print progress). The argument \code{verbose} can additionally be set to obtain even more granular details (set to \code{TRUE} to print warnings and other low-level messages). 
}
\value{
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
the conditional density (properly normalized if \code{proper=TRUE} or raw if \code{proper=FALSE})
}
\item{f.unadjusted}{ 
the raw conditional density or \code{NA} if \code{proper=FALSE}
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
...
}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>
}
\note{
Those in need of confidence intervals should use the \code{\link{plot.bkcde}} function (or simply \code{\link{plot}}) with the \code{ci} argument set to TRUE (the \code{\link{predict.bkcde}} function (or simply \code{\link{predict}}) can be used to predict the conditional density at new values but does not provide confidence intervals). Note that both \code{x.eval} and \code{y.eval} can be provided to the \code{\link{plot.bkcde}} function (will provide identical fitted values to using \code{\link{predict.bkcde}}), and the intervals can be saved in an object for further processing using the option \code{plot.behavior="data"} or \code{plot.behavior="plot-data"}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sub.cv}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}
}
\examples{
library(bkcde)
## Example 1 - a small sample of data and a conditional Beta
## distribution
set.seed(42)
n <- 100
x <- runif(n,-.25,.25)
y <- rbeta(n,1+x,1.5+x)
f.yx <- bkcde(x=x, y=y)
par(mfrow=c(2,2),cex=.75)
plot(f.yx,persp=FALSE,ci=TRUE,ci.preplot=TRUE)
plot(f.yx,ci=TRUE,ci.preplot=TRUE,expand=1.75)
summary(f.yx)
\dontrun{
## Example 2 - A nonlinear Gaussian conditional relationship with
## a large sample of data, heteroskedasticity, and sub-sampled
## cross-validation (may take a few minutes to run, bump up to 1e+06
## if you are curious and have a few minutes to go for a coffee)
set.seed(42)
n <- 1e+05
n.grid <- 50
x <- runif(n,-2,2)
y <- rnorm(n,mean=x^3,sd=1+abs(x))
## The default proper=FALSE is fine given the large sample size as we really only 
## need the bandwidths, but we manually set to TRUE for plot()...
f.yx <- bkcde(x=x,y=y,n.grid=n.grid,cv="sub",progress=TRUE)
par(mfrow=c(1,2),cex=.75)
x.seq <- sort(unique(f.yx$x.eval))
y.seq <- sort(unique(f.yx$y.eval))
dgp.mat <- matrix(dnorm(f.yx$y.eval,mean=f.yx$x.eval^3,sd=1+abs(f.yx$x.eval)/2),n.grid,n.grid)
persp(y=y.seq,
      x=x.seq,
      z=dgp.mat,
      xlab="X",
      ylab="Y",
      zlab="f(y|x)",
      zlim=range(dgp.mat,f.yx$f[f.yx$f>0]),
      theta=120,
      phi=45,
      main="True",
      ticktype="detailed",
      expand=0.75)
plot(f.yx,
     zlim=range(dgp.mat,f.yx$f),
     theta=120,
     phi=45,
     n.grid=n.grid,
     main="Estimate",
     progress=TRUE,
     proper=TRUE,
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
