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
      cv = c("full","sub"),
      degree.max = 3, 
      degree.min = 0, 
      degree = NULL, 
      fitted.cores = 12,
      ksum.cores = 1, 
      n.grid = 10,
      n.integrate = 1000, 
      n.sub = 100,
      nmulti = 3, 
      optim.degree.cores = NULL, 
      optim.nmulti.cores = NULL, 
      penalty.cutoff = .Machine$double.xmin, 
      penalty.method = c("smooth", "constant", "trim"), 
      poly.raw = FALSE, 
      progress = FALSE,
      proper.cores = 12, 
      proper = TRUE, 
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
  \item{cv}{
  the type of cross-validation to use. The default is "full"
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
  the number of grid points to use for evaluation data grid when x.eval and y.eval are not provided. The default is 10
  }
  \item{n.integrate}{
  the number of points to use for numerical integration. The default is 1000
}
\item{n.sub}{
  the number of observations to use when \code{cv="sub"}. The default is 100
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
  a logical value indicating whether to use proper normalization. The default is TRUE
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
%%  ~~ If necessary, more details than the description above ~~
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
Jeffrey S. Racine <racinej@mcmaster.ca>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fast.optim}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}
}
\examples{
library(bkcde)
set.seed(42)
n <- 100
x <- runif(n,-.25,.25)
y <- rbeta(n,1+x,1.5+x)
f.yx <- bkcde(x=x, y=y, proper=FALSE)
summary(f.yx)
par(mfrow=c(2,2))
plot(f.yx,plot.3D=FALSE,ci=TRUE,ci.preplot=TRUE)
plot(f.yx,plot.3D=TRUE,ci=TRUE,ci.preplot=TRUE,expand=1.75)
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
