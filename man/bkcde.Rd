\name{bkcde}
\alias{bkcde}
\alias{bkcde.call}
\alias{bkcde.default}

\title{
Boundary Corrected Polynomial Adaptive Conditional Kernel Density Estimation
}
\description{
\code{\link{bkcde}} (\sQuote{Boundary Kernel Conditional Density Estimator}) is a function that estimates the conditional density of a response variable given a predictor variable using a boundary corrected polynomial adaptive kernel density estimator. The function is designed to be fast and memory efficient by using parallel processing and vectorized operations. Key parameter choices such as the order of the polynomial and bandwidths tuned to the data being analyzed are obtained by a data-driven method (likelihood or least-squares cross-validation). The function is also designed to be flexible by allowing the user to specify the degree of the polynomial and the bandwidths of the kernel functions manually, if so desired. The function is \dQuote{overloaded} with arguments that allow the user to control the number of cores used for parallel processing, the number of grid points used for evaluation, the number of points used for numerical integration, the number of multi-starts used for optimization, the number of resamples used for sub-sampled cross-validation, the method used to penalize negative density values, and whether to use raw or orthogonal polynomials. The function returns an object of class \code{\link{bkcde}} that can be used to plot the conditional density, predict the conditional density at new values, and summarize the conditional density.
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
      bwmethod = c("cv.ml","cv.ls"),
      bwscaling = FALSE,
      cv = c("auto","full","sub"),
      cv.auto.threshold = 5000,
      cv.only = FALSE,
      degree.max = 3, 
      degree.min = 0, 
      degree = NULL, 
      fitted.cores = detectCores(),
      n.grid = 25,
      n.integrate = 100, 
      n.sub = 300,
      nmulti = 3, 
      optim.degree.cores = NULL, 
      optim.ksum.cores = 1, 
      optim.nmulti.cores = NULL, 
      optim.sf.y.lb=0.25,
      optim.sf.x.lb=0.5,
      penalty.cutoff = .Machine$double.xmin, 
      penalty.method = c("smooth", "constant", "trim"), 
      poly.raw = FALSE, 
      progress = FALSE,
      proper.cores = detectCores(),
      proper = FALSE, 
      resamples = 10,
      verbose = FALSE, 
      ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
  a bandwidth vector of length 2, where h[1] is the bandwidth of the kernel function for the response variable and h[2] is the bandwidth of the kernel function for the predictor variable. If h is NULL, the function will estimate the optimal bandwidths using the \code{bkcde.optim} function.
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
\item{bwmethod}{
  the cross-validation method to use (likelihood or least-squares cross-validation)
}
\item{bwscaling}{
  a logical value that when set to TRUE the supplied bandwidths are interpreted as \sQuote{scale factors} 
}
  \item{cv}{
  the type of cross-validation to use. The default is \code{"auto"}. The options are \code{"full"} (full sample cross-validation), \code{"sub"} (sub-sampled cross-validation via \code{\link{sub.cv}}), and \code{"auto"} (automatic selection of \code{"full"} or \code{"sub"} based on the sample size). When \code{cv="auto"} and the sample size is less than cv.auto.threshold, \code{"full"} cross-validation is used, otherwise \code{"sub"} cross-validation is used
  }
  \item{cv.auto.threshold}{
  the threshold sample size used to determine whether to use \code{"full"} or \code{"sub"} cross-validation when \code{cv="auto"}. The default is 5000
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
  the number of cores to use for the fitted values. The default is \code{\link[parallel]{detectCores}}
}
  \item{n.grid}{
  the number of grid points to use for evaluation data grid when x.eval and y.eval are not provided. The default is 25
  }
  \item{n.integrate}{
  the number of points to use for numerical integration. The default is 100
}
\item{n.sub}{
  the number of observations to use when \code{cv="sub"}. The default is 300
  }
  \item{nmulti}{
  the number of multi-starts to use for optimization. The default is 3
}
  \item{optim.degree.cores}{
  manually set the number of cores to use for the optimization of the degree and override the default that depends on \code{\link[parallel]{detectCores}}, otherwise set to balance the load across the cores. The default is NULL
}
  \item{optim.ksum.cores}{
  the number of cores to use for the kernel sum when conducting optimization. The default is 1
}
  \item{optim.nmulti.cores}{
  manually set the number of cores to use for the optimization of the multi-starts and override the default that depends on \code{\link[parallel]{detectCores}},  otherwise set to balance the load across the cores. The default is NULL
}
  \item{optim.sf.y.lb}{
  the lower bound for the scale factor for y when conducting search (must be positive)
}
  \item{optim.sf.x.lb}{
  the lower bound for the scale factor for x when conducting search (must be positive)
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
  the number of cores to use for the proper normalization. The default is \code{\link[parallel]{detectCores}}
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
  additional arguments to be passed to the \code{bkcde.optim} function
}
}
\details{
A few words on computational cost may be in order. This function is designed to be fast and memory efficient by using parallel processing and vectorized operations. This function is also intended to be sufficiently flexible to allow the user to do essentially anything they want, but constructed to protect the user from unnecessary computation. 

Data-driven bandwidth selection is achieved through numerical optimization of a delete-one likelihood or least-squares function (cross-validation). By construction, this involves computational order \eqn{O(n^2)} calculations (\eqn{n} denotes the sample size) for each function evaluation (i.e., each \eqn{(h_y,h_x)} pair). Sub-sampling (Racine (1993)) can be invoked to reduce the computational overhead (see \code{\link{sub.cv}} and the option \code{cv="sub"}) allowing the user to conduct data-driven bandwidth selection for very large samples (say \code{1e+06} or higher). 

Then there is the question of the default function evaluation points \code{x.eval} and \code{y.eval}. One might be tempted to use all sample realizations (which you can by calling \code{bkcde(x=x,y=y,x.eval=x,y.eval=y)} but if the user simply wishes to plot the density, this can result in substantial wasted computation. Furthermore (more details below), if the option \code{proper=TRUE} is invoked this can lead to substantial wasted computational overhead when evaluating the estimated density on all sample realizations. The default is set to evaluate the density on a uniform grid of \code{n.grid=25} points for \code{x.eval} and \code{y.eval} that span the range of the sample realizations resulting in an evaluation grid of \eqn{25\times 25=625} points by default, which can be modified by the user. The goal of avoiding wasted computation is also why the default is to use  sub-sampled cross-validation when the sample size exceeds \code{cv.auto.threshold=5000} (the user can of course override this default).

Next, a few words on computational cost when rendering the resultant estimate proper may be in order (the default is \code{proper=FALSE}). When \code{proper=TRUE}, for each unique value of \code{x} in \code{x.eval} the conditional density is normalized so that a) the integral of the conditional density returned is equal to 1 and b) the conditional density returned is non-negative (it may be so regardless of whether \code{proper=TRUE} or \code{proper=FALSE}, and if polynomial \code{degree=0} the estimate is guaranteed to be non-negative, while if \code{y.lb=-Inf} and \code{y.ub=Inf} the estimate is guaranteed to integrate to 1 over its support regardless of the polynomial degree). 

To render an estimate proper, numerical integration becomes necessary hence a sequence of uniformly spaced points for \code{y} of length \code{n.integrate=1000} is automatically generated (the integral is computed over this sequence, and the length of this sequence can be changed by the user - its range is discussed shortly). When finite bounds \eqn{[a,b]} for \code{y}  exist  the sequence extends across this range (i.e., \code{seq(a,b,length=1000)}). Otherwise, the function \code{\link{extendrange}} is invoked with argument \code{f=2} to generate the sequence and \code{ceiling(n.integrate/(2*f+1))} evaluation points (200 by default) will lie inside the sample range, the remaining (800 by default) lying outside of this range (\code{f} denotes the fraction by which the range is extended; see \code{\link{extendrange}} for further details). For each unique value of \code{x} in \code{x.eval} the conditional density is estimated on the uniform grid, negative values set to zero, and the integral of the conditional density is re-computed and used to normalize the conditional density returned ensuring it is non-negative possessing total probability one. Summaries of the integrals prior to and after normalization are returned in the output object and can be accessed using the \code{\link{summary}} function. Note that when a large number of unique values of \code{x} are used for \code{x.eval} (e.g., all sample realizations in a large sample) this procedure can involve significant computational overhead. The setting for \code{proper} can be toggled on or off at any time, and when constructing confidence intervals the user can choose to use the raw estimate or the proper estimate (which may be indistinguishable). The \code{proper} argument is set to \code{FALSE} by default, though a warning is always issued should any negative density values be produced (though integration to one is of course not guaranteed and the user ought to confirm that results with and without \code{proper=FALSE} are comparable if choosing the latter). The argument \code{progress=TRUE} can be set to obtain more granular details of progress. The argument \code{verbose=TRUE} can additionally be set to obtain even more granular details. 

In addition to computing the conditional density, the conditional mean and cumulative distribution function are also computed. When \code{proper=TRUE} this involves no additional computational cost since they are constructed from the same set of numerical integrals required to render the conditional density proper. The estimated conditional density is denoted by \eqn{\hat f(y|x)}, and when \code{proper=TRUE} the conditional mean is computed as \eqn{\hat E[Y|X=x] = \int y \hat f(y|x) dy} and the cumulative distribution function is computed as \eqn{\hat F[Y=y|X=x] = \int_{-\infty}^y \hat f(t|x) dt} for all \eqn{x} in \code{x.eval} (so if \code{x.eval} and \code{y.eval} are set to the sample realizations \code{x} these objects can be computed for each sample realization \eqn{(x_i,y_i),i=1,\dots,n}). When \code{proper=FALSE} they are computed using the local polynomial smooth for the sample realizations \code{y} and the doubly truncated CDF kernel for the mean and cumulative conditional density, respectively.  The conditional mean and cumulative distribution function are returned in the output object \code{foo} and can be accessed as \code{foo$g} and \code{foo$F}, respectively. Note that both the conditional mean and distribution use the same smoothing parameters as for the conditional density, which will not be optimal in the square error sense, however, they will all be analytic counterparts to one another (which would not occur if different bandwidths and polynomial orders were chosen separately for each object). See the demos named \code{beta} and \code{normal} for some illustrations involving nonlinear DGPs and large sample sizes (i.e., run the command \code{demo(beta)} or \code{demo(normal)}).

If desired, first derivative estimates for \eqn{g(x)}, \eqn{f(y|x)}, and \eqn{F(y|x)} with respect to \eqn{x} can be obtained by setting the arguments \code{degree.min=1}, \code{degree.max=1}, and \code{poly.raw=TRUE} in the \code{\link{bkcde}} function (the estimates are extracted from the derivative coefficient of the local linear polynomial fit). The first derivative of the conditional mean is computed as \eqn{\frac{\partial}{\partial x} \hat E[Y|X=x] = \int y \frac{\partial}{\partial x} \hat f(y|x) dy} and the first derivative of the cumulative distribution function is computed as \eqn{\frac{\partial}{\partial x} \hat F[Y=y|X=x] = \int_{-\infty}^y \frac{\partial}{\partial x} \hat f(t|x) dt} for all \eqn{x} in \code{x.eval}. The first derivative of the conditional density is computed as \eqn{\frac{\partial}{\partial x} \hat f(y|x)} for all \eqn{x} in \code{x.eval} and \eqn{y} in \code{y.eval}. The first derivative of the conditional mean and cumulative distribution function are returned in the output object \code{foo} and can be accessed as \code{foo$g1} and \code{foo$F1}, respectively. The first derivative of the conditional density is returned in the output object \code{foo} and can be accessed as \code{foo$f1}. The first derivative of the conditional density is computed using the local polynomial smooth for the sample realizations \code{y} and the doubly truncated CDF kernel for the mean and cumulative conditional density, respectively. Be aware that this package is designed to produce optimal conditional density estimates, so the user should be cautious when interpreting derivative results (the derivatives will be undersmoothed relative to the MSE-optimal degree of smoothing that would be appropriate for optimal estimation of derivatives). 

}
\value{
\item{convergence.mat}{
a matrix of convergence values for each multi-start and degree reported by \code{\link{optim}}
}
\item{convergence.vec}{
a vector of convergence values for each multi-start reported by \code{\link{optim}}
}
\item{convergence}{
a scalar convergence value reported by \code{\link{optim}}
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
\item{g}{
the conditional expectation of the response variable given the predictor variable, \eqn{E[Y|X=x]}
}
\item{F}{
the cumulative distribution function of the conditional density \eqn{F[Y=y|X=x]}
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
the conditional density  \eqn{f[Y=y|X=x]} (properly normalized if \code{proper=TRUE} or raw if \code{proper=FALSE})
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
\item{optim.ksum.cores}{
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
a scalar value reported by \code{\link{optim}}
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
Racine, J.S. (1993), \dQuote{An Efficient Cross-Validation Algorithm For Window Width Selection for Nonparametric Kernel Regression,} Communications in Statistics, October, Volume 22, Issue 4, pages 1107-1114.
}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>
}
\note{
Those in need of confidence intervals should use the \code{\link{plot.bkcde}} function (or simply \code{\link{plot}}) with the \code{ci} argument set to TRUE (the \code{\link{predict.bkcde}} function (or simply \code{\link{predict}}) can be used to predict the conditional density at new values but does not provide confidence intervals). Note that both \code{x.eval} and \code{y.eval} can be provided to the \code{\link{plot.bkcde}} function (will provide identical fitted values to using \code{\link{predict.bkcde}}), and the intervals can be saved in an object for further processing using the option \code{plot.behavior="data"} or \code{plot.behavior="plot-data"}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sub.cv}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}, \code{\link[np:npcdens]{npcdens}}, \code{\link[lpcde:lpcde]{lpcde}}
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
f.yx <- bkcde(x=x,y=y,n.grid=n.grid,progress=TRUE)
par(mfrow=c(2,2),cex=.75)
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
      main="True PDF",
      ticktype="detailed",
      expand=0.75)
plot(f.yx,
     zlim=range(dgp.mat,f.yx$f),
     theta=120,
     phi=45,
     n.grid=n.grid,
     main="Estimated PDF",
     progress=TRUE,
     expand=0.75)
persp(y=y.seq,
      x=x.seq,
      z=matrix(f.yx$F,n.grid,n.grid),
      xlab="X",
      ylab="Y",
      zlab="F(y|x)",
      theta=120,
      phi=45,
      main="Estimated CDF",
      ticktype="detailed",
      expand=0.75)  
plot(x.seq,f.yx$g[1:n.grid],
     xlab="X",
     ylab="E[Y|X=x]",
     main="Mean",
     type="l")
lines(x.seq,x.seq^3,col=2,lty=2)
legend("topleft",legend=c("Estimated","True"),col=1:2,lty=1:2,bty="n")
summary(f.yx)
## We want derivatives of g(x), f(y|x), and F(y|x) with respect to x, so we fit
## a model with poly.raw=TRUE and degree.min=1, degree.max=1 thereby generating
## a local polynomial model of degree 1, and the derivative estimates are taken
## from the derivative coefficient of these local models. We take a very large
## sample (10 million observations) to ensure that the derivative estimates are
## accurate. We plot the estimated mean function g(x) and true mean function,
## likewise we plot the estimated conditional density f(y|x) and the estimated
## conditional distribution F(y|x). We also plot the estimated derivative of
## g(x), g1(x), the derivative of f(y|x), f1(y|x), and finally the estimated
## derivative of F(y|x), F1(y|x). We conduct 5 multistarts and 25 resamples
## using sub-sample cross-validation as otherwise the computation time would be
## totally impractical.
library(bkcde)
set.seed(42)
n <- 1e+07
n.grid <- 50
x <- runif(n,0,1)
y <- rnorm(n,mean=2*sin(4*pi*x),sd=1+abs(x))
f.yx <- bkcde(x=x,
              y=y,
              cv="sub",
              n.grid=n.grid,
              n.sub=500,
              nmulti=5,
              resamples=25,
              progress=TRUE,
              poly.raw=TRUE,
              degree.min=1,
              degree.max=1)
par(mfrow=c(2,3))
x.seq <- f.yx$x.eval[1:f.yx$n.grid]
plot(x.seq,f.yx$g[1:f.yx$n.grid],
     ylim=range(f.yx$g[1:f.yx$n.grid],2*sin(4*pi*x.seq)),
     xlab="x",ylab="g(x)",type="l")
lines(x.seq,2*sin(4*pi*x.seq),lty=2,col=2)
legend("topright",legend=c("g(x)","2*sin(4*pi*x)"),lty=1:2,col=1:2,bty="n")
persp(unique(f.yx$x.eval),
      unique(f.yx$y.eval),
      matrix(f.yx$f,f.yx$n.grid,f.yx$n.grid),
      theta=120,phi=45,
      xlab="x",
      ylab="y",
      zlab="f(y|x)",
      ticktype="detailed",
      shade=TRUE)
persp(unique(f.yx$x.eval),
      unique(f.yx$y.eval),
      matrix(f.yx$F,f.yx$n.grid,f.yx$n.grid),
      theta=120,phi=45,
      xlab="x",
      ylab="y",
      zlab="F(y|x)",
      ticktype="detailed",
      shade=TRUE)
plot(x.seq,f.yx$g1[1:f.yx$n.grid],
     ylim=range(f.yx$g1[1:f.yx$n.grid],8*pi*cos(4*pi*x.seq)),
     xlab="x",ylab="g1(x)",type="l")
lines(x.seq,8*pi*cos(4*pi*x.seq),lty=2,col=2)
legend("topright",legend=c("g1(x)","8*pi*cos(4*pi*x)"),lty=1:2,col=1:2,bty="n")
persp(unique(f.yx$x.eval),
      unique(f.yx$y.eval),
      matrix(f.yx$f1,f.yx$n.grid,f.yx$n.grid),
      theta=120,phi=45,
      xlab="x",
      ylab="y",
      zlab="f1(y|x)",
      ticktype="detailed",
      shade=TRUE)
persp(unique(f.yx$x.eval),
      unique(f.yx$y.eval),
      matrix(f.yx$F1,f.yx$n.grid,f.yx$n.grid),
      theta=120,phi=45,
      xlab="x",
      ylab="y",
      zlab="F1(y|x)",
      ticktype="detailed",
      shade=TRUE)
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
