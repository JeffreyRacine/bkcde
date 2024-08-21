\name{plot.bkcde}
\alias{plot.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S3 Plot Method for Boundary Corrected Polynomial Adaptive Kernel Conditional Density Estimation
}
\description{
Plot method for \code{bkcde} objects which also provides for bootstrapped confidence intervals that include i) pointwise, ii) Bonferroni, iii) Simultaneous, and iv) all. Plotted values can be saved for further processing using the argument \code{plot.behavior}. Both 2D and 3D plots are supported.
}
\usage{
\method{plot}{bkcde}(x, 
     B = 3999, 
     alpha = 0.05, 
     ci = FALSE, 
     ci.bias.correct = TRUE, 
     ci.cores = NULL, 
     ci.method = c("Pointwise", "Bonferroni", "Simultaneous", "all"), 
     ci.preplot = TRUE,
     fitted.cores = NULL,
     n.grid = NULL,
     persp = TRUE, 
     phi = NULL, 
     plot.behavior = c("plot", "plot-data", "data"), 
     plot.unadjusted = FALSE,
     progress = FALSE, 
     proper = NULL, 
     proper.cores = NULL, 
     sub = NULL, 
     theta = NULL, 
     type = NULL, 
     x.eval = NULL, 
     xlab = NULL, 
     y.eval = NULL,
     ylab = NULL, 
     ylim = NULL, 
     zlab = NULL, 
     zlim = NULL, 
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  object of class \code{bkcde}
}
  \item{B}{
 number of bootstrap samples
}
  \item{alpha}{
 significance level
}
  \item{ci}{
  logical, whether to plot confidence intervals
}
  \item{ci.bias.correct}{
 logical, whether to bias correct the confidence intervals
}
  \item{ci.cores}{
 number of cores to use for bootstrapping
}
  \item{ci.method}{
 method for confidence intervals
}
\item{ci.preplot}{
 logical, whether to generate an initial plot and display it before bootstrapping confidence intervals
}
\item{fitted.cores}{
 number of cores to use for fitted values
 }
\item{n.grid}{
 number of grid points
 }
  \item{persp}{
 logical, whether to plot in 3D using \code{\link{persp}} (default) or 2D using \code{\link{plot}}
}
  \item{phi}{
 phi angle for 3D plot
}
  \item{plot.behavior}{
 behavior of plot
}
\item{plot.unadjusted}{
 logical, whether to plot unadjusted density
}
  \item{progress}{
 logical, whether to display progress when estimating and bootstrapping
}
  \item{proper}{
 logical, whether to ensure estimated density (and bootstrap) is non-negative and integrates to 1
}
  \item{proper.cores}{
 number of cores to use for proper
}
  \item{sub}{
 subtitle for plot
}
  \item{theta}{
 theta angle for 3D plot
}
  \item{type}{
  (character) type of plot
}
  \item{x.eval}{
 x value(s) to evaluate at
}
  \item{xlab}{
 x label
}
\item{y.eval}{
 y value(s) to evaluate at
 }
  \item{ylab}{
 y label
}
  \item{ylim}{
 y limits
}
  \item{zlab}{
 z label
}
  \item{zlim}{
 z limits
}
  \item{\dots}{
  additional arguments
}
}
\details{
When an object created by a call to \code{bkcde} is to be plotted, it was created using either the option \code{proper=TRUE} or \code{proper=FALSE}. This will be the default in the call to \code{plot.bkcde} unless the user specifies otherwise (in \code{plot.bkcde} the default is \code{proper=NULL}, i.e., a null value). If the user specifies \code{proper=TRUE} or \code{proper=FALSE} in the call to \code{plot.bkcde} this will override whichever was used in the call to \code{bkcde}. If the user specifies \code{ci=TRUE} in the call to \code{plot.bkcde}, the same rule applies.

If \code{proper=TRUE} is specified in the call to \code{plot.bkcde} (or the call to \code{bkcde} and is not overwritten in the call to \code{plot.bkcde}) and \code{ci=FALSE} is specified in the call to \code{plot.bkcde}, then the option \code{plot.unadjusted=TRUE} can be used to plot the unadjusted density estimate along with the adjusted one. This is useful for comparing the two. However, if the option \code{ci=TRUE} is specified then the two density estimates are not plotted as the confidence intervals are based on the object corresponding to the value of \code{proper} used in the call to \code{bkcde} (or \code{plot.bkcde} if specified). The exception is when the option \code{ci=TRUE} is specified and \code{ci.preplot=TRUE} is used (default) along with \code{proper=TRUE}. In this case, the two density estimates are plotted first while the bootstrap confidence intervals are being constructed, however the subsequent plot containing the confidence intervals will only contain the density estimate corresponding to the value of \code{proper}.

\code{plot.bkcde} tries to avoid recomputing unnecessary objects by assessing what is being plotted and what is already available in the \code{bkcde} object. For example, if the user is plotting a 3D plot (default) and does not change options manually in the plot call (e.g., \code{n.grid}, \code{proper}, etc.) then the values used in the call to \code{bkcde} will be used directly. However, if the user changes the number of grid points for the 3D plot, then the estimate will be recomputed using the new grid points. Similarly, if the user changes the value of \code{proper} in the call to \code{plot.bkcde} then the estimate be recomputed using the new value of \code{proper}. This is done to avoid unnecessary recomputation and to speed up the plotting process, as the following example illustrates:

\preformatted{
  set.seed(42)
  n <- 1e+04
  x <- runif(n,-.25,.25)
  y <- rbeta(n,1+x,1.5+x)
  f.yx <- bkcde(h=c(.5,.5), degree=1, x=x, y=y, n.grid=25, proper=TRUE)
  plot(f.yx)
}

For 2D objects (not the default), the same rules apply provided that the values of \code{x.eval} and \code{y.eval} are those appropriate for a 2D object (i.e., all elements of \code{x.eval} must be identical), and \code{n.grid} is set to a value appropriate for a 2D plot, i.e., \code{n.grid=length(x.eval)}. The following example illustrates this:

\preformatted{
  set.seed(42)
  n <- 1e+04
  x <- runif(n,-.25,.25)
  y <- rbeta(n,1+x,1.5+x)
  x.eval <- rep(0,100)
  y.eval <- seq(min(y),max(y),length=100)
  f.yx <- bkcde(h=c(.5,.5), degree=1, x=x, y=y, 
                x.eval=x.eval, y.eval=y.eval, 
                n.grid = length(x.eval),proper=TRUE)
  plot(f.yx,persp=FALSE)
}
  
Otherwise, \code{plot.bkcde} will need to recompute the estimates on the fly which can be time consuming (see the option \code{progress} in such cases if you wish to be informed of which stage computation is at - this is particularly helpful when \code{ci=TRUE} which invokes bootstrapping that can be computationally costly, particularly for larger sample sizes). 

For those interested in alternative plotting frameworks, you might consider using RGL (OpenGL) to create interactive 3D plots that you can rotate, zoom and otherwise manipulate dynamically using your mouse or cursor. The \pkg{rgl} package must be installed to use the following code (click on the object with your mouse and drag to rotate, use the mouse wheel to zoom in and out, and use the right mouse button to pan):

\preformatted{
  library(rgl)
  set.seed(42)
  n <- 1e+06
  x <- runif(n,0,1)
  y <- rnorm(n,sin(2*pi*x))
  f.yx <- bkcde(x=x,y=y,cv="sub",n.grid=50,progress=TRUE)
  z <- matrix(f.yx$f,50,50)
  open3d()
  par3d(windowRect=c(900,100,900+640,100+640))
  view3d(theta = 120, phi = 45, fov = 80)
  persp3d(x=sort(unique(f.yx$x.eval)),
          y=sort(unique(f.yx$y.eval)),
          z=z,
          xlab="x",
          ylab="y",
          zlab="f(y|x)",
          ticktype="detailed",      
          border="red",
          color=topo.colors(1000)[(1000-1)*(z-min(z))/(max(z)-min(z))+1],
          alpha=.7,
          back="lines")
  grid3d(c("x", "y+", "z"))
}


}
\value{
...
}
\references{
...
}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>
}
\note{
Note that the default number of bootstrap replications, grid size, and level of confidence are set to ensure a sufficient number of replications for the Bonferroni confidence intervals. Of course, fewer can be used if desired, particularly if only pointwise intervals are needed.

Note that if you want \dQuote{quick and dirty} plots of the conditional PDF and CDF you can use the following code (the conditional PDF, and conditional mean, are only computed when proper=TRUE):

\preformatted{
  f.yx <- bkcde(x=x,y=y,proper=TRUE)
  par(mfrow=c(1,2))
  persp(unique(f.yx$x.eval),unique(f.yx$y.eval),matrix(f.yx$f,f.yx$n.grid,f.yx$n.grid),
        theta=120,phi=45,xlab="x",ylab="y",zlab="f(y|x)",main="Conditional PDF")
  persp(unique(f.yx$x.eval),unique(f.yx$y.eval),matrix(f.yx$F,f.yx$n.grid,f.yx$n.grid),
        theta=120,phi=45,xlab="x",ylab="y",zlab="F(y|x)",main="Conditional CDF")
}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bkcde}}, \code{\link{sub.cv}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}
}
\examples{
library(bkcde)
## Example 1 - cross-validated bandwidths and polynomial order
set.seed(42)
n <- 100
x <- runif(n,-.25,.25)
y <- rbeta(n,1+x,1.5+x)
f.yx <- bkcde(x=x, y=y)
foo <- plot(f.yx,ci=TRUE,ci.method="all",plot.behavior="plot-data",proper=FALSE)
plot(f.yx,ci=TRUE,proper=FALSE,ci.method="Simultaneous")
\dontrun{
## Example 2 - large sample, manual bandwidths and polynomial order
n <- 100000
x <- rnorm(n)
y <- rnorm(n,mean=x)
f.yx.large <- bkcde(h=c(0.5,1),degree=0,x=x,y=y)
plot(f.yx.large,n.grid=50)
summary(f.yx.large)
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
