\name{plot.bkcde}
\alias{plot.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S3 Plot Method for Boundary Corrected Polynomial Adaptive Kernel Conditional Density Estimation
}
\description{
Plot method for \code{bkcde} objects.
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
     ksum.cores = NULL, 
     phi = NULL, 
     plot.2D.n.grid = 100, 
     plot.2D.y.grid = NULL, 
     plot.3D = TRUE, 
     plot.3D.n.grid = 10, 
     plot.3D.x.grid = NULL, 
     plot.3D.y.grid = NULL, 
     plot.behavior = c("plot", "plot-data", "data"), 
     progress = TRUE, 
     proper = NULL, 
     proper.cores = NULL, 
     sub = NULL, 
     theta = NULL, 
     type = NULL, 
     x.eval = NULL, 
     xlab = NULL, 
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
  \item{ksum.cores}{
 number of cores to use for kernel sum
}
  \item{phi}{
 phi angle for 3D plot
}
  \item{plot.2D.n.grid}{
 number of grid points for 2D plot
}
  \item{plot.2D.y.grid}{
 y grid points for 2D plot
}
  \item{plot.3D}{
 logical, whether to plot in 3D
}
  \item{plot.3D.n.grid}{
 number of grid points for 3D plot
}
  \item{plot.3D.x.grid}{
 x grid points for 3D plot
}
  \item{plot.3D.y.grid}{
 y grid points for 3D plot
}
  \item{plot.behavior}{
 behavior of plot
}
  \item{progress}{
 logical, whether to display progress bar for bootstrapping
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
 x value to evaluate at
}
  \item{xlab}{
 x label
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
Jeffrey S. Racine <racinej@mcmaster.ca>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bkcde}}, \code{\link{fast.optim}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}
}
\examples{
library(bkcde)
## Example 1 - cross-validated bandwidths and polynomial order
set.seed(42)
n <- 250
x <- runif(n,-.25,.25)
y <- rbeta(n,1+x,1.5+x)
f.yx <- bkcde(x=x, y=y, proper=FALSE)
foo <- plot(f.yx,ci=TRUE,ci.method="all",plot.behavior="plot-data",proper=FALSE)
plot(f.yx,ci=TRUE,plot.3D=TRUE,proper=FALSE,ci.method="Simultaneous")
## Example 2 - large sample, manual bandwidths and polynomial order
n <- 100000
x <- rnorm(n)
y <- rnorm(n,mean=x)
f.yx.large <- bkcde(h=c(0.5,1),degree=0,x=x,y=y,proper=FALSE)
plot(f.yx.large,plot.3D.n.grid=50)
summary(f.yx.large)
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
