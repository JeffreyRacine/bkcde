\name{bkcde.default}
\alias{bkcde.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default bkcde Function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bkcde.default(h = NULL, 
              x = NULL, 
              y = NULL, 
              x.eval = NULL, 
              y.eval = NULL, 
              x.lb = NULL, 
              y.lb = NULL, 
              x.ub = NULL, 
              y.ub = NULL, 
              degree.cores = NULL, 
              degree.max = 5, 
              degree.min = 0, 
              degree = 0, 
              ksum.cores = 1, 
              n.integrate = 1000, 
              nmulti.cores = NULL, 
              nmulti = 5, 
              poly.raw = FALSE, 
              proper = TRUE, 
              verbose = FALSE,
              ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
  a bandwidth vector of length 2, where \code{h[1]} is the bandwidth for the response variable and \code{h[2]} is the bandwidth for the predictor variable
}
  \item{x}{
  a vector of predictor variable values
}
  \item{y}{
  a vector of response variable values
}
  \item{x.eval}{
  a vector of predictor variable values at which to evaluate the density estimate
}
  \item{y.eval}{
  a vector of response variable values at which to evaluate the density estimate
}
  \item{x.lb}{
  the lower bound of the predictor variable
}
  \item{y.lb}{
  the lower bound of the response variable
}
  \item{x.ub}{
  the upper bound of the predictor variable
}
  \item{y.ub}{
  the upper bound of the response variable
}
  \item{degree.cores}{
  the number of cores to use in parallel processing for the degree optimization
}
  \item{degree.max}{
  the maximum degree of the polynomial to use in the density estimate
}
  \item{degree.min}{
  the minimum degree of the polynomial to use in the density estimate
}
  \item{degree}{
  the degree of the polynomial to use in the density estimate
}
  \item{ksum.cores}{
  the number of cores to use in parallel processing for the kernel sum
}
  \item{n.integrate}{
  the number of points to use in the numerical integration of the density estimate
}
  \item{nmulti.cores}{
  the number of cores to use in parallel processing for the multi-dimensional optimization
}
  \item{nmulti}{
  the number of multi-dimensional optimization runs to use in the density estimate
}
  \item{poly.raw}{
  a logical indicating whether to use raw or orthogonal polynomials in the density estimate
}
  \item{proper}{
  a logical indicating whether to normalize the density estimate to integrate to 1
}
  \item{verbose}{
  a logical indicating whether to print innocuous warnings and suggestions
}
  \item{\dots}{
  additional arguments to be passed to the optimization function
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(bkcde)

set.seed(42)
n <- 250
n.grid <- 100

x <- rnorm(n)
x.eval <- 0
x.eval.grid <- rep(x.eval,n.grid)
y <- rnorm(n,mean=x)
y.grid <- seq(min(y),max(y),length=n.grid)
dgp <- dnorm(y.grid,mean=x.eval)

par(mfrow=c(2,2))

bkcde.out <- bkcde(x=x,y=y,
                   x.eval=x.eval.grid,
                   y.eval=y.grid,
                   proper=TRUE)

sqrt(sum((bkcde.out$f - dgp)^2))

summary(bkcde.out)

plot(bkcde.out,ylim=range(dgp,bkcde.out$f),main="Proper (Corrected) BKCDE")
lines(y.grid,dgp,col=3)

## Create a perspective plot, ensure estimate is proper at each unique x value
## using predict() function which will invoke proper=TRUE since it was the
## option used in the function call (default)

x.seq <- seq(min(x),max(x),length=25)
y.seq <- seq(min(y),max(y),length=length(x.seq))
data.mat.mat <- matrix(NA,nrow=length(x.seq),ncol=length(x.seq))
for(i in 1:length(x.seq)){
  data.mat.mat[i,] <- predict(bkcde.out, newdata=data.frame(y=y.seq,x=rep(x.seq[i],length(x.seq))))
}
persp(x=x.seq,y=y.seq,z=data.mat.mat,xlab="x",ylab="y",zlab="f(y|x)",theta=120,phi=25,ticktype="detailed",main="Proper (Corrected) BKCDE")

bkcde.out <- bkcde(x=x,y=y,
                   x.eval=x.eval.grid,
                   y.eval=y.grid,
                   proper=FALSE)

sqrt(sum((bkcde.out$f - dgp)^2))

summary(bkcde.out)

plot(bkcde.out,ylim=range(dgp,bkcde.out$f),main="Improper (Potentially) BKCDE")
lines(y.grid,dgp,col=3)

## Create a perspective plot, but since the estimate is improper, the estimate
## can potentially take on negative values

x.seq <- seq(min(x),max(x),length=25)
y.seq <- seq(min(y),max(y),length=length(x.seq))
data.mat.mat <- matrix(NA,nrow=length(x.seq),ncol=length(x.seq))
for(i in 1:length(x.seq)){
  data.mat.mat[i,] <- predict(bkcde.out, newdata=data.frame(y=y.seq,x=rep(x.seq[i],length(x.seq))))
}
persp(x=x.seq,y=y.seq,z=data.mat.mat,xlab="x",ylab="y",zlab="f(y|x)",theta=120,phi=25,ticktype="detailed",main="Improper (Potentially) BKCDE")
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
