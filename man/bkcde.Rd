\name{bkcde}
\alias{bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boundary Kernel Polynomial Adaptive Conditional Density Estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bkcde(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
\link{bkcde.default}
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
