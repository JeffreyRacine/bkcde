\name{bkcde-package}
\alias{bkcde-package}
\docType{package}
\title{
Boundary Corrected Polynomial Adaptive Kernel Conditional Density Estimation
}
\description{
Boundary-corrected polynomial adaptive kernel conditional density estimation
for bounded or partially bounded responses.
}
\details{
This package provides a variety of functions to support estimation of bounded response polynomial adaptive conditional kernel density estimates.
}
\author{
Jeffrey S. Racine
}
\references{
...
}
\keyword{package}

%% Uncomment below to imitate parts of library(help = bkcde)
\seealso{
\code{\link{bkcde}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}
}
\examples{
library(bkcde)
set.seed(42)
n <- 100
x <- runif(n,-.25,.25)
s1 <- 1
s2 <- 1.25
y <- rbeta(n,s1+x,s2+x)
f.yx <- bkcde(h=c(0.15,0.10),
              degree=0,
              x=x,
              y=y,
              fitted.cores=1,
              proper.cores=1)
summary(f.yx)
}
