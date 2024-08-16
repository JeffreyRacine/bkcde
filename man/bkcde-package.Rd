\name{bkcde-package}
\alias{bkcde-package}
\docType{package}
\title{
\packageTitle{bkcde}
}
\description{
\packageDescription{bkcde}
}
\details{
This package provides a variety of functions to support estimation of bounded response polynomial adaptive conditional kernel density estimates.
}
\author{
\packageAuthor{bkcde}

Maintainer: \packageMaintainer{bkcde}
}
\references{
...
}
\keyword{package}

%% Uncomment below to imitate parts of library(help = bkcde)
\section{The \file{DESCRIPTION} File}{\packageDESCRIPTION{bkcde}}
\section{Documentation Index}{\packageIndices{bkcde}}

\seealso{
\code{\link{bkcde}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}
}
\examples{
library(bkcde)
set.seed(42)
n <- 250
x <- runif(n,-.25,.25)
s1 <- 1
s2 <- 1.25
y <- rbeta(n,s1+x,s2+x)
f.yx <- bkcde(x=x, y=y)
summary(f.yx)
plot(f.yx,persp=TRUE,ci=TRUE,ci.method="Bonferroni")
}
