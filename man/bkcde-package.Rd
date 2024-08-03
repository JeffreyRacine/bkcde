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
%%  ~~ Literature or other references for background information ~~
}
\keyword{package}

%% Uncomment below to imitate parts of library(help = bkcde)
%\section{The \file{DESCRIPTION} File}{\packageDESCRIPTION{bkcde}}
%\section{Documentation Index}{\packageIndices{bkcde}}

\seealso{
%%  ~~ Optional links to other man pages, e.g. ~~
%%  ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
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
f.yx <- bkcde(x=x, y=y, proper=FALSE)
summary(f.yx)
plot(f.yx,plot.3D=TRUE,ci=TRUE,proper=FALSE,ci.method="Bonferroni")
}
