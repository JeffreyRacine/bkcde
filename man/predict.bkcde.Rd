\name{predict.bkcde}
\alias{predict.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S3 Predict Method for Boundary Corrected Polynomial Adaptive Kernel Conditional Density Estimation
}
\description{
This function predicts the conditional density of a response variable given a covariate using a boundary corrected polynomial adaptive kernel conditional density estimation model.
}
\usage{
\method{predict}{bkcde}(object, newdata, proper = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{bkcde} returned by the function \code{\link{bkcde}}.
}
  \item{newdata}{
A data frame containing the covariate and response variable for which the conditional density is to be predicted. The data frame must contain columns named \code{x} and \code{y}.
}
  \item{proper}{
A logical value indicating whether the conditional density should be rendered proper. If \code{NULL}, the value of \code{proper} used in the \code{bkcde} function that created the \code{object} is used.
}
  \item{\dots}{
Additional arguments to be passed to the \code{bkcde} function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{f}{
  A vector of predicted conditional density values.
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
\code{\link{bkcde}}, \code{\link{sub.cv}}, \code{\link{plot.bkcde}}, \code{\link{predict.bkcde}}
}
\examples{
library(bkcde)
set.seed(42)
n <- 250
x <- runif(n,-.25,.25)
y <- rbeta(n,1+x,1.5+x)
f.yx <- bkcde(x=x, y=y, proper=FALSE)
x.eval <- c(0,.1)
y.eval <- c(.25,0.5)
## Predicted density values
predict(f.yx, data.frame(x=x.eval, y=y.eval))
## True density values
dbeta(y.eval, 1+x.eval, 1.5+x.eval)
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
