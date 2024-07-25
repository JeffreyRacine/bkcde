\name{fitted.bkcde}
\alias{fitted.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An S3 Method to Extract the Fitted Values from a \code{bkcde} Object
}
\description{
This function extracts the fitted values from a \code{bkcde} object.
}
\usage{
fitted.bkcde(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{bkcde}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of fitted values.
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
n <- 100
n.grid <- 100
x <- rnorm(n)
x.eval <- 0
x.eval.grid <- rep(x.eval,n.grid)
y <- rnorm(n,mean=x)
y.grid <- seq(min(y),max(y),length=n.grid)
dgp <- dnorm(y.grid,mean=x.eval)
bkcde.out <- bkcde(x=x,y=y,
                   x.eval=x.eval.grid,
                   y.eval=y.grid,
                   proper=TRUE)
plot(y.grid,fitted(bkcde.out))
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
