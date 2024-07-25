\name{integrate.trapezoidal}
\alias{integrate.trapezoidal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative Trapezoidal Integration
}
\description{
Integrates a function using the trapezoidal rule.
}
\usage{
integrate.trapezoidal(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a domain vector (can be sorted or unsorted)
}
  \item{y}{
  a range vector
}
}
\details{
This function computes the integral of a function using the trapezoidal rule. An edge correction is applied to the first and last intervals to account for the fact that the trapezoidal rule is not exact at the edges.
}
\value{
This function returns a vector of the same length as the input vectors, with the integral of the function at each point.
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
x <- seq(0, 2*pi, length=100)
y <- sin(x)
par(mfrow=c(2,1))
plot(x, y, ylab="sin(x)",type="l")
abline(h=0,lty=2)
plot(x,integrate.trapezoidal(x, y),ylab="Cumulative Integral",type="l")
} 
