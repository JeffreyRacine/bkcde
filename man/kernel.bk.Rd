\name{kernel.bk}
\alias{kernel.bk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boundary Kernel Function
}
\description{
The boundary kernel function is a doubly truncated Gaussian kernel function that is used to estimate the density of a random variable that can explicitly take on values in a bounded interval. The bounds can be the empirical support of the data or can be user-specified.
}
\usage{
kernel.bk(x, X, h, a = -Inf, b = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  point at which the kernel is to be evaluated
}
  \item{X}{
  vector of data points
}
  \item{h}{
  bandwidth
}
  \item{a}{
  lower bound of the interval
}
  \item{b}{
  upper bound of the interval
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a vector of length X containing the values of the kernel function at the points in x-X
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
x <- rnorm(10)
bk <- kernel.bk(x[1], x, 1, -1, 1)
}