\name{integrate.trapezoidal}
\alias{integrate.trapezoidal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
integrate.trapezoidal(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x, y) 
{
    n <- length(x)
    rank.x <- rank(x)
    order.x <- order(x)
    y <- y[order.x]
    x <- x[order.x]
    int.vec <- numeric(length(x))
    cx <- x[2] - x[1]
    ca <- (y[2] - y[1])/cx
    cb <- (y[n] - y[n - 1])/cx
    cf <- cx^2/12 * (cb - ca)
    if (!is.finite(cf)) 
        cf <- 0
    int.vec[1] <- 0
    int.vec[2:n] <- cumsum((x[2:n] - x[2:n - 1]) * (y[2:n] + 
        y[2:n - 1])/2)
    return(int.vec[rank.x] - cf)
  }
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
