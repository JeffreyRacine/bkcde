\name{summary.bkcde}
\alias{summary.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.bkcde(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (object, ...) 
{
    if (!inherits(object, "bkcde")) 
        stop("object must be of class bkcde in summary.bkcde()")
    cat("Call:\n")
    cat("bkcde(h=", object$h, ", x, y, x.eval, y.eval, y.lb=", 
        object$y.lb, ", y.ub=", object$y.ub, ", x.lb=", object$x.lb, 
        ", x.ub=", object$x.ub, ", degree=", object$degree, ")\n")
    cat("\n")
    cat("Number of sample realizations: ", length(object$y), 
        "\n")
    cat("Number of evaluation points: ", length(object$y.eval), 
        "\n")
    cat("Bandwidths: h.y = ", object$h[1], ", h.x = ", object$h[2], 
        "\n")
    cat("Degree of local polynomial: ", object$degree, "\n")
    cat("Elapsed time: ", object$secs.elapsed, " seconds\n")
    cat("Integral of estimate prior to adjustment: ", object$f.yx.integral, 
        "\n")
    cat("\n")
    invisible()
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
