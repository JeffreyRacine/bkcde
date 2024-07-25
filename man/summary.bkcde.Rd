\name{summary.bkcde}
\alias{summary.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An S3 Method for Summarizing bkcde Objects
}
\description{
This function provides a summary of the results of a \code{bkcde} object.
}
\usage{
summary.bkcde(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class \code{bkcde}.
}
  \item{\dots}{
  Additional arguments to be passed to \code{summary.bkcde}.
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
        ", x.ub=", object$x.ub, ", degree=", object$degree, ")\n", 
        sep = "")
    cat("\n")
    cat("Number of sample realizations: ", length(object$y), 
        "\n", sep = "")
    cat("Number of evaluation points: ", length(object$y.eval), 
        "\n", sep = "")
    cat("Bandwidths: h.y = ", object$h[1], ", h.x = ", object$h[2], 
        "\n", sep = "")
    cat("Degree of local polynomial: ", object$degree, "\n", 
        sep = "")
    if (!is.null(object$f.yx.integral)) 
        cat("Integral of estimate (prior to adjustment): ", object$f.yx.integral, 
            "\n", sep = "")
    cat("Number of cores used in parallel processing for kernel sum: ", 
        object$ksum.cores, "\n", sep = "")
    cat("Number of cores used in parallel processing for degree selection: ", 
        object$degree.cores, "\n", sep = "")
    cat("Number of cores used in parallel processing for multistart optimization: ", 
        object$nmulti.cores, "\n", sep = "")
    cat("Total number of cores used in parallel processing: ", 
        object$ksum.cores * object$degree.cores * object$nmulti.cores, 
        "\n", sep = "")
    cat("Elapsed time (total): ", formatC(object$secs.elapsed, 
        format = "f", digits = 2), " seconds\n", sep = "")
    cat("Optimization and estimation time: ", formatC(object$secs.estimate + 
        sum(object$secs.optim.mat), format = "f", digits = 2), 
        " seconds\n", sep = "")
    cat("Optimization and estimation time per core: ", formatC((object$secs.estimate + 
        sum(object$secs.optim.mat))/(object$ksum.cores * object$degree.cores * 
        object$nmulti.cores), format = "f", digits = 2), " seconds/core\n", 
        sep = "")
    cat("Parallel efficiency: ", formatC(object$secs.elapsed/(object$secs.estimate + 
        sum(object$secs.optim.mat)), format = "f", digits = 2), 
        " (allow for overhead and blocking, ideal = ", formatC(1/(object$ksum.cores * 
            object$degree.cores * object$nmulti.cores), format = "f", 
            digits = 2), ")\n", sep = "")
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
