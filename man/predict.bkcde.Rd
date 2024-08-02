\name{predict.bkcde}
\alias{predict.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.bkcde(object, newdata, proper = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{proper}{
%%     ~~Describe \code{proper} here~~
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
function (object, newdata, proper = NULL, ...) 
{
    if (!inherits(object, "bkcde")) 
        stop("object must be of class bkcde in predict.bkcde()")
    if (!is.data.frame(newdata)) 
        stop("newdata must be a data frame in predict.bkcde()")
    if (!all(names(newdata) \%in\% c("x", "y"))) 
        stop("newdata must contain columns x and y in predict.bkcde()")
    if (!is.null(proper) & !is.logical(proper)) 
        stop("proper must be logical in predict.bkcde()")
    if (is.null(proper)) 
        proper <- object$proper
    return(bkcde(h = object$h, x = object$x, y = object$y, x.eval = newdata$x, 
        y.eval = newdata$y, y.lb = object$y.lb, y.ub = object$y.ub, 
        x.lb = object$x.lb, x.ub = object$x.ub, proper = proper, 
        degree = object$degree, ...)$f)
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
