\name{EssDee}
\alias{EssDee}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EssDee(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (y) 
{
    if (any(dim(as.matrix(y)) == 0)) 
        return(0)
    sd.vec <- apply(as.matrix(y), 2, sd)
    IQR.vec <- apply(as.matrix(y), 2, IQR)/qnorm(0.25, lower.tail = F) * 
        2
    mad.vec <- apply(as.matrix(y), 2, mad)
    a <- apply(cbind(sd.vec, IQR.vec, mad.vec), 1, function(x) max(x))
    if (any(a <= 0)) 
        warning(paste("variable ", which(a <= 0), " appears to be constant", 
            sep = ""))
    a <- apply(cbind(sd.vec, IQR.vec, mad.vec), 1, function(x) min(x[x > 
        0]))
    return(a)
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
