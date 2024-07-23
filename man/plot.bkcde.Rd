\name{plot.bkcde}
\alias{plot.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.bkcde(x, ci = FALSE, ci.method = c("all", "Pointwise", "Bonferroni", "Simultaneous"), ci.bias.correct = TRUE, alpha = 0.05, B = 9999, mc.cores = NULL, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{ci}{
%%     ~~Describe \code{ci} here~~
}
  \item{ci.method}{
%%     ~~Describe \code{ci.method} here~~
}
  \item{ci.bias.correct}{
%%     ~~Describe \code{ci.bias.correct} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{mc.cores}{
%%     ~~Describe \code{mc.cores} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (x, ci = FALSE, ci.method = c("all", "Pointwise", "Bonferroni", 
    "Simultaneous"), ci.bias.correct = TRUE, alpha = 0.05, B = 9999, 
    mc.cores = NULL, plot = TRUE, ...) 
{
    if (!inherits(x, "bkcde")) 
        stop("x must be of class bkcde in plot.bkcde()")
    if (!is.logical(ci)) 
        stop("ci must be logical in plot.bkcde()")
    ci.method <- match.arg(ci.method)
    if (alpha < 0 | alpha > 1) 
        stop("alpha must lie in [0,1] in plot.bkcde()")
    if (B < 1) 
        stop("B must be at least 1 in plot.bkcde()")
    if (!is.null(mc.cores)) 
        if (mc.cores < 1) 
            stop("mc.cores must be at least 1 in plot.bkcde()")
    ci.pw.lb <- ci.pw.ub <- ci.bf.lb <- ci.bf.ub <- ci.sim.lb <- ci.sim.ub <- bias.vec <- NULL
    secs.start <- Sys.time()
    if (ci) {
        suppressPackageStartupMessages(library(parallel))
        if (is.null(mc.cores)) 
            mc.cores <- detectCores()
        boot.mat <- t(mcmapply(function(b) {
            ii <- sample(1:length(x$y), replace = TRUE)
            bkcde(h = x$h, x = x$x[ii], y = x$y[ii], x.eval = x$x.eval, 
                y.eval = x$y.eval, y.lb = x$y.lb, y.ub = x$y.ub, 
                x.lb = x$x.lb, x.ub = x$x.ub, degree = x$degree)$f
        }, 1:B, mc.cores = mc.cores))
        if (ci.bias.correct) {
            bias.vec <- colMeans(boot.mat) - x$f
            boot.mat <- sweep(boot.mat, 2, bias.vec, "-")
        }
        ci.pw.lb <- apply(boot.mat, 2, quantile, probs = alpha/2)
        ci.pw.ub <- apply(boot.mat, 2, quantile, probs = 1 - 
            alpha/2)
        ci.bf.lb <- apply(boot.mat, 2, quantile, probs = alpha/(2 * 
            length(x$y.eval)))
        ci.bf.ub <- apply(boot.mat, 2, quantile, probs = 1 - 
            alpha/(2 * length(x$y.eval)))
        suppressPackageStartupMessages(require(MCPAN))
        ci.SCS <- SCSrank(boot.mat, conf.level = 1 - alpha)$conf.int
        ci.sim.lb <- ci.SCS[, 1]
        ci.sim.ub <- ci.SCS[, 2]
        if (ci.method == "Pointwise") {
            ylim <- range(c(x$f, ci.pw.lb, ci.pw.ub))
        }
        else if (ci.method == "Bonferroni") {
            ylim <- range(c(x$f, ci.bf.lb, ci.bf.ub))
        }
        else if (ci.method == "Simultaneous") {
            ylim <- range(c(x$f, ci.pw.lb, ci.pw.ub, ci.bf.lb, 
                ci.bf.ub))
        }
        else {
            ylim <- range(c(x$f, ci.pw.lb, ci.pw.ub, ci.bf.lb, 
                ci.bf.ub, ci.sim.lb, ci.sim.ub))
        }
    }
    else {
        ylim <- range(x$f)
    }
    if (plot) {
        plot(x$y.eval, x$f, sub = paste("(degree = ", x$degree, 
            ", h.y = ", round(x$h[1], 4), ", h.x = ", round(x$h[2], 
                4), ", n = ", length(x$y), ")", sep = ""), ylim = ylim, 
            ylab = "f(y|x)", xlab = "y|x", type = "l", panel.first = grid(lty = 1))
        if (ci & ci.method == "Pointwise") {
            lines(x$y.eval, ci.pw.lb, lty = 2)
            lines(x$y.eval, ci.pw.ub, lty = 2)
            legend("topright", legend = c("f(y|x)", paste(100 * 
                (1 - alpha), "\% ", ci.method, " CIs", sep = "")), 
                lty = c(1, 2), bty = "n")
        }
        else if (ci & ci.method == "Bonferroni") {
            lines(x$y.eval, ci.bf.lb, lty = 2)
            lines(x$y.eval, ci.bf.ub, lty = 2)
            legend("topright", legend = c("f(y|x)", paste(100 * 
                (1 - alpha), "\% ", ci.method, " CIs", sep = "")), 
                lty = c(1, 2), bty = "n")
        }
        else if (ci & ci.method == "Simultaneous") {
            lines(x$y.eval, ci.sim.lb, lty = 2)
            lines(x$y.eval, ci.sim.ub, lty = 2)
            legend("topright", legend = c("f(y|x)", paste(100 * 
                (1 - alpha), "\% ", ci.method, " CIs", sep = "")), 
                lty = c(1, 2), bty = "n")
        }
        else if (ci & ci.method == "all") {
            lines(x$y.eval, ci.pw.lb, lty = 2)
            lines(x$y.eval, ci.pw.ub, lty = 2)
            lines(x$y.eval, ci.sim.lb, lty = 3)
            lines(x$y.eval, ci.sim.ub, lty = 3)
            lines(x$y.eval, ci.bf.lb, lty = 4)
            lines(x$y.eval, ci.bf.ub, lty = 4)
            legend("topright", legend = c("f(y|x)", paste(100 * 
                (1 - alpha), "\% Pointwise CIs", sep = ""), paste(100 * 
                (1 - alpha), "\% Simultaneous CIs", sep = ""), 
                paste(100 * (1 - alpha), "\% Bonferroni CIs", 
                  sep = "")), lty = 1:4, bty = "n")
        }
    }
    else {
        return(list(f = x$f, bias.vec = bias.vec, ci.pw.lb = ci.pw.lb, 
            ci.pw.ub = ci.pw.ub, ci.bf.lb = ci.bf.lb, ci.bf.ub = ci.bf.ub, 
            ci.sim.lb = ci.sim.lb, ci.sim.ub = ci.sim.ub, secs.elapsed = as.numeric(difftime(Sys.time(), 
                secs.start, units = "secs"))))
    }
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
