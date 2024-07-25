\name{plot.bkcde}
\alias{plot.bkcde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An S3 Method for Plotting Bivariate Kernel Conditional Density Estimates
}
\description{
This function plots bivariate kernel conditional density estimates and optionally their confidence intervals.
}
\usage{
plot.bkcde(x, 
           ci = FALSE, 
           ci.method = c("all", "Pointwise", "Bonferroni", "Simultaneous"), 
           ci.bias.correct = TRUE, 
           alpha = 0.05, B = 9999, 
           plot.cores = NULL, 
           plot = TRUE, 
           sub = NULL, 
           ylim = NULL, 
           ylab = NULL, 
           xlab = NULL, 
           type = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  An object of class \code{bkcde}.
}
  \item{ci}{
  A logical indicating whether to compute and plot confidence intervals.
}
  \item{ci.method}{
  A character string indicating the method to use for computing confidence intervals.  Possible values are \code{"all"}, \code{"Pointwise"}, \code{"Bonferroni"}, and \code{"Simultaneous"}.
}
  \item{ci.bias.correct}{
  A logical indicating whether to bias correct the confidence intervals.
}
  \item{alpha}{
  A numeric value indicating the significance level.
}
  \item{B}{
  An integer indicating the number of bootstrap samples to use for computing confidence intervals.
}
  \item{plot.cores}{
  An integer indicating the number of cores to use for parallel computation of confidence intervals.  If \code{NULL}, the number of cores is automatically detected.
}
  \item{plot}{
  A logical indicating whether to plot the bivariate kernel conditional density estimate or return a data list.
}
  \item{sub}{
  A character string to be used as a subtitle for the plot.
}
  \item{ylim}{
  A numeric vector of length 2 indicating the y-axis limits.
}
  \item{ylab}{
  A character string to be used as the y-axis label.
}
  \item{xlab}{
  A character string to be used as the x-axis label.
}
  \item{type}{
  A character string indicating the type of plot to use.  Possible values are \code{"l"} for lines and \code{"p"} for points.
}
  \item{\dots}{
  Additional arguments to be passed to the \code{plot} function.
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
function (x, ci = FALSE, ci.method = c("all", "Pointwise", "Bonferroni", 
    "Simultaneous"), ci.bias.correct = TRUE, alpha = 0.05, B = 9999, 
    plot.cores = NULL, plot = TRUE, sub = NULL, ylim = NULL, 
    ylab = NULL, xlab = NULL, type = NULL, ...) 
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
    if (!is.null(plot.cores)) 
        if (plot.cores < 1) 
            stop("plot.cores must be at least 1 in plot.bkcde()")
    ci.pw.lb <- ci.pw.ub <- ci.bf.lb <- ci.bf.ub <- ci.sim.lb <- ci.sim.ub <- bias.vec <- NULL
    secs.start <- Sys.time()
    if (ci) {
        suppressPackageStartupMessages(library(parallel))
        if (is.null(plot.cores)) 
            plot.cores <- detectCores()
        boot.mat <- t(mcmapply(function(b) {
            ii <- sample(1:length(x$y), replace = TRUE)
            bkcde(h = x$h, x = x$x[ii], y = x$y[ii], x.eval = x$x.eval, 
                y.eval = x$y.eval, y.lb = x$y.lb, y.ub = x$y.ub, 
                x.lb = x$x.lb, x.ub = x$x.ub, degree = x$degree)$f
        }, 1:B, mc.cores = plot.cores))
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
            if (is.null(ylim)) 
                ylim <- range(c(x$f, ci.pw.lb, ci.pw.ub))
        }
        else if (ci.method == "Bonferroni") {
            if (is.null(ylim)) 
                ylim <- range(c(x$f, ci.bf.lb, ci.bf.ub))
        }
        else if (ci.method == "Simultaneous") {
            if (is.null(ylim)) 
                ylim <- range(c(x$f, ci.pw.lb, ci.pw.ub, ci.bf.lb, 
                  ci.bf.ub))
        }
        else {
            if (is.null(ylim)) 
                ylim <- range(c(x$f, ci.pw.lb, ci.pw.ub, ci.bf.lb, 
                  ci.bf.ub, ci.sim.lb, ci.sim.ub))
        }
    }
    else {
        if (is.null(ylim)) 
            ylim <- range(x$f)
    }
    if (plot) {
        if (is.null(sub)) 
            sub <- paste("(degree = ", x$degree, ", h.y = ", 
                round(x$h[1], 3), ", h.x = ", round(x$h[2], 3), 
                ", n = ", length(x$y), ")", sep = "")
        if (is.null(ylab)) 
            ylab <- "f(y|x)"
        if (is.null(xlab)) 
            xlab <- paste("y|x=", x$x.eval[1], sep = "")
        if (is.null(type)) 
            type <- "l"
        plot(x$y.eval, x$f, sub = sub, ylim = ylim, ylab = ylab, 
            xlab = xlab, type = type, panel.first = grid(lty = 1), 
            ...)
        if (ci & ci.method == "Pointwise") {
            lines(x$y.eval, ci.pw.lb, lty = 2)
            lines(x$y.eval, ci.pw.ub, lty = 2)
            legend("topright", legend = c("Estimated f(y|x)", 
                paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                  sep = "")), lty = c(1, 2), bty = "n")
        }
        else if (ci & ci.method == "Bonferroni") {
            lines(x$y.eval, ci.bf.lb, lty = 2)
            lines(x$y.eval, ci.bf.ub, lty = 2)
            legend("topright", legend = c("Estimated f(y|x)", 
                paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                  sep = "")), lty = c(1, 2), bty = "n")
        }
        else if (ci & ci.method == "Simultaneous") {
            lines(x$y.eval, ci.sim.lb, lty = 2)
            lines(x$y.eval, ci.sim.ub, lty = 2)
            legend("topright", legend = c("Estimated f(y|x)", 
                paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                  sep = "")), lty = c(1, 2), bty = "n")
        }
        else if (ci & ci.method == "all") {
            lines(x$y.eval, ci.pw.lb, lty = 2)
            lines(x$y.eval, ci.pw.ub, lty = 2)
            lines(x$y.eval, ci.sim.lb, lty = 3)
            lines(x$y.eval, ci.sim.ub, lty = 3)
            lines(x$y.eval, ci.bf.lb, lty = 4)
            lines(x$y.eval, ci.bf.ub, lty = 4)
            legend("topright", legend = c("Estimated f(y|x)", 
                paste(100 * (1 - alpha), "\% Pointwise CIs", sep = ""), 
                paste(100 * (1 - alpha), "\% Simultaneous CIs", 
                  sep = ""), paste(100 * (1 - alpha), "\% Bonferroni CIs", 
                  sep = "")), lty = 1:4, bty = "n")
        }
    }
    else {
        return(list(f = x$f, bias.vec = bias.vec, ci.pw.lb = ci.pw.lb, 
            ci.pw.ub = ci.pw.ub, ci.bf.lb = ci.bf.lb, ci.bf.ub = ci.bf.ub, 
            ci.sim.lb = ci.sim.lb, ci.sim.ub = ci.sim.ub, secs.elapsed = as.numeric(difftime(Sys.time(), 
                secs.start, units = "secs")), plot.cores = plot.cores))
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
