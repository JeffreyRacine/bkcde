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
\method{plot}{bkcde}(x, 
     B = 3999, 
     alpha = 0.05, 
     ci = FALSE, 
     ci.bias.correct = TRUE, 
     ci.cores = NULL, 
     ci.method = c("Pointwise", "Bonferroni", "Simultaneous", "all"), 
     ci.progress = TRUE, 
     ksum.cores = NULL, 
     phi = NULL, 
     plot.2D.n.grid = 100, 
     plot.2D.y.grid = NULL, 
     plot.3D = FALSE, 
     plot.3D.n.grid = 10, 
     plot.3D.x.grid = NULL, 
     plot.3D.y.grid = NULL, 
     plot.behavior = c("plot", "plot-data", "data"), 
     proper = NULL, 
     proper.cores = NULL, 
     sub = NULL, 
     theta = NULL, 
     type = NULL, 
     x.eval = NULL, 
     xlab = NULL, 
     ylab = NULL, 
     ylim = NULL, 
     zlab = NULL, 
     zlim = NULL, 
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{ci}{
%%     ~~Describe \code{ci} here~~
}
  \item{ci.bias.correct}{
%%     ~~Describe \code{ci.bias.correct} here~~
}
  \item{ci.cores}{
%%     ~~Describe \code{ci.cores} here~~
}
  \item{ci.method}{
%%     ~~Describe \code{ci.method} here~~
}
  \item{ci.progress}{
%%     ~~Describe \code{ci.progress} here~~
}
  \item{ksum.cores}{
%%     ~~Describe \code{ksum.cores} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{plot.2D.n.grid}{
%%     ~~Describe \code{plot.2D.n.grid} here~~
}
  \item{plot.2D.y.grid}{
%%     ~~Describe \code{plot.2D.y.grid} here~~
}
  \item{plot.3D}{
%%     ~~Describe \code{plot.3D} here~~
}
  \item{plot.3D.n.grid}{
%%     ~~Describe \code{plot.3D.n.grid} here~~
}
  \item{plot.3D.x.grid}{
%%     ~~Describe \code{plot.3D.x.grid} here~~
}
  \item{plot.3D.y.grid}{
%%     ~~Describe \code{plot.3D.y.grid} here~~
}
  \item{plot.behavior}{
%%     ~~Describe \code{plot.behavior} here~~
}
  \item{proper}{
%%     ~~Describe \code{proper} here~~
}
  \item{proper.cores}{
%%     ~~Describe \code{proper.cores} here~~
}
  \item{sub}{
%%     ~~Describe \code{sub} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{x.eval}{
%%     ~~Describe \code{x.eval} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{zlab}{
%%     ~~Describe \code{zlab} here~~
}
  \item{zlim}{
%%     ~~Describe \code{zlim} here~~
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
function (x, B = 3999, alpha = 0.05, ci = FALSE, ci.bias.correct = TRUE, 
    ci.cores = NULL, ci.method = c("Pointwise", "Bonferroni", 
        "Simultaneous", "all"), ci.progress = TRUE, ksum.cores = NULL, 
    phi = NULL, plot.2D.n.grid = 100, plot.2D.y.grid = NULL, 
    plot.3D = FALSE, plot.3D.n.grid = 10, plot.3D.x.grid = NULL, 
    plot.3D.y.grid = NULL, plot.behavior = c("plot", "plot-data", 
        "data"), proper = NULL, proper.cores = NULL, sub = NULL, 
    theta = NULL, type = NULL, x.eval = NULL, xlab = NULL, ylab = NULL, 
    ylim = NULL, zlab = NULL, zlim = NULL, ...) 
{
    if (!inherits(x, "bkcde")) 
        stop("x must be of class bkcde in plot.bkcde()")
    if (!is.logical(ci)) 
        stop("ci must be logical in plot.bkcde()")
    if (!is.logical(ci.progress)) 
        stop("ci.progress must be logical in plot.bkcde()")
    ci.method <- match.arg(ci.method)
    plot.behavior <- match.arg(plot.behavior)
    if (alpha <= 0 | alpha >= 1) 
        stop("alpha must lie in (0,1) in plot.bkcde()")
    if (B < 1) 
        stop("B must be at least 1 in plot.bkcde()")
    if (!is.null(ci.cores)) 
        if (ci.cores < 1) 
            stop("ci.cores must be at least 1 in plot.bkcde()")
    ci.pw.lb <- ci.pw.ub <- ci.bf.lb <- ci.bf.ub <- ci.sim.lb <- ci.sim.ub <- bias.vec <- NULL
    if (!is.null(proper) & !is.logical(proper)) 
        stop("proper must be logical in plot.bkcde()")
    if (plot.3D & !is.null(x.eval)) 
        warning("x.eval passed but ignored in plot.bkcde() when plot.3D = TRUE", 
            immediate. = TRUE)
    if (!is.null(plot.3D.x.grid) & !is.null(plot.3D.y.grid) & 
        length(plot.3D.x.grid) != length(plot.3D.y.grid)) 
        stop("length of plot.3D.x.grid must be equal to length of plot.3D.y.grid in plot.bkcde()")
    if (is.null(proper)) 
        proper <- x$proper
    if (!plot.3D & is.null(x.eval)) 
        x.eval <- median(x$x.eval)
    if (is.null(ksum.cores)) {
        ksum.cores <- x$ksum.cores
    }
    else {
        ksum.cores <- ksum.cores
    }
    if (is.null(proper.cores)) {
        proper.cores <- 1
    }
    else {
        proper.cores <- proper.cores
    }
    secs.start <- Sys.time()
    if (plot.3D) {
        if (is.null(plot.3D.x.grid)) {
            x.grid <- seq(min(x$x.eval), max(x$x.eval), length = plot.3D.n.grid)
        }
        else {
            x.grid <- plot.3D.x.grid
            plot.3D.n.grid <- length(x.grid)
        }
        if (is.null(plot.3D.y.grid)) {
            y.grid <- seq(min(x$y.eval), max(x$y.eval), length = plot.3D.n.grid)
        }
        else {
            y.grid <- plot.3D.y.grid
            plot.3D.n.grid <- length(y.grid)
        }
        if (length(unique(x.grid)) == 1) 
            stop("only one unique x.eval value, cannot deploy persp() in plot.bkcde() (perhaps call bkcde() with non-unique x.eval OR provide plot.3D.x.grid?)")
        if (length(unique(y.grid)) == 1) 
            stop("only one unique y.eval value, cannot deploy persp() in plot.bkcde() (perhaps call bkcde() with non-unique y.eval OR plot.3D.y.grid?)")
        data.grid <- expand.grid(x.grid, y.grid)
        x.plot.eval <- data.grid$Var1
        y.plot.eval <- data.grid$Var2
        x.fitted <- predict(x, newdata = data.frame(x = x.plot.eval, 
            y = y.plot.eval), proper = proper, ksum.cores = ksum.cores, 
            proper.cores = proper.cores, ...)
        predict.mat <- matrix(x.fitted, plot.3D.n.grid, plot.3D.n.grid)
        if (plot.behavior != "data") {
            if (is.null(theta)) 
                theta <- 120
            if (is.null(phi)) 
                phi <- 45
            if (is.null(xlab)) 
                xlab <- "x"
            if (is.null(ylab)) 
                ylab <- "y"
            if (is.null(zlab)) 
                zlab <- "f(y|x)"
            persp.lim(x = x.grid, y = y.grid, z = predict.mat, 
                xlab = xlab, ylab = ylab, zlab = zlab, theta = theta, 
                phi = phi, ticktype = "detailed", ylim = ylim, 
                zlim = zlim, ...)
        }
    }
    else if (!plot.3D) {
        predict.mat <- NULL
        if (is.null(plot.2D.y.grid)) {
            y.plot.eval <- y.grid <- seq(min(x$y.eval), max(x$y.eval), 
                length = plot.2D.n.grid)
        }
        else {
            y.plot.eval <- y.grid <- plot.2D.y.grid
            plot.2D.n.grid <- length(y.grid)
        }
        x.plot.eval <- x.grid <- rep(x.eval, length(y.plot.eval))
        x.fitted <- predict(x, newdata = data.frame(x = x.plot.eval, 
            y = y.plot.eval), proper = proper, ksum.cores = ksum.cores, 
            ...)
        if (plot.behavior != "data") {
            if (is.null(sub)) 
                sub <- paste("(degree = ", x$degree, ", h.y = ", 
                  round(x$h[1], 3), ", h.x = ", round(x$h[2], 
                    3), ", n = ", length(x$y), ")", sep = "")
            if (is.null(ylab)) 
                ylab <- "f(y|x)"
            if (is.null(xlab)) 
                xlab <- paste("y|x=", round(x.eval, digits = 2), 
                  sep = "")
            if (is.null(type)) 
                type <- "l"
            plot(y.plot.eval[order(y.plot.eval)], x.fitted[order(y.plot.eval)], 
                sub = sub, ylim = ylim, ylab = ylab, xlab = xlab, 
                type = type, panel.first = grid(lty = 1), ...)
        }
    }
    if (ci) {
        if (is.null(ci.cores)) 
            ci.cores <- detectCores()
        mcmapply.progress <- function(..., progress = TRUE) {
            if (progress) {
                pbmcmapply(...)
            }
            else {
                mcmapply(...)
            }
        }
        boot.mat <- t(mcmapply.progress(function(b) {
            ii <- sample(1:length(x$y), replace = TRUE)
            bkcde(h = x$h, x = x$x[ii], y = x$y[ii], x.eval = x.plot.eval, 
                y.eval = y.plot.eval, y.lb = x$y.lb, y.ub = x$y.ub, 
                x.lb = x$x.lb, x.ub = x$x.ub, proper = proper, 
                proper.cores = proper.cores, ksum.cores = ksum.cores, 
                degree = x$degree)$f
        }, 1:B, mc.cores = ci.cores, progress = ci.progress))
        if (ci.bias.correct) {
            bias.vec <- colMeans(boot.mat) - x.fitted
            boot.mat <- sweep(boot.mat, 2, bias.vec, "-")
            if (proper) 
                boot.mat <- pmax(boot.mat, 0)
        }
        ci.pw.lb <- apply(boot.mat, 2, quantile, probs = alpha/2)
        ci.pw.ub <- apply(boot.mat, 2, quantile, probs = 1 - 
            alpha/2)
        ci.bf.lb <- apply(boot.mat, 2, quantile, probs = alpha/(2 * 
            length(y.plot.eval)))
        ci.bf.ub <- apply(boot.mat, 2, quantile, probs = 1 - 
            alpha/(2 * length(y.plot.eval)))
        ci.SCS <- SCSrank(boot.mat, conf.level = 1 - alpha)$conf.int
        ci.sim.lb <- ci.SCS[, 1]
        ci.sim.ub <- ci.SCS[, 2]
        cat("\r                                                                                             ")
    }
    else {
        if (is.null(ylim)) 
            ylim <- NULL
    }
    if (plot.behavior != "data") {
        if (plot.3D) {
            if (ci.method == "Pointwise") {
                if (is.null(zlim)) 
                  zlim <- range(c(x.fitted, ci.pw.lb, ci.pw.ub))
            }
            else if (ci.method == "Bonferroni") {
                if (is.null(zlim)) 
                  zlim <- range(c(x.fitted, ci.bf.lb, ci.bf.ub))
            }
            else if (ci.method == "Simultaneous") {
                if (is.null(zlim)) 
                  zlim <- range(c(x.fitted, ci.sim.lb, ci.sim.ub))
            }
            else {
                if (is.null(zlim)) 
                  zlim <- range(c(x.fitted, ci.pw.lb, ci.pw.ub, 
                    ci.bf.lb, ci.bf.ub, ci.sim.lb, ci.sim.ub))
            }
            if (ci & ci.method == "Pointwise") {
                persp.lim(x = x.grid, y = y.grid, z = matrix(ci.pw.lb, 
                  plot.3D.n.grid, plot.3D.n.grid), xlab = "", 
                  ylab = "", zlab = "", theta = theta, phi = phi, 
                  ticktype = "detailed", border = "grey", col = NA, 
                  lty = 2, ylim = ylim, zlim = zlim, ...)
                par(new = TRUE)
                persp.lim(x = x.grid, y = y.grid, z = predict.mat, 
                  xlab = xlab, ylab = ylab, zlab = zlab, theta = theta, 
                  phi = phi, ticktype = "detailed", ylim = ylim, 
                  zlim = zlim, ...)
                par(new = TRUE)
                persp.lim(x = x.grid, y = y.grid, z = matrix(ci.pw.ub, 
                  plot.3D.n.grid, plot.3D.n.grid), xlab = "", 
                  ylab = "", zlab = "", theta = theta, phi = phi, 
                  ticktype = "detailed", border = "grey", col = NA, 
                  lty = 2, ylim = ylim, zlim = zlim, ...)
                legend("topright", legend = c("Estimated f(y|x)", 
                  paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                    sep = "")), lty = c(1, 2), bty = "n")
            }
            else if (ci & ci.method == "Bonferroni") {
                persp.lim(x = x.grid, y = y.grid, z = matrix(ci.bf.lb, 
                  plot.3D.n.grid, plot.3D.n.grid), xlab = "", 
                  ylab = "", zlab = "", theta = theta, phi = phi, 
                  ticktype = "detailed", border = "grey", col = NA, 
                  lty = 2, ylim = ylim, zlim = zlim, ...)
                par(new = TRUE)
                persp.lim(x = x.grid, y = y.grid, z = predict.mat, 
                  xlab = xlab, ylab = ylab, zlab = zlab, theta = theta, 
                  phi = phi, ticktype = "detailed", ylim = ylim, 
                  zlim = zlim, ...)
                par(new = TRUE)
                persp.lim(x = x.grid, y = y.grid, z = matrix(ci.bf.ub, 
                  plot.3D.n.grid, plot.3D.n.grid), xlab = "", 
                  ylab = "", zlab = "", theta = theta, phi = phi, 
                  ticktype = "detailed", border = "grey", col = NA, 
                  lty = 2, ylim = ylim, zlim = zlim, ...)
                legend("topright", legend = c("Estimated f(y|x)", 
                  paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                    sep = "")), lty = c(1, 2), bty = "n")
            }
            else if (ci & ci.method == "Simultaneous") {
                persp.lim(x = x.grid, y = y.grid, z = matrix(ci.sim.lb, 
                  plot.3D.n.grid, plot.3D.n.grid), xlab = "", 
                  ylab = "", zlab = "", theta = theta, phi = phi, 
                  ticktype = "detailed", border = "grey", col = NA, 
                  lty = 2, ylim = ylim, zlim = zlim, ...)
                par(new = TRUE)
                persp.lim(x = x.grid, y = y.grid, z = predict.mat, 
                  xlab = xlab, ylab = ylab, zlab = zlab, theta = theta, 
                  phi = phi, ticktype = "detailed", ylim = ylim, 
                  zlim = zlim, ...)
                par(new = TRUE)
                persp.lim(x = x.grid, y = y.grid, z = matrix(ci.sim.ub, 
                  plot.3D.n.grid, plot.3D.n.grid), xlab = "", 
                  ylab = "", zlab = "", theta = theta, phi = phi, 
                  ticktype = "detailed", border = "grey", col = NA, 
                  lty = 2, ylim = ylim, zlim = zlim, ...)
                legend("topright", legend = c("Estimated f(y|x)", 
                  paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                    sep = "")), lty = c(1, 2), bty = "n")
            }
            else if (ci & ci.method == "all") {
                warning("plotting all confidence intervals in plot.bkcde() may be visually overwhelming, ignored (but you can retrieve the ci data via plot=FALSE)", 
                  immediate. = TRUE)
            }
        }
        else {
            if (ci.method == "Pointwise") {
                if (is.null(ylim)) 
                  ylim <- range(c(x.fitted, ci.pw.lb, ci.pw.ub))
            }
            else if (ci.method == "Bonferroni") {
                if (is.null(ylim)) 
                  ylim <- range(c(x.fitted, ci.bf.lb, ci.bf.ub))
            }
            else if (ci.method == "Simultaneous") {
                if (is.null(ylim)) 
                  ylim <- range(c(x.fitted, ci.sim.lb, ci.sim.ub))
            }
            else {
                if (is.null(ylim)) 
                  ylim <- range(c(x.fitted, ci.pw.lb, ci.pw.ub, 
                    ci.bf.lb, ci.bf.ub, ci.sim.lb, ci.sim.ub))
            }
            plot(y.plot.eval[order(y.plot.eval)], x.fitted[order(y.plot.eval)], 
                sub = sub, ylim = ylim, ylab = ylab, xlab = xlab, 
                type = type, panel.first = grid(lty = 1), ...)
            if (ci & ci.method == "Pointwise") {
                lines(y.plot.eval[order(y.plot.eval)], ci.pw.lb[order(y.plot.eval)], 
                  lty = 2)
                lines(y.plot.eval[order(y.plot.eval)], ci.pw.ub[order(y.plot.eval)], 
                  lty = 2)
                legend("topright", legend = c("Estimated f(y|x)", 
                  paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                    sep = "")), lty = c(1, 2), bty = "n")
            }
            else if (ci & ci.method == "Bonferroni") {
                lines(y.plot.eval[order(y.plot.eval)], ci.bf.lb[order(y.plot.eval)], 
                  lty = 2)
                lines(y.plot.eval[order(y.plot.eval)], ci.bf.ub[order(y.plot.eval)], 
                  lty = 2)
                legend("topright", legend = c("Estimated f(y|x)", 
                  paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                    sep = "")), lty = c(1, 2), bty = "n")
            }
            else if (ci & ci.method == "Simultaneous") {
                lines(y.plot.eval[order(y.plot.eval)], ci.sim.lb[order(y.plot.eval)], 
                  lty = 2)
                lines(y.plot.eval[order(y.plot.eval)], ci.sim.ub[order(y.plot.eval)], 
                  lty = 2)
                legend("topright", legend = c("Estimated f(y|x)", 
                  paste(100 * (1 - alpha), "\% ", ci.method, " CIs", 
                    sep = "")), lty = c(1, 2), bty = "n")
            }
            else if (ci & ci.method == "all") {
                lines(y.plot.eval[order(y.plot.eval)], ci.pw.lb[order(y.plot.eval)], 
                  lty = 2)
                lines(y.plot.eval[order(y.plot.eval)], ci.pw.ub[order(y.plot.eval)], 
                  lty = 2)
                lines(y.plot.eval[order(y.plot.eval)], ci.sim.lb[order(y.plot.eval)], 
                  lty = 3)
                lines(y.plot.eval[order(y.plot.eval)], ci.sim.ub[order(y.plot.eval)], 
                  lty = 3)
                lines(y.plot.eval[order(y.plot.eval)], ci.bf.lb[order(y.plot.eval)], 
                  lty = 4)
                lines(y.plot.eval[order(y.plot.eval)], ci.bf.ub[order(y.plot.eval)], 
                  lty = 4)
                legend("topright", legend = c("Estimated f(y|x)", 
                  paste(100 * (1 - alpha), "\% Pointwise CIs", 
                    sep = ""), paste(100 * (1 - alpha), "\% Simultaneous CIs", 
                    sep = ""), paste(100 * (1 - alpha), "\% Bonferroni CIs", 
                    sep = "")), lty = 1:4, bty = "n")
            }
        }
    }
    if (plot.behavior != "plot") {
        return(list(bias.vec = bias.vec, ci.bf.lb = ci.bf.lb, 
            ci.bf.ub = ci.bf.ub, ci.cores = ci.cores, ci.pw.lb = ci.pw.lb, 
            ci.pw.ub = ci.pw.ub, ci.sim.lb = ci.sim.lb, ci.sim.ub = ci.sim.ub, 
            f = x.fitted, proper.cores = proper.cores, secs.elapsed = as.numeric(difftime(Sys.time(), 
                secs.start, units = "secs")), x.grid = x.grid, 
            x = x.plot.eval, y.grid = y.grid))
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
