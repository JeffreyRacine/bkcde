## Utility functions for bkcde package

## integrate.trapezoidal() computes the cumulative integral at each sample
## realization using the trapezoidal rule and the cumsum() function as we need
## to compute this in a computationally efficient manner to render some
## estimators "proper" (non-negative and integrating to 1). Vectors are paired
## naturally but need not be ordered (the function will take care of that for
## you).  The function also includes a correction term to ensure the integral is
## correct at the boundary. Optimized version supplanted previous
## (microbenchmark indicates 50% computation time using diff etc.)

integrate.trapezoidal <- function(x, y) {
  n <- length(x)
  order.x <- order(x)
  x <- x[order.x]
  y <- y[order.x]
  dx <- diff(x)
  dy <- diff(y)
  cx <- dx[1]
  ca <- dy[1] / cx
  cb <- dy[n-1] / cx
  cf <- cx^2 / 12 * (cb - ca)
  if (!is.finite(cf)) cf <- 0
  int.vec <- c(0, cumsum(dx * (y[-n] + y[-1]) / 2))
  int.vec <- int.vec - cf
  int.vec[order(order.x)]
}

## get.integral.weights() returns the vector of weights such that the integral
## can be computed as a dot product sum(y * w). This allows pre-computation
## for fixed grids, avoiding repeated overhead in optimization loops.

get.integral.weights <- function(x) {
  n <- length(x)
  # Ensure x is sorted for the weights calculation, though we assume the input
  # grid (y.seq) in optim is already sorted.
  if(is.unsorted(x)) x <- sort(x)
  
  dx <- diff(x)
  w <- numeric(n)
  
  # Standard trapezoidal weights: dx[i]/2 for ends, (dx[i-1] + dx[i])/2 for interior
  w[1] <- dx[1]/2
  w[n] <- dx[n-1]/2
  if(n > 2) {
    w[2:(n-1)] <- (dx[1:(n-2)] + dx[2:(n-1)])/2
  }
  
  # Boundary correction terms matching integrate.trapezoidal
  # Correction CF = (cx^2 / 12) * ( (y[n]-y[n-1])/cx - (y[2]-y[1])/cx )
  # We subtract CF from the integral.
  # CF = (cx/12) * y[n] - (cx/12) * y[n-1] - (cx/12) * y[2] + (cx/12) * y[1]
  # So we adjust weights:
  # w[n] -= cx/12
  # w[n-1] += cx/12
  # w[2] += cx/12
  # w[1] -= cx/12
  
  cx <- dx[1]
  c_term <- cx/12
  
  w[1] <- w[1] - c_term
  w[2] <- w[2] + c_term
  w[n-1] <- w[n-1] + c_term
  w[n] <- w[n] - c_term
  
  return(w)
}

## NZD() is the "No Zero Divide" (NZD) function (so e.g., 0/0 = 0) based on
## accepted coding practice for a variety of languages

NZD <- function(a) {
  eps <- .Machine$double.eps
  if (length(a) == 1) {
    if (a >= 0) { if (a < eps) return(eps) } else { if (a > -eps) return(-eps) }
    return(a)
  }
  idx <- which(abs(a) < eps)
  if (length(idx) > 0) a[idx] <- ifelse(a[idx] >= 0, eps, -eps)
  a
}

## NZD_pos() is a faster version of NZD() for non-negative inputs

NZD_pos <- function(a) {
  eps <- .Machine$double.eps
  if (length(a) == 1) return(if (a < eps) eps else a)
  idx <- which(a < eps)
  if (length(idx) > 0) a[idx] <- eps
  a
}

## EssDee() returns a robust measure of spread (it can accept both vectors and
## matrices)

EssDee <- function(y){
  if(any(dim(as.matrix(y)) == 0)) return(0)
  sd.vec <- apply(as.matrix(y),2,sd)
  IQR.vec <- apply(as.matrix(y),2,IQR)/qnorm(.25,lower.tail=F)*2
  mad.vec <- apply(as.matrix(y),2,mad)
  a <- apply(cbind(sd.vec,IQR.vec,mad.vec),1, function(x) max(x))
  if(any(a<=0)) warning(paste("variable ",which(a<=0)," appears to be constant",sep=""))
  a <- apply(cbind(sd.vec,IQR.vec,mad.vec),1, function(x) min(x[x>0]))  
  return(a)
}

## Simple function used to find the mode of the degree vector (may not be
## unique, when non-unique the user must determine which one to be used)

find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

## bkcde.bin.data() pre-bins data for linear-binning optimization
## This is called once before optimization to avoid re-binning on every
## objective function evaluation

bkcde.bin.data <- function(x, y, x.lb, x.ub, y.lb, y.ub, n.binned) {
  x.grid <- seq(x.lb, x.ub, length.out = n.binned)
  y.grid <- seq(y.lb, y.ub, length.out = n.binned)
  delta.x <- x.grid[2] - x.grid[1]
  delta.y <- y.grid[2] - y.grid[1]
  
  x.idx <- pmin(pmax(ceiling((x - x.lb) / delta.x), 1), n.binned)
  y.idx <- pmin(pmax(ceiling((y - y.lb) / delta.y), 1), n.binned)
  counts <- as.matrix(table(factor(x.idx, levels=1:n.binned), 
                            factor(y.idx, levels=1:n.binned)))
  
  active <- which(counts > 0, arr.ind = TRUE)
  return(list(
    x.act = x.grid[active[,1]], 
    y.act = y.grid[active[,2]], 
    w.act = as.numeric(counts[active]),
    n.obs = length(x)
  ))
}

## Allow for progress to be displayed while running in parallel, use pbmcapply()
## instead of mcmapply() and pbmclapply() instead of mclapply() to display a
## progress bar. The progress bar is only displayed if progress=TRUE in the
## function calls.

mcmapply.progress <- function(...,progress=TRUE) {
  if(progress) {
    pbmcmapply(...)
  } else {
    mcmapply(...)
  }
}

mclapply.progress <- function(...,progress=TRUE) {
  if(progress) {
    pbmclapply(...)
  } else {
    mclapply(...)
  }
}

## persp() does not allow ylim and zlim to be set to NULL, so this function
## allows for these options to be passable thereby avoiding unnecessary
## duplication of code in plot.bkcde()

persp.lim <- function(x=x,y=y,z=z,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=NULL,zlim=NULL,...) {
  if(is.null(ylim) & is.null(zlim)) {
    persp(x=x,y=y,z=z,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype=ticktype,...)
  } else if(is.null(ylim) & !is.null(zlim)) {
    persp(x=x,y=y,z=z,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype=ticktype,zlim=zlim,...)
  } else if(!is.null(ylim) & is.null(zlim)) {
    persp(x=x,y=y,z=z,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype=ticktype,ylim=ylim,...)
  } else {
    persp(x=x,y=y,z=z,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype=ticktype,ylim=ylim,zlim=zlim,...)
  }
}
