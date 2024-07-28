## Description: This file contains the R functions used in the simulation study
## for implementing the boundary kernel conditional density estimator with
## cross-validated bandwidth selection.  The functions are used in the
## simulation study in the file mc.R, licensed under GPL-3.0-or-later; written
## by Jeffrey Racine <racinej@mcmaster.ca>

## mclapply() and mcmapply() from the parallel package are used throughout
## instead of lapply() and mapply() to allow for multi-core computing, and when
## ksum.cores, degree.cores, or nmulti.cores are set to a value greater than 1,
## the number of cores used in the parallel processing is set to the value of
## the respective argument.  But when these are set to 1, the number of cores
## used in the parallel processing is set to 1, i.e., serial processing occurs
## exactly as if lapply() and mapply() were being used.

## The functions are as follows:

## This function will compute the cumulative integral at each sample realization
## using the trapezoidal rule and the cumsum() function as we need to compute
## this in a computationally efficient manner to render some estimators "proper"
## (non-negative and integrating to 1).

integrate.trapezoidal <- function(x,y) {
  n <- length(x)
  rank.x <- rank(x)
  order.x <- order(x)
  y <- y[order.x]
  x <- x[order.x]
  int.vec <- numeric(length(x))
  ## Use a correction term at the boundary: -cx^2/12*(f'(b)-f'(a)),
  ## check for NaN case
  cx  <- x[2]-x[1]
  ca <- (y[2]-y[1])/cx
  cb <- (y[n]-y[n-1])/cx
  cf <- cx^2/12*(cb-ca)
  if(!is.finite(cf)) cf <- 0
  int.vec[1] <- 0
  int.vec[2:n] <- cumsum((x[2:n]-x[2:n-1])*(y[2:n]+y[2:n-1])/2)
  return(int.vec[rank.x]-cf)
}

## This is the "No Zero Divide" (NZD) function (so e.g., 0/0 = 0) based on
## accepted coding practice for a variety of languages

NZD <- function(a) {
  ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
}

## This function returns a robust measure of spread (it can accept both vectors
## and matrices)

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

## This is the doubly truncated Gaussian boundary kernel function from Racine et
## al 2024

kernel.bk <- function(x,X,h,a=-Inf,b=Inf) {
  ## Checking for bounds involves a bit of overhead (20%), so here we presume a
  ## check is performed outside of this function - make sure this is the case!
  ## ifelse(X < a | X > b, 0, dnorm((x-X)/h)/(h*(pnorm((b-x)/h)-pnorm((a-x)/h))))
  dnorm((x-X)/h)/(h*(pnorm((b-x)/h)-pnorm((a-x)/h)))
}

## This is a likelihood function that supports constant, smooth, and trim
## approaches for dealing with density estimates (delete-one) that may be
## improper and negative in particular. Note we use the smallest non-zero
## normalized floating-point number (a power of the radix, i.e., double.base ^
## double.min.exp, normally 2.225074e-308) as the cutoff value for the penalty.
## There are some slightly smaller numbers possible, i.e. log(4.450148e-324) [1]
##  -744.4401, log(4.450148e-325) = -Inf, while log(.Machine$double.xmin) =
## -708.3964 (testing on Apple Silicon M2 Max), but this gives a very short
## range around zero and otherwise provides a smooth function that, to my way of
## thinking, seems to sensibly handles negative delete-one values.

log.likelihood <- function(delete.one.values,
                           penalty.method=c("smooth","constant","trim"),
                           penalty.cutoff=.Machine$double.xmin,
                           verbose=FALSE,
                           degree=degree,
                           h=h) {
  penalty.method <- match.arg(penalty.method)
  if(penalty.cutoff <= 0) stop("penalty.cutoff must be positive in log.likelihood()")
  likelihood.vec <- numeric(length(delete.one.values))
  cutoff.val <- penalty.cutoff
  log.cutoff <- log(cutoff.val)
  if(penalty.method=="constant") {
    ## My "traditional" method, a constant penalty for negative delete-one
    ## values (so the log likelihood function is non smooth, and all negative
    ## values receive identical penalties)
    likelihood.vec[delete.one.values > cutoff.val] <- log(delete.one.values[delete.one.values > cutoff.val])
    likelihood.vec[delete.one.values <= cutoff.val] <- log.cutoff
    if(verbose & any(0 < delete.one.values & delete.one.values < cutoff.val)) warning("delete-one density lies in constant cutoff zone in log.likelihood() [degree = ",
                                                                                      degree,
                                                                                      ", ",
                                                                                      length(delete.one.values[0 < delete.one.values & delete.one.values < cutoff.val]),
                                                                                      " element(s), h.y = ",
                                                                                      round(h[1],4),
                                                                                      ", h.x = ",
                                                                                      round(h[2],4),
                                                                                      "]",
                                                                                      immediate. = TRUE)
  } else if(penalty.method=="smooth") {
    ## A smooth penalty for negative delete-one values (so the log likelihood
    ## function is smooth except for an extremely narrow range at zero, and
    ## negative values receive a penalty that increases as the value becomes
    ## more negative)
    likelihood.vec[delete.one.values > cutoff.val] <- log(delete.one.values[delete.one.values > cutoff.val])
    likelihood.vec[delete.one.values < -cutoff.val] <- -log(abs(delete.one.values[delete.one.values < -cutoff.val]))+2*log.cutoff
    likelihood.vec[-cutoff.val < delete.one.values & delete.one.values < cutoff.val] <- log.cutoff
    if(verbose & any(-cutoff.val < delete.one.values & delete.one.values < cutoff.val)) warning("delete-one density lies in smooth cutoff zone in log.likelihood() [degree = ",
                                                                                                degree,
                                                                                                ", ",
                                                                                                length(delete.one.values[-cutoff.val < delete.one.values & delete.one.values < cutoff.val]),
                                                                                                " element(s), h.y = ",
                                                                                                round(h[1],4),
                                                                                                ", h.x = ",
                                                                                                round(h[2],4),
                                                                                                "]",
                                                                                                immediate. = TRUE)
  } else if(penalty.method=="trim") {
    ## A trim penalty for negative delete-one values (so the log likelihood
    ## ignores negative values so can be shorter than the vector passed in)
    likelihood.vec[delete.one.values > cutoff.val] <- log(delete.one.values[delete.one.values>cutoff.val])
    likelihood.vec[delete.one.values <= cutoff.val] <- NA
    likelihood.vec <- likelihood.vec[!is.na(likelihood.vec)]
  }
  return(likelihood.vec)
}

## This is the leave-one-out likelihood cross-validation function that supports
## local polynomial estimation of degree p (raw polynomials are the default, but
## orthogonal polynomials can be used as well and appear to provide identical
## results for modest p.max)

bkcde.loo <- function(h=NULL,
                      x=NULL,
                      y=NULL,
                      y.lb=NULL,
                      y.ub=NULL,
                      x.lb=NULL,
                      x.ub=NULL,
                      poly.raw=FALSE,
                      degree=0,
                      ksum.cores=1,
                      penalty.method=NULL,
                      penalty.cutoff=NULL,
                      verbose=FALSE) {
  ## Perform some argument checking
  if(y.lb>=y.ub) stop("y.lb must be less than y.ub in bkcde.loo()")
  if(x.lb>=x.ub) stop("x.lb must be less than x.ub in bkcde.loo()")
  if(is.null(x)) stop("must provide x in bkcde.loo()")
  if(is.null(y)) stop("must provide y in bkcde.loo()")
  if(!is.logical(poly.raw)) stop("poly.raw must be logical in bkcde.loo()")
  if(ksum.cores < 1) stop("ksum.cores must be at least 1 in bkcde.loo()")
  if(is.null(penalty.method)) stop("must provide penalty.method in bkcde.loo()")
  if(is.null(penalty.cutoff)) stop("must provide penalty.cutoff in bkcde.loo()")
  if(degree < 0 | degree >= length(y)) stop("degree must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde.loo()")
  if(degree==0) {
    ## For degree 0 don't invoke the overhead associated with lm.wfit(), just
    ## compute the delete-one estimate \hat f_{-i}(y|x) as efficiently as
    ## possible
    f.loo <- as.numeric(mcmapply(function(i){kernel.bk.x<-kernel.bk(x[i],x[-i],h[2],x.lb,x.ub);mean(kernel.bk(y[i],y[-i],h[1],y.lb,y.ub)*kernel.bk.x)/NZD(mean(kernel.bk.x))},1:length(y),mc.cores=ksum.cores))
  } else {
    X.poly <- poly(x,raw=poly.raw,degree=degree)
    X <- cbind(1,X.poly)
    ## For degree > 0 we use, e.g., lm(y~I(x^2)) and fitted values from the
    ## regression to compute the delete-one estimate \hat f_{-i}(y|x) rather
    ## than the intercept term from lm(y-I(x[i]-X)^2), which produce identical
    ## results for raw polynomials. Note singular design matrices are permitted
    ## but will return NA for coefficients, so we need to check for this (should
    ## not occur with orthogonal polynomials)
    f.loo <- as.numeric(mcmapply(function(i){beta.hat<-coef(lm.wfit(x=X[-i,,drop=FALSE],y=kernel.bk(y[i],y[-i],h[1],y.lb,y.ub),w=NZD(kernel.bk(x[i],x[-i],h[2],x.lb,x.ub))));beta.hat[!is.na(beta.hat)]%*%t(X[i,!is.na(beta.hat), drop = FALSE])},1:length(y),mc.cores=ksum.cores))
  }
  return(sum(log.likelihood(f.loo,penalty.method=penalty.method,penalty.cutoff=penalty.cutoff,verbose=verbose,degree=degree,h=h)))
}

## This function computes the conditional density \hat f(y|x) where, if no
## bandwidth is provided, then likelihood cross-validation is used to select the
## bandwidth via numerical optimization with 5 restarts by default to maximize
## the likelihood and (hopefully) avoid local optima. This function supports
## local polynomial orders [0,1,...,n-1] where n is the number of sample
## realizations (raw polynomials are the default, but orthogonal polynomials can
## be used as well and appear to provide identical results for modest p.max)

bkcde <- function(...) UseMethod("bkcde")

bkcde.default <- function(h=NULL,
                          x=NULL,
                          y=NULL,
                          x.eval=NULL,
                          y.eval=NULL,
                          x.lb=NULL,
                          y.lb=NULL,
                          x.ub=NULL,
                          y.ub=NULL,
                          degree.cores=NULL,
                          degree.max=5,
                          degree.min=0,
                          degree=0,
                          ksum.cores=1,
                          n.integrate=1000,
                          nmulti.cores=NULL,
                          nmulti=5,
                          penalty.method=c("smooth","constant","trim"),
                          penalty.cutoff=.Machine$double.xmin,
                          poly.raw=FALSE,
                          proper=TRUE,
                          verbose=FALSE,
                          ...) {
  ## Perform some argument checking, in this function parallel processing takes
  ## place over the number of multistarts, so ideally the number of cores
  ## requested would be equal to the number of multistarts (this is particularly
  ## useful to avoid local optima in the optimization of the bandwidths)
  if(is.null(x)) stop("must provide x in bkcde()")
  if(is.null(y)) stop("must provide y in bkcde()")
  if(is.null(x.eval)) stop("must provide x.eval in bkcde()")
  if(is.null(y.eval)) stop("must provide y.eval in bkcde()")
  if(is.null(y.lb)) y.lb <- min(y)
  if(is.null(y.ub)) y.ub <- max(y)
  if(is.null(x.lb)) x.lb <- min(x)
  if(is.null(x.ub)) x.ub <- max(x)
  if(any(y<y.lb) | any(y>y.ub)) stop("y must lie in [y.lb,y.ub] in bkcde()")
  if(any(y.eval<y.lb) | any(y.eval>y.ub)) stop("y.eval must lie in [y.lb,y.ub] in bkcde()")
  if(any(x<x.lb) | any(x>x.ub)) stop("x must lie in [x.lb,x.ub] in bkcde()")
  if(any(x.eval<x.lb) | any(x.eval>x.ub)) stop("x.eval must lie in [x.lb,x.ub] in bkcde()")
  if(y.lb>=y.ub) stop("y.lb must be less than y.ub in bkcde()")
  if(x.lb>=x.ub) stop("x.lb must be less than x.ub in bkcde()")
  if(!is.logical(poly.raw)) stop("poly.raw must be logical in bkcde()")
  if(!is.logical(proper)) stop("proper must be logical in bkcde()")
  if(!is.logical(verbose)) stop("verbose must be logical in bkcde()")
  if(nmulti < 1) stop("nmulti must be at least 1 in bkcde()")
  if(n.integrate < 1) stop("n.integrate must be at least 1 in bkcde()")
  if(degree < 0 | degree >= length(y)) stop("degree must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.min < 0 | degree.min >= length(y)) stop("degree.min must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.max < 0 | degree.max >= length(y)) stop("degree.max must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.min > degree.max) stop("degree.min must be <= degree.max in bkcde()")
  if(ksum.cores < 1) stop("ksum.cores must be at least 1 in bkcde()")
  if(is.null(degree.cores)) degree.cores <- degree.max-degree.min+1
  if(is.null(nmulti.cores)) nmulti.cores <- nmulti
  penalty.method <- match.arg(penalty.method)
  if(penalty.cutoff <= 0) stop("penalty.cutoff must be positive in bkcde()")
  secs.start.total <- Sys.time()
  if(is.null(h)) {
    optim.out <- bkcde.optim(x=x,
                             y=y,
                             y.lb=y.lb,
                             y.ub=y.ub,
                             x.lb=x.lb,
                             x.ub=x.ub,
                             poly.raw=poly.raw,
                             degree.min=degree.min,
                             degree.max=degree.max,
                             nmulti=nmulti,
                             ksum.cores=ksum.cores,
                             degree.cores=degree.cores,
                             nmulti.cores=nmulti.cores,
                             penalty.method=penalty.method,
                             penalty.cutoff=penalty.cutoff,
                             verbose=verbose,
                             ...)
    h <- optim.out$par
    h.mat <- optim.out$par.mat
    degree <- optim.out$degree
    degree.mat <- optim.out$degree.mat
    value <- optim.out$value
    value.vec <- optim.out$value.vec
    value.mat <- optim.out$value.mat
    convergence <- optim.out$convergence
    convergence.vec <- optim.out$convergence.vec
    convergence.mat <- optim.out$convergence.mat
    secs.optim <- optim.out$secs.optim
    secs.optim.mat <- optim.out$secs.optim.mat
  } else {
    h.mat <- NULL
    degree.mat <- NULL
    value <- NULL
    value.vec <- NULL
    value.mat <- NULL
    convergence <- NULL
    convergence.vec <- NULL
    convergence.mat <- NULL
    secs.optim <- NULL
    secs.optim.mat <- NULL
  }
  secs.start.estimate <- Sys.time()
  ## Compute the conditional density estimate
  if(degree == 0) {
    ## For degree 0 don't invoke the overhead associated with lm.wfit(), just
    ## compute the estimate \hat f(y|x) as efficiently as possible
    f.yx <- as.numeric(mcmapply(function(i){kernel.bk.x<-kernel.bk(x.eval[i],x,h[2],x.lb,x.ub);mean(kernel.bk(y.eval[i],y,h[1],y.lb,y.ub)*kernel.bk.x)/NZD(mean(kernel.bk.x))},1:length(y.eval),mc.cores=ksum.cores))
  } else {
    ## Choice of raw or orthogonal polynomials
    X.poly <- poly(x,raw=poly.raw,degree=degree)
    X <- cbind(1,X.poly)
    ## For degree > 0 we use, e.g., lm(y~I(x^2)) and fitted values from the
    ## regression to estimate \hat f(y|x) rather than the intercept term from
    ## lm(y-I(x[i]-X)^2), which produce identical results for raw polynomials
    f.yx <- as.numeric(mcmapply(function(i){beta.hat<-coef(lm.wfit(x=X,y=kernel.bk(y.eval[i],y,h[1],y.lb,y.ub),w=NZD(kernel.bk(x.eval[i],x,h[2],x.lb,x.ub))));beta.hat[!is.na(beta.hat)]%*%t(cbind(1,predict(X.poly,x.eval[i]))[,!is.na(beta.hat),drop = FALSE])},1:length(y.eval),mc.cores=ksum.cores))
  }
  if(proper & degree > 0) {
    ## Ensure the estimate is proper - this is peculiar to this implementation
    ## following Cattaneo et al 2023 (i.e., at a scalar evaluation point only).
    ## Note if degree is 0 then the estimate will be proper by definition so
    ## there is no need for an adjustment.
    if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
    if(is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(y.lb,extendrange(y,f=10)[2],length=n.integrate)
    if(!is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(extendrange(y,f=10)[1],y.ub,length=n.integrate)
    if(!is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(extendrange(y,f=10)[1],extendrange(y,f=10)[2],length=n.integrate)
    ## Presume estimation at single X evaluation point, again peculiar to this
    ## implementation
    K <- kernel.bk(x.eval[1],x,h[2],x.lb,x.ub)
    X.poly <- poly(x,raw=poly.raw,degree=degree)
    X <- cbind(1,X.poly)
    X.eval <- cbind(1,predict(X.poly,x.eval[1]))
    ## For degree > 0 we use, e.g., lm(y~I(x^2)) and fitted values from the
    ## regression to estimate \hat f(y|x) rather than the intercept term from
    ## lm(y-I(x[i]-X)^2), which produce identical results for raw polynomials
    f.seq <- as.numeric(mcmapply(function(i){beta.hat<-coef(lm.wfit(x=X,y=kernel.bk(y.seq[i],y,h[1],y.lb,y.ub),w=NZD(K)));beta.hat[!is.na(beta.hat)]%*%t(X.eval[,!is.na(beta.hat),drop = FALSE])},1:n.integrate,mc.cores=ksum.cores))
    ## If proper = TRUE, ensure the final result is proper (i.e., non-negative
    ## and integrates to 1, non-negativity of f.yx is already ensured above)
    if(verbose & any(f.yx < 0)) warning("negative density estimate reset to 0 via option proper=TRUE in bkcde() [degree = ",
                                        degree,
                                        ", ",
                                        length(f.yx[f.yx<0]),
                                        " element(s), h.y = ",
                                        round(h[1],4),
                                        ", h.x = ",
                                        round(h[2],4),
                                        "]",
                                        immediate. = TRUE)
    f.yx[f.yx < 0] <- 0
    f.seq[f.seq < 0] <- 0
    int.f.seq <- integrate.trapezoidal(y.seq,f.seq)[length(y.seq)]
    f.yx <- f.yx/int.f.seq
  } else {
    ## Issue warning if the estimate is improper (here, negative)
    int.f.seq <- NULL
    if(verbose & any(f.yx < 0)) warning("negative density estimate encountered, consider option proper=TRUE in bkcde() [degree = ",
                                        degree,
                                        ", ", 
                                        length(f.yx[f.yx<0]), 
                                        " element(s), h.y = ",
                                        round(h[1],4),
                                        ", h.x = ",
                                        round(h[2],4),
                                        "]",
                                        immediate. = TRUE)
  }
  return.list <- list(convergence.mat=convergence.mat,
                      convergence.vec=convergence.vec,
                      convergence=convergence,
                      degree.cores=degree.cores,
                      degree.mat=degree.mat,
                      degree.max=degree.max,
                      degree.min=degree.min,
                      degree=degree,
                      f.yx.integral=int.f.seq,
                      f=f.yx,
                      h.mat=h.mat,
                      h=h,
                      ksum.cores=ksum.cores,
                      nmulti.cores=nmulti.cores,
                      proper=proper,
                      secs.elapsed=as.numeric(difftime(Sys.time(),secs.start.total,units="secs")),
                      secs.estimate=as.numeric(difftime(Sys.time(),secs.start.estimate,units="secs")),
                      secs.optim.mat=secs.optim.mat,
                      value.mat=value.mat,
                      value.vec=value.vec,
                      value=value,
                      x.eval=x.eval,
                      x.lb=x.lb,
                      x.ub=x.ub,
                      x=x,
                      y.eval=y.eval,
                      y.lb=y.lb,
                      y.ub=y.ub,
                      y=y)
  class(return.list) <- "bkcde"
  return(return.list)
}

## This function conducts numerical optimization for bandwidth selection in
## bkcde() using the optim() function with the L-BFGS-B method which allows box
## constraints, that is each variable can be given a lower and/or upper bound
## (bandwidths must be positive so this is necessary).

bkcde.optim <- function(x=x,
                        y=y,
                        y.lb=y.lb,
                        y.ub=y.ub,
                        x.lb=x.lb,
                        x.ub=x.ub,
                        poly.raw=poly.raw,
                        degree.min=degree.min,
                        degree.max=degree.max,
                        nmulti=nmulti,
                        ksum.cores=ksum.cores,
                        degree.cores=degree.cores,
                        nmulti.cores=nmulti.cores,
                        penalty.method=penalty.method,
                        penalty.cutoff=penalty.cutoff,
                        verbose=verbose,
                        ...) {
  ## Conduct some argument checking
  if(degree.min < 0 | degree.max >= length(y)) stop("degree.min must lie in [0,1,...,",
                                                    length(y)-1,
                                                    "] (i.e., [0,1,dots, n-1]) in bkcde.optim()")
  if(degree.max < 0 | degree.max >= length(y)) stop("degree.max must lie in [0,1,...,",
                                                    length(y)-1,
                                                    "] (i.e., [0,1,dots, n-1]) in bkcde.optim()")
  if(degree.max < degree.min) stop("degree.max must be >= degree.min in bkcde.optim()")
  if(missing(x)) stop("must provide x in bkcde.optim()")
  if(missing(y)) stop("must provide y in bkcde.optim()")
  if(missing(y.lb)) stop("must provide y.lb in bkcde.optim()")
  if(missing(y.ub)) stop("must provide y.ub in bkcde.optim()")    
  if(missing(x.lb)) stop("must provide x.lb in bkcde.optim()")
  if(missing(x.ub)) stop("must provide x.ub in bkcde.optim()")
  if(!is.logical(poly.raw)) stop("poly.raw must be logical in bkcde.optim()")
  ## Get the sample size which we use to initialize the bandwidths using some
  ## common rules of thumb, set search bounds for bandwidths
  n <- length(y)
  lower <- 0.1*c(EssDee(y),EssDee(x))*n^{-1/6}
  upper <- 1000*c(EssDee(y),EssDee(x))
  ## Here we conduct optimization over all models in parallel each having
  ## degree p in [degree.min,degree.max]
  degree.return <- mclapply(degree.min:degree.max, function(p) {
    ## Here we run the optimization for each model over nmulti multistarts in
    ## parallel
    nmulti.return <- mclapply(1:nmulti, function(i) {
      if(i==1) {
        init <- c(EssDee(y),EssDee(x))*n^{-1/6}
      } else {
        init <- runif(2,0.5,5)*c(EssDee(y),EssDee(x))*n^{-1/6}
      }
      st <- system.time(optim.return <- optim(par=init,
                                              fn=bkcde.loo,
                                              x=x,
                                              y=y,
                                              y.lb=y.lb,
                                              y.ub=y.ub,
                                              x.lb=x.lb,
                                              x.ub=x.ub,
                                              poly.raw=poly.raw,
                                              degree=p,
                                              ksum.cores=ksum.cores,
                                              penalty.method=penalty.method,
                                              penalty.cutoff=penalty.cutoff,
                                              verbose=verbose,
                                              lower=lower,
                                              upper=upper,
                                              method="L-BFGS-B",
                                              control=list(fnscale = -1)))
      optim.return$secs.optim <- st["elapsed"]
      optim.return$degree <- p
      optim.return
    },mc.cores = nmulti.cores)
    optim.out <- nmulti.return[[which.max(sapply(nmulti.return, function(x) x$value))]]
    optim.out$value.vec <- sapply(nmulti.return, function(x) x$value)
    optim.out$degree.vec <- sapply(nmulti.return, function(x) x$degree)
    optim.out$convergence.vec <- sapply(nmulti.return, function(x) x$convergence)
    optim.out$secs.optim.vec <- sapply(nmulti.return, function(x) x$secs.optim)
    optim.out
  },mc.cores = degree.cores)
  ## Return object with largest likelihood function over all models and
  ## multistarts, padded with additional information
  output.return <- degree.return[[which.max(sapply(degree.return, function(x) x$value))]]
  output.return$par.mat <- t(sapply(degree.return, function(x) x$par))
  output.return$value.vec <- sapply(degree.return, function(x) x$value)
  output.return$value.mat <- t(sapply(degree.return, function(x) x$value.vec))
  output.return$convergence.vec <- sapply(degree.return, function(x) x$convergence)
  output.return$convergence.mat <- t(sapply(degree.return, function(x) x$convergence.vec))
  output.return$degree.mat <- t(sapply(degree.return, function(x) x$degree.vec))
  output.return$secs.optim <- sapply(degree.return, function(x) x$secs.optim)
  output.return$secs.optim.mat <- t(sapply(degree.return, function(x) x$secs.optim.vec))
  return(output.return)
}

## The following S3 function is used to plot the results of the boundary kernel
## CDE along with bootstrap confidence intervals generated as either pointwise
## or Bonferroni corrected intervals. A handful of options are available,
## including returning the confidence intervals (pointwise, Bonferroni and
## simultaneous) and estimates.

plot.bkcde <- function(x,
                       ci = FALSE, 
                       ci.method = c("all","Pointwise","Bonferroni","Simultaneous"), 
                       ci.bias.correct = TRUE,
                       alpha = 0.05, 
                       B = 9999, 
                       plot.cores = NULL,
                       plot = TRUE,
                       sub = NULL,
                       ylim = NULL,
                       ylab = NULL,
                       xlab = NULL,
                       type = NULL,
                       ...) {
  if(!inherits(x,"bkcde")) stop("x must be of class bkcde in plot.bkcde()")
  if(!is.logical(ci)) stop("ci must be logical in plot.bkcde()")
  ci.method <- match.arg(ci.method)
  if(alpha <= 0 | alpha >= 1) stop("alpha must lie in (0,1) in plot.bkcde()")
  if(B < 1) stop("B must be at least 1 in plot.bkcde()")
  if(!is.null(plot.cores)) if(plot.cores < 1) stop("plot.cores must be at least 1 in plot.bkcde()")
  ci.pw.lb <- ci.pw.ub <- ci.bf.lb <- ci.bf.ub <- ci.sim.lb <- ci.sim.ub <- bias.vec <- NULL
  secs.start <- Sys.time()
  if(ci) {
    if(is.null(plot.cores)) plot.cores <- detectCores()
    boot.mat <- t(mcmapply(function(b){
      ii <- sample(1:length(x$y),replace=TRUE)
      bkcde(h=x$h,
            x=x$x[ii],
            y=x$y[ii],
            x.eval=x$x.eval,
            y.eval=x$y.eval,
            y.lb=x$y.lb,
            y.ub=x$y.ub,
            x.lb=x$x.lb,
            x.ub=x$x.ub,
            proper=x$proper,
            degree=x$degree)$f
    },1:B,mc.cores=plot.cores))
    if(ci.bias.correct) {
      bias.vec <- colMeans(boot.mat) - x$f
      boot.mat <- sweep(boot.mat,2,bias.vec,"-")
    }
    ci.pw.lb <- apply(boot.mat, 2, quantile, probs = alpha / 2)
    ci.pw.ub <- apply(boot.mat, 2, quantile, probs = 1 - alpha / 2)
    ci.bf.lb <- apply(boot.mat, 2, quantile, probs = alpha / (2 * length(x$y.eval)))
    ci.bf.ub <- apply(boot.mat, 2, quantile, probs = 1 - alpha / (2 * length(x$y.eval)))
    ci.SCS <- SCSrank(boot.mat, conf.level=1-alpha)$conf.int
    ci.sim.lb <- ci.SCS[,1]
    ci.sim.ub <- ci.SCS[,2]
    if(ci.method == "Pointwise") {
      if(is.null(ylim)) ylim <-  range(c(x$f,ci.pw.lb,ci.pw.ub))
    } else if(ci.method == "Bonferroni") {
      if(is.null(ylim)) ylim <-  range(c(x$f,ci.bf.lb,ci.bf.ub))
    } else if(ci.method == "Simultaneous") {
      if(is.null(ylim)) ylim <-  range(c(x$f,ci.pw.lb,ci.pw.ub,ci.bf.lb,ci.bf.ub))
    } else {
      if(is.null(ylim)) ylim <-  range(c(x$f,ci.pw.lb,ci.pw.ub,ci.bf.lb,ci.bf.ub,ci.sim.lb,ci.sim.ub))
    }
  } else {
    if(is.null(ylim)) ylim <-  range(x$f)
  }
  if(plot) {
    if(is.null(sub)) sub <- paste("(degree = ",x$degree,", h.y = ",round(x$h[1],3), ", h.x = ",round(x$h[2],3),", n = ",length(x$y),")",sep="")
    if(is.null(ylab)) ylab <- "f(y|x)"
    if(is.null(xlab)) xlab <- paste("y|x=",x$x.eval[1],sep="")
    if(is.null(type)) type <- "l"
    plot(x$y.eval,x$f,
         sub=sub,
         ylim=ylim,
         ylab=ylab,
         xlab=xlab,
         type=type,
         panel.first=grid(lty=1),
         ...)
    if(ci & ci.method == "Pointwise") {
      lines(x$y.eval,ci.pw.lb,lty=2)
      lines(x$y.eval,ci.pw.ub,lty=2)
      legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
    } else if(ci & ci.method == "Bonferroni") {
      lines(x$y.eval,ci.bf.lb,lty=2)
      lines(x$y.eval,ci.bf.ub,lty=2)
      legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
    } else if(ci & ci.method == "Simultaneous") {
      lines(x$y.eval,ci.sim.lb,lty=2)
      lines(x$y.eval,ci.sim.ub,lty=2)
      legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
    } else if(ci & ci.method == "all") {
      lines(x$y.eval,ci.pw.lb,lty=2)
      lines(x$y.eval,ci.pw.ub,lty=2)
      lines(x$y.eval,ci.sim.lb,lty=3)
      lines(x$y.eval,ci.sim.ub,lty=3)
      lines(x$y.eval,ci.bf.lb,lty=4)
      lines(x$y.eval,ci.bf.ub,lty=4)
      legend("topright",
             legend=c("Estimated f(y|x)",
                      paste(100*(1-alpha),"% Pointwise CIs",sep=""),
                      paste(100*(1-alpha),"% Simultaneous CIs",sep=""),
                      paste(100*(1-alpha),"% Bonferroni CIs",sep="")
             ),lty=1:4,bty="n")
    }
  } else {
    return(list(f=x$f,
                bias.vec=bias.vec,
                ci.pw.lb=ci.pw.lb,
                ci.pw.ub=ci.pw.ub,
                ci.bf.lb=ci.bf.lb,
                ci.bf.ub=ci.bf.ub,
                ci.sim.lb=ci.sim.lb,
                ci.sim.ub=ci.sim.ub,
                secs.elapsed=as.numeric(difftime(Sys.time(),secs.start,units="secs")),
                plot.cores=plot.cores))
  }
}

fitted.bkcde <- function(object, ...) {
  if(!inherits(object,"bkcde")) stop("object must be of class bkcde in fitted.bkcde()")
  return(object$f)
}

predict.bkcde <- function(object, newdata, ...) {
  if(!inherits(object,"bkcde")) stop("object must be of class bkcde in predict.bkcde()")
  if(!is.data.frame(newdata)) stop("newdata must be a data frame in predict.bkcde()")
  if(!all(names(newdata) %in% c("x","y"))) stop("newdata must contain columns x and y in predict.bkcde()")
  return(bkcde(h=object$h,
               x=object$x,
               y=object$y,
               x.eval=newdata$x,
               y.eval=newdata$y,
               y.lb=object$y.lb,
               y.ub=object$y.ub,
               x.lb=object$x.lb,
               x.ub=object$x.ub,
               proper=object$proper,
               degree=object$degree)$f)
}

summary.bkcde <- function(object, ...) {
  if(!inherits(object,"bkcde")) stop("object must be of class bkcde in summary.bkcde()")
  cat("Call:\n")
  cat("bkcde(h=",object$h,", x, y, x.eval, y.eval, y.lb=",object$y.lb,", y.ub=",object$y.ub,", x.lb=",object$x.lb,", x.ub=",object$x.ub,", degree=",object$degree,")\n",sep="")
  cat("\n")
  cat("Number of sample realizations: ",length(object$y),"\n",sep="")
  cat("Number of evaluation points: ",length(object$y.eval),"\n",sep="")
  cat("Bandwidths: h.y = ",object$h[1],", h.x = ",object$h[2],"\n",sep="")
  cat("Degree of local polynomial: ",object$degree,"\n",sep="")
  if(!is.null(object$f.yx.integral)) cat("Integral of estimate (prior to adjustment): ",object$f.yx.integral,"\n",sep="")
  cat("Number of cores used in parallel processing for kernel sum: ",object$ksum.cores,"\n",sep="")
  cat("Number of cores used in parallel processing for degree selection: ",object$degree.cores,"\n",sep="")
  cat("Number of cores used in parallel processing for multistart optimization: ",object$nmulti.cores,"\n",sep="")
  cat("Total number of cores used in parallel processing: ",object$ksum.cores*object$degree.cores*object$nmulti.cores,"\n",sep="")
  cat("Elapsed time (total): ",formatC(object$secs.elapsed,format="f",digits=2)," seconds\n",sep="")
  cat("Optimization and estimation time: ",formatC(object$secs.estimate+sum(object$secs.optim.mat),format="f",digits=2)," seconds\n",sep="")
  cat("Optimization and estimation time per core: ",formatC((object$secs.estimate+sum(object$secs.optim.mat))/(object$ksum.cores*object$degree.cores*object$nmulti.cores),format="f",digits=2)," seconds/core\n",sep="")
  cat("Parallel efficiency: ",formatC(object$secs.elapsed/(object$secs.estimate+sum(object$secs.optim.mat)),format="f",digits=2),
      " (allow for overhead and blocking, ideal = ",formatC(1/(object$ksum.cores*object$degree.cores*object$nmulti.cores),format="f",digits=2),")\n",sep="")
  cat("\n")
  invisible()
}
