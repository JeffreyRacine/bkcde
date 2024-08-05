## Description: This file contains the R functions used in the simulation study
## for implementing the boundary kernel conditional density estimator with
## cross-validated bandwidth and polynomial order selection.  The functions are
## used in the simulation study in the file mc.R, licensed under
## GPL-3.0-or-later; written by Jeffrey Racine <racinej@mcmaster.ca>

## mclapply() and mcmapply() from the parallel package are used throughout
## instead of lapply() and mapply() to allow for multi-core computing, and when
## ksum.cores, optim.degree.cores, or optim.nmulti.cores are set to a value greater than 1,
## the number of cores used in the parallel processing is set to the value of
## the respective argument.  But when these are set to 1, the number of cores
## used in the parallel processing is set to 1, i.e., serial processing occurs
## exactly as if lapply() and mapply() were being used. If verbose=TRUE is
## enabled, it appears warnings are most likely to appear immediately running
## things in serial mode.

## The functions are briefly described below.  The functions include
## integrate.trapezoidal(), NZD(), EssDee(), kernel.bk(), log.likelihood(),
## bkcde.loo(), bkcde(), bkcde.default(), bkcde.optim(), plot.bkcde(), and
## SCSrank() (the last from the MCSPAN package [Multiple contrast tests and
## simultaneous confidence intervals]).

## integrate.trapezoidal() computes the cumulative integral at each sample
## realization using the trapezoidal rule and the cumsum() function as we need
## to compute this in a computationally efficient manner to render some
## estimators "proper" (non-negative and integrating to 1). Vectors are paired
## naturally but need not be ordered (the function will take care of that for
## you).  The function also includes a correction term to ensure the integral is
## correct at the boundary.

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

## NZD() is the "No Zero Divide" (NZD) function (so e.g., 0/0 = 0) based on
## accepted coding practice for a variety of languages

NZD <- function(a) {
  ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
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

## kernel.bk() is the doubly truncated Gaussian boundary kernel function from
## Racine et al 2024

kernel.bk <- function(x,X,h,a=-Inf,b=Inf) {
  ## Checking for bounds involves a bit of overhead (20%), so here we presume a
  ## check is performed outside of this function - make sure this is the case!
  ## ifelse(X < a | X > b, 0, dnorm((x-X)/h)/(h*(pnorm((b-x)/h)-pnorm((a-x)/h))))
  dnorm((x-X)/h)/(h*(pnorm((b-x)/h)-pnorm((a-x)/h)))
}

## log.likelihood() returns a likelihood function that supports constant,
## smooth, and trim approaches for dealing with density estimates (delete-one)
## that may be improper and, in particular, negative. Note we use the smallest
## non-zero normalized floating-point number (a power of the radix, i.e.,
## double.base ^ double.min.exp, normally 2.225074e-308) as the cutoff value for
## the penalty. There are some slightly smaller numbers possible, i.e.
## log(4.450148e-324) = -744.4401, log(4.450148e-325) = -Inf, while
## log(.Machine$double.xmin) = -708.3964 (testing on Apple Silicon M2 Max), but
## this gives a very short range around zero where penalties will be constant
## (but appropriate) and otherwise provides a smooth function that, to my way of
## thinking, seems to sensibly handles negative delete-one values. Since
## log(0)=-Inf and since we need to return finite values to any numerical
## optimization function, we need to trap this case hence checking whether the
## vector of delete-one values exceeds 0 or not (we check that it is greater
## than .Machine$double.xmin).  You can create a vector containing negative
## values then plot it to see how the penalty functions behave for different
## choices (they are identical for positive values). If you don't approve of
## this hack simply use the two-lines in the "constant" case.

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
                                                                                      round(h[1],5),
                                                                                      ", h.x = ",
                                                                                      round(h[2],5),
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
                                                                                                round(h[1],5),
                                                                                                ", h.x = ",
                                                                                                round(h[2],5),
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

## bkcde.loo() is the leave-one-out likelihood cross-validation function that
## supports local polynomial estimation of degree p (raw polynomials or
## orthogonal polynomials can be used and appear to provide identical results
## for modest p.max)

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
    f.loo <- as.numeric(mcmapply(function(i){kernel.bk.x<-kernel.bk(x[i],x[-i],h[2],x.lb,x.ub);mean(kernel.bk(y[i],y[-i],h[1],y.lb,y.ub)*kernel.bk.x)/NZD(mean(kernel.bk.x))},1:length(y),mc.cores=ksum.cores))
  } else {
    X.poly <- poly(x,raw=poly.raw,degree=degree)
    X <- cbind(1,X.poly)
    f.loo <- as.numeric(mcmapply(function(i){beta.hat<-coef(lm.wfit(x=X[-i,,drop=FALSE],y=kernel.bk(y[i],y[-i],h[1],y.lb,y.ub),w=NZD(kernel.bk(x[i],x[-i],h[2],x.lb,x.ub))));beta.hat[!is.na(beta.hat)]%*%t(X[i,!is.na(beta.hat), drop = FALSE])},1:length(y),mc.cores=ksum.cores))
  }
  return(sum(log.likelihood(f.loo,penalty.method=penalty.method,penalty.cutoff=penalty.cutoff,verbose=verbose,degree=degree,h=h)))
}

## bckde() and bkcde.default() compute the conditional density \hat f(y|x)
## where, if no bandwidth is provided, then likelihood cross-validation is used
## to select the bandwidths and polynomial order via numerical optimization with
## 5 restarts by default and polynomial orders 0,1,...,5 by default. Restarting
## is used in an attempt to maximize the likelihood and (hopefully) avoid local
## optima. This function supports local polynomial orders [0,1,...,n-1] where n
## is the number of sample realizations (raw polynomials or orthogonal
## polynomials can be used and appear to provide identical results for modest
## degree.max)

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
                          degree.max=3,
                          degree.min=0,
                          degree=NULL,
                          fitted.cores=12,
                          ksum.cores=1,
                          n.grid=10,
                          n.integrate=1000,
                          nmulti=3,
                          optim.degree.cores=NULL,
                          optim.nmulti.cores=NULL,
                          penalty.cutoff=.Machine$double.xmin,
                          penalty.method=c("smooth","constant","trim"),
                          poly.raw=FALSE,
                          proper.cores=12,
                          proper=TRUE,
                          verbose=FALSE,
                          ...) {
  ## Perform some argument checking. In this function parallel processing takes
  ## place over the number of multistarts, so ideally the number of cores
  ## requested would be equal to the number of multistarts (this is particularly
  ## useful to avoid local optima in the optimization of the bandwidths)
  if(is.null(x)) stop("must provide x in bkcde()")
  if(is.null(y)) stop("must provide y in bkcde()")
  if(length(x) != length(y)) stop("length of x must be equal to length of y in bkcde()")
  if(!is.null(x.eval) & is.null(y.eval) & length(x.eval) != length(y)) stop("length of x.eval must be equal to length of y in bkcde() when y.eval is NULL")
  if(!is.null(x.eval) & !is.null(y.eval) & length(x.eval) != length(y.eval)) stop("length of x.eval must be equal to length of y.eval in bkcde() when x.eval and y.eval are not NULL")
  if(!is.null(y.eval) & is.null(x.eval)) stop("must provide x.eval in bkcde() when y.eval is not NULL")
  ## We set x.eval and y.eval to short sequences if they are not provided to
  ## avoid excessive computation with large samples when only x and y are
  ## provided (essentially what would be the default for plot.3D.n.grid in
  ## plot.bkcde()). Of course these can be changed, and plot will override them
  ## if desired, as will predict.bkcde().
  if(is.null(x.eval) & is.null(y.eval)) {
    data.grid <- expand.grid(seq(min(x),max(x),length=n.grid),seq(min(y),max(y),length=n.grid))
    x.eval <- data.grid$Var1
    y.eval <- data.grid$Var2
  }
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
  if(!is.null(h) & is.null(degree)) stop("must provide degree in bkcde() when h is not NULL")
  if(!is.null(h) & length(h) != 2) stop("h must be a vector of length 2 in bkcde()")
  if(!is.null(degree) && (degree < 0 | degree >= length(y))) stop("degree must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.min < 0 | degree.min >= length(y)) stop("degree.min must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.max < 0 | degree.max >= length(y)) stop("degree.max must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.min > degree.max) stop("degree.min must be <= degree.max in bkcde()")
  if(ksum.cores < 1) stop("ksum.cores must be at least 1 in bkcde()")
  if(proper.cores < 1) stop("proper.cores must be at least 1 in bkcde()")
  if(fitted.cores < 1) stop("fitted.cores must be at least 1 in bkcde()")
  if(!is.null(optim.degree.cores) && optim.degree.cores < 1) stop("optim.degree.cores must be at least 1 in bkcde()")
  if(!is.null(optim.nmulti.cores) && optim.nmulti.cores < 1) stop("optim.nmulti.cores must be at least 1 in bkcde()")
  if(is.null(optim.degree.cores)) optim.degree.cores <- degree.max-degree.min+1
  if(is.null(optim.nmulti.cores)) optim.nmulti.cores <- nmulti
  penalty.method <- match.arg(penalty.method)
  if(penalty.cutoff <= 0) stop("penalty.cutoff must be positive in bkcde()")
  secs.start.total <- Sys.time()
  ## If no bandwidth is provided, then likelihood cross-validation is used to
  ## obtain the bandwidths and polynomial order (use ksum.cores,
  ## optim.degree.cores, optim.nmulti.cores)
  if(is.null(h)) {
    optim.out <- bkcde.optim(x=x,
                             y=y,
                             y.lb=y.lb,
                             y.ub=y.ub,
                             x.lb=x.lb,
                             x.ub=x.ub,
                             degree.max=degree.max,
                             degree.min=degree.min,
                             ksum.cores=ksum.cores,
                             nmulti=nmulti,
                             optim.degree.cores=optim.degree.cores,
                             optim.nmulti.cores=optim.nmulti.cores,
                             penalty.cutoff=penalty.cutoff,
                             penalty.method=penalty.method,
                             poly.raw=poly.raw,
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
  ## Compute the fitted conditional density estimate (use fitted.cores)
  if(degree == 0) {
    ## For degree 0 don't invoke the overhead associated with lm.wfit(), just
    ## compute the estimate \hat f(y|x) as efficiently as possible
    f.yx <- as.numeric(mcmapply(function(i){kernel.bk.x<-kernel.bk(x.eval[i],x,h[2],x.lb,x.ub);mean(kernel.bk(y.eval[i],y,h[1],y.lb,y.ub)*kernel.bk.x)/NZD(mean(kernel.bk.x))},1:length(y.eval),mc.cores=fitted.cores))
  } else {
    ## Choice of raw or orthogonal polynomials
    X.poly <- poly(x,raw=poly.raw,degree=degree)
    X <- cbind(1,X.poly)
    ## For degree > 0 we use, e.g., lm(y~I(x^2)) and fitted values from the
    ## regression to estimate \hat f(y|x) rather than the intercept term from
    ## lm(y-I(x[i]-X)^2), which produce identical results for raw polynomials
    f.yx <- as.numeric(mcmapply(function(i){beta.hat<-coef(lm.wfit(x=X,y=kernel.bk(y.eval[i],y,h[1],y.lb,y.ub),w=NZD(kernel.bk(x.eval[i],x,h[2],x.lb,x.ub))));beta.hat[!is.na(beta.hat)]%*%t(cbind(1,predict(X.poly,x.eval[i]))[,!is.na(beta.hat),drop = FALSE])},1:length(y.eval),mc.cores=fitted.cores))
  }
  ## Ensure the estimate is proper (use proper.cores over unique(x.eval) which
  ## could be < # proper.cores allocated)
  if(proper) {
    ## Create a sequence of values along an appropriate grid to compute the integral.
    if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
    if(is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(y.lb,extendrange(y,f=10)[2],length=n.integrate)
    if(!is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(extendrange(y,f=10)[1],y.ub,length=n.integrate)
    if(!is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(extendrange(y,f=10)[1],extendrange(y,f=10)[2],length=n.integrate)
    ## Note that we have a conditional density f(y|x) with potentially repeated
    ## x values, so for each unique x in x.eval we ensure f(y|x) is proper
    ## (avoid unnecessary computation, particularly when x.eval contains a
    ## constant or x.eval is taken from expand.grid() and contains a repeated
    ## sequence of identical values). Test for unique values of x.eval to reduce
    ## potential computation. We use mclapply to return the list of integrals
    ## evaluated on y.seq for all unique x.eval values (this is done in parallel
    ## and can save substantial time)
    x.eval.unique <- unique(x.eval)
    int.f.seq.pre.neg <- numeric()
    int.f.seq <- numeric()
    int.f.seq.post <- numeric()
    ## We test for only 1 unique value of x.eval to avoid parallel processing in
    ## the outer mcmapply call and invoke fitting the mcmapply sequence of
    ## f(y|x) values with proper.cores
    proper.out <- mclapply(1:length(x.eval.unique),function(j) {
      K <- kernel.bk(x.eval.unique[j],x,h[2],x.lb,x.ub)
      if(degree == 0) {
        f.seq <- as.numeric(mcmapply(function(i){mean(kernel.bk(y.seq[i],y,h[1],y.lb,y.ub)*K)/NZD(mean(K))},1:n.integrate,mc.cores=ifelse(length(x.eval.unique)>1,1,proper.cores)))
      } else {
        X.poly <- poly(x,raw=poly.raw,degree=degree)
        X <- cbind(1,X.poly)
        X.eval <- cbind(1,predict(X.poly,x.eval.unique[j]))
        f.seq <- as.numeric(mcmapply(function(i){beta.hat<-coef(lm.wfit(x=X,y=kernel.bk(y.seq[i],y,h[1],y.lb,y.ub),w=NZD(K)));beta.hat[!is.na(beta.hat)]%*%t(X.eval[,!is.na(beta.hat),drop = FALSE])},1:n.integrate,mc.cores=ifelse(length(x.eval.unique)>1,1,proper.cores)))
      }
      ## Compute integral of f.seq including any possible negative values
      int.f.seq.pre.neg[j]<- integrate.trapezoidal(y.seq,f.seq)[length(y.seq)]
      ## Set any possible negative f.seq values to 0
      f.seq[f.seq < 0] <- 0
      ## Compute integral of f.seq after setting any possible negative values to 0
      int.f.seq[j] <- integrate.trapezoidal(y.seq,f.seq)[length(y.seq)]
      ## Compute integral of f.seq after setting any possible negative values to 0
      ## and correcting to ensure final estimate integrates to 1
      int.f.seq.post[j] <- integrate.trapezoidal(y.seq,f.seq/int.f.seq[j])[length(y.seq)]
      return(list(int.f.seq.pre.neg=int.f.seq.pre.neg[j],
                  int.f.seq=int.f.seq[j],
                  int.f.seq.post=int.f.seq.post[j]))
    },mc.cores = ifelse(length(x.eval.unique)>1,proper.cores,1))
    ## Now gather the results, correct for negative entries then divide elements
    ## of f.xy by the corresponding integral (one for each x.eval.unique) to
    ## ensure the estimate is proper
    if(verbose & any(f.yx < 0)) warning("negative density estimate reset to 0 via option proper=TRUE in bkcde() [degree = ",
                                        degree,
                                        ", j = ",
                                        length(isTRUE(f.yx < 0)),
                                        " element(s), h.y = ",
                                        round(h[1],5),
                                        ", h.x = ",
                                        round(h[2],5),
                                        "]",
                                        immediate. = TRUE)
    f.yx[f.yx < 0] <- 0
    int.f.seq.pre.neg <- sapply(proper.out, function(x) x$int.f.seq.pre.neg)
    int.f.seq <- sapply(proper.out, function(x) x$int.f.seq)
    int.f.seq.post <- sapply(proper.out, function(x) x$int.f.seq.post)
    f.yx <- f.yx/int.f.seq[match(x.eval, x.eval.unique)]
    ## As a summary measure report the mean of the integrals
    int.f.seq.pre.neg <- mean(int.f.seq.pre.neg)
    int.f.seq <- mean(int.f.seq)
    int.f.seq.post <- mean(int.f.seq.post)
  } else {
    int.f.seq.pre.neg <- NA
    int.f.seq <- NA
    int.f.seq.post <- NA
    if(verbose & any(f.yx < 0)) warning("negative density estimate encountered, consider option proper=TRUE in bkcde() [degree = ",
                                        degree,
                                        ", ", 
                                        length(isTRUE(f.yx < 0)),
                                        " element(s), h.y = ",
                                        round(h[1],5),
                                        ", h.x = ",
                                        round(h[2],5),
                                        "]",
                                        immediate. = TRUE)
  }
  return.list <- list(convergence.mat=convergence.mat,
                      convergence.vec=convergence.vec,
                      convergence=convergence,
                      degree.mat=degree.mat,
                      degree.max=degree.max,
                      degree.min=degree.min,
                      degree=degree,
                      fitted.cores=fitted.cores,
                      f.yx.integral.post=int.f.seq.post,
                      f.yx.integral.pre.neg=int.f.seq.pre.neg,
                      f.yx.integral=int.f.seq,
                      f=f.yx,
                      h.mat=h.mat,
                      h=h,
                      ksum.cores=ksum.cores,
                      optim.degree.cores=optim.degree.cores,
                      optim.nmulti.cores=optim.nmulti.cores,
                      proper.cores=proper.cores,
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

## bkcde.optim() conducts numerical optimization for bandwidth selection and
## polynomial order in bkcde() using the optim() function with the L-BFGS-B
## method which allows box constraints, that is each variable can be given a
## lower and/or upper bound (bandwidths must be positive so this is necessary).

bkcde.optim <- function(x=x,
                        y=y,
                        y.lb=y.lb,
                        y.ub=y.ub,
                        x.lb=x.lb,
                        x.ub=x.ub,
                        degree.max=degree.max,
                        degree.min=degree.min,
                        ksum.cores=ksum.cores,
                        nmulti=nmulti,
                        optim.degree.cores=optim.degree.cores,
                        optim.nmulti.cores=optim.nmulti.cores,
                        penalty.cutoff=penalty.cutoff,
                        penalty.method=penalty.method,
                        poly.raw=poly.raw,
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
  ## Here we conduct optimization over all models (i.e., polynomial orders) in
  ## parallel each having degree p in [degree.min,degree.max]
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
    },mc.cores = optim.nmulti.cores)
    optim.out <- nmulti.return[[which.max(sapply(nmulti.return, function(x) x$value))]]
    optim.out$value.vec <- sapply(nmulti.return, function(x) x$value)
    optim.out$degree.vec <- sapply(nmulti.return, function(x) x$degree)
    optim.out$convergence.vec <- sapply(nmulti.return, function(x) x$convergence)
    optim.out$secs.optim.vec <- sapply(nmulti.return, function(x) x$secs.optim)
    optim.out
  },mc.cores = optim.degree.cores)
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

## persp() does not allow ylim and zlim to be set to NULL, so to be passable
## this function avoids unnecessary duplication of code in plot.bkcde()

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

## plot.bkcde() is used to plot the results of the boundary kernel CDE along
## with bootstrap confidence intervals generated as either pointwise,
## Bonferroni, or simultaneous intervals. A simple bootstrap mean correction is
## applied (there are likely better ways of doing this outside of a bootstrap
## procedure). A handful of options are available, including returning the
## confidence intervals (pointwise, Bonferroni and simultaneous) and estimates.
## Note that a large number of bootstrap replications should be used,
## particularly if Bonferroni corrections are requested, and in such cases the
## number of bootstrap replications should increase with the grid size (i.e., the
## number of evaluation points). The function uses the SCSrank() function from
## the MCPAN package to compute simultaneous confidence intervals.

plot.bkcde <- function(x,
                       B = 3999,
                       alpha = 0.05,
                       ci = FALSE,
                       ci.bias.correct = TRUE,
                       ci.cores = NULL,
                       ci.method = c("Pointwise","Bonferroni","Simultaneous","all"),
                       ci.preplot = TRUE,
                       ci.progress = TRUE,
                       fitted.cores = NULL,
                       ksum.cores = NULL,
                       phi = NULL,
                       plot.2D.n.grid = 100,
                       plot.2D.y.grid = NULL,
                       plot.3D = TRUE,
                       plot.3D.n.grid = 10,
                       plot.3D.x.grid = NULL,
                       plot.3D.y.grid = NULL,
                       plot.behavior = c("plot","plot-data","data"),
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
                       ...) {
  if(!inherits(x,"bkcde")) stop("x must be of class bkcde in plot.bkcde()")
  if(!is.logical(ci)) stop("ci must be logical in plot.bkcde()")
  if(!is.logical(ci.progress)) stop("ci.progress must be logical in plot.bkcde()")
  if(!is.logical(ci.preplot)) stop("ci.preplot must be logical in plot.bkcde()")
  ## Note that the Bonferroni method places a restriction on the smallest number
  ## of bootstrap replications such that (alpha/(2*plot.2D.n.grid))*(B+1) or
  ## (alpha/(2*plot.3D.n.grid^2))*(B+1) is a positive integer - this is on the
  ## user to ensure. The simultaneous method will almost certainly require fewer
  ## bootstrap replications than the Bonferroni method. The pointwise method
  ## requires the fewest, alpha*(B+1) must be a positive integer. So, for the
  ## Bonferroni bounds the smallest number of bootstrap replications using the
  ## default (alpha=0.05, plot.2D.n.grid=100 or alpha=0.05, plot.3D.n.grid=10 is
  ## 3999. For non-default values the user should adjust B accordingly.
  ci.method <- match.arg(ci.method)
  plot.behavior <- match.arg(plot.behavior)
  if(alpha <= 0 | alpha >= 1) stop("alpha must lie in (0,1) in plot.bkcde()")
  if(B < 1) stop("B must be at least 1 in plot.bkcde()")
  if(!is.null(ci.cores) && ci.cores < 1) stop("ci.cores must be at least 1 in plot.bkcde()")
  ci.pw.lb <- ci.pw.ub <- ci.bf.lb <- ci.bf.ub <- ci.sim.lb <- ci.sim.ub <- bias.vec <- NULL
  if(!is.null(proper) & !is.logical(proper)) stop("proper must be logical in plot.bkcde()")
  if(plot.3D & !is.null(x.eval)) warning("x.eval passed but ignored in plot.bkcde() when plot.3D = TRUE",immediate. = TRUE)
  if(!is.null(plot.3D.x.grid) & !is.null(plot.3D.y.grid) & length(plot.3D.x.grid) != length(plot.3D.y.grid)) stop("length of plot.3D.x.grid must be equal to length of plot.3D.y.grid in plot.bkcde()")
  if(plot.2D.n.grid < 2) stop("plot.2D.n.grid must be at least 2 in plot.bkcde()")
  if(plot.3D.n.grid < 2) stop("plot.3D.n.grid must be at least 2 in plot.bkcde()")
  if(ci.preplot==FALSE & ci==FALSE & ci.method != "data") stop("ci.preplot must be TRUE when ci is TRUE and ci.method is not 'data' in plot.bkcde()")
  if(is.null(proper)) proper <- x$proper
  if(!plot.3D & is.null(x.eval)) x.eval <- median(x$x.eval)
  if(is.null(ksum.cores)) ksum.cores <- x$ksum.cores
  ## proper.cores and fitted.cores are used in the predict() function and only
  ## in the bootstrap if ci=TRUE and ci.cores>1
  if(is.null(proper.cores)) proper.cores <- detectCores()
  if(is.null(fitted.cores)) fitted.cores <- detectCores()
  secs.start <- Sys.time()
  ## For the user, whether ci=TRUE or not get the estimate plotted asap
  ## otherwise they are faced with a blank screen
  if(plot.3D) {
    ## Plot 3D
    if(is.null(plot.3D.x.grid)) {
      x.grid <- seq(min(x$x.eval),max(x$x.eval),length=plot.3D.n.grid)
      y.grid <- seq(min(x$y.eval),max(x$y.eval),length=plot.3D.n.grid)
    } else {
      x.grid <- plot.3D.x.grid
      y.grid <- plot.3D.y.grid
      plot.3D.n.grid <- length(x.grid)
    }
    if(length(unique(x.grid))==1) stop("only one unique x.eval value, cannot deploy persp() in plot.bkcde() (perhaps call bkcde() with non-unique x.eval OR provide plot.3D.x.grid?)")
    if(length(unique(y.grid))==1) stop("only one unique y.eval value, cannot deploy persp() in plot.bkcde() (perhaps call bkcde() with non-unique y.eval OR plot.3D.y.grid?)")
    data.grid <- expand.grid(x.grid,y.grid)
    x.plot.eval <- data.grid$Var1
    y.plot.eval <- data.grid$Var2
    x.fitted <- predict(x,newdata=data.frame(x=x.plot.eval,y=y.plot.eval),proper=proper,ksum.cores=ksum.cores,fitted.cores=fitted.cores,proper.cores=proper.cores,...)
    predict.mat <- matrix(x.fitted,plot.3D.n.grid,plot.3D.n.grid)
    if(is.null(theta)) theta <- 120
    if(is.null(phi)) phi <- 45
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- "y"
    if(is.null(zlab)) zlab <- "f(y|x)"
    if(ci.preplot & plot.behavior != "data") {
      persp.lim(x=x.grid,y=y.grid,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,...)    
    }
  } else if(!plot.3D) {
    ## Plot 2D
    predict.mat <- NULL
    if(is.null(plot.2D.y.grid)) {
      y.plot.eval <- y.grid <- seq(min(x$y.eval),max(x$y.eval),length=plot.2D.n.grid)
    } else {
      y.plot.eval <- y.grid <- plot.2D.y.grid
      plot.2D.n.grid <- length(y.grid)
    }
    x.plot.eval <- x.grid <- rep(x.eval,length(y.plot.eval))
    x.fitted <- predict(x,newdata=data.frame(x=x.plot.eval,y=y.plot.eval),proper=proper,ksum.cores=ksum.cores,...)
    if(is.null(sub)) sub <- paste("(degree = ",x$degree,", h.y = ",round(x$h[1],3), ", h.x = ",round(x$h[2],3),", n = ",length(x$y),")",sep="")
    if(is.null(ylab)) ylab <- "f(y|x)"
    if(is.null(xlab)) xlab <- paste("y|x=",round(x.eval,digits=2),sep="")
    if(is.null(type)) type <- "l"
    if(ci.preplot & plot.behavior != "data") {
      plot(y.plot.eval[order(y.plot.eval)],x.fitted[order(y.plot.eval)],
           sub=sub,
           ylim=ylim,
           ylab=ylab,
           xlab=xlab,
           type=type,
           panel.first=grid(lty=1),
           ...)
    }
  }
  if(ci) {
    ## All processing goes into computing the matrix of bootstrap estimates, so
    ## once this is done it makes sense to then generate all three types of
    ## confidence intervals
    if(is.null(ci.cores)) ci.cores <- detectCores()
    mcmapply.progress <- function(...,progress=TRUE) {
      if(progress) {
        pbmcmapply(...)
      } else {
        mcmapply(...)
      }
    }
    boot.mat <- t(mcmapply.progress(function(b){
      ii <- sample(1:length(x$y),replace=TRUE)
      bkcde(h=x$h,
            x=x$x[ii],
            y=x$y[ii],
            x.eval=x.plot.eval,
            y.eval=y.plot.eval,
            y.lb=x$y.lb,
            y.ub=x$y.ub,
            x.lb=x$x.lb,
            x.ub=x$x.ub,
            fitted.cores=ifelse(ci.cores>1,1,fitted.cores),
            proper=proper,
            proper.cores=ifelse(ci.cores>1,1,proper.cores),
            ksum.cores=ksum.cores,
            degree=x$degree)$f
    },1:B,mc.cores=ci.cores,progress=ci.progress))
    if(ci.bias.correct) {
      bias.vec <- colMeans(boot.mat) - x.fitted
      boot.mat <- sweep(boot.mat,2,bias.vec,"-")
      if(proper) boot.mat <- pmax(boot.mat,0)
    }
    ci.pw.lb <- apply(boot.mat, 2, quantile, probs = alpha / 2)
    ci.pw.ub <- apply(boot.mat, 2, quantile, probs = 1 - alpha / 2)
    ci.bf.lb <- apply(boot.mat, 2, quantile, probs = alpha / (2 * length(y.plot.eval)))
    ci.bf.ub <- apply(boot.mat, 2, quantile, probs = 1 - alpha / (2 * length(y.plot.eval)))
    ci.SCS <- SCSrank(boot.mat, conf.level=1-alpha)$conf.int
    ci.sim.lb <- ci.SCS[,1]
    ci.sim.ub <- ci.SCS[,2]
    cat("\r                                                                                             ")
  } else {
    if(is.null(ylim)) ylim <- NULL
  }
  if(plot.behavior != "data") {
    if(plot.3D & ci) {
      ## Plot 3D again with zlim set for the confidence intervals (not checking for if(is.null(zlim)) yet)
      if(ci.method == "Pointwise") {
        if(is.null(zlim)) zlim <-  range(c(x.fitted,ci.pw.lb,ci.pw.ub))
      } else if(ci.method == "Bonferroni") {
        if(is.null(zlim)) zlim <-  range(c(x.fitted,ci.bf.lb,ci.bf.ub))
      } else if(ci.method == "Simultaneous") {
        if(is.null(zlim)) zlim <-  range(c(x.fitted,ci.sim.lb,ci.sim.ub))
      } else {
        if(is.null(zlim)) zlim <-  range(c(x.fitted,ci.pw.lb,ci.pw.ub,ci.bf.lb,ci.bf.ub,ci.sim.lb,ci.sim.ub))
      }
      ## Unlike plot() persp() does accept a null ylim argument so we need to check...
      if(ci.method == "Pointwise") {
        ## First lower, then plot, then upper (surfaces)
        persp.lim(x=x.grid,y=y.grid,z=matrix(ci.pw.lb,plot.3D.n.grid,plot.3D.n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.grid,y=y.grid,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.grid,y=y.grid,z=matrix(ci.pw.ub,plot.3D.n.grid,plot.3D.n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci.method == "Bonferroni") {
        persp.lim(x=x.grid,y=y.grid,z=matrix(ci.bf.lb,plot.3D.n.grid,plot.3D.n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.grid,y=y.grid,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.grid,y=y.grid,z=matrix(ci.bf.ub,plot.3D.n.grid,plot.3D.n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci.method == "Simultaneous") {
        persp.lim(x=x.grid,y=y.grid,z=matrix(ci.sim.lb,plot.3D.n.grid,plot.3D.n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.grid,y=y.grid,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.grid,y=y.grid,z=matrix(ci.sim.ub,plot.3D.n.grid,plot.3D.n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci & ci.method == "all") {
        warning("plotting all confidence intervals in plot.bkcde() may be visually overwhelming, ignored (but you can retrieve the ci data via plot=FALSE)",immediate. = TRUE)
      }
    } else if(ci) {
      ## Plot 2D again with ylim set for the confidence intervals
      if(ci.method == "Pointwise") {
        if(is.null(ylim)) ylim <-  range(c(x.fitted,ci.pw.lb,ci.pw.ub))
      } else if(ci.method == "Bonferroni") {
        if(is.null(ylim)) ylim <-  range(c(x.fitted,ci.bf.lb,ci.bf.ub))
      } else if(ci.method == "Simultaneous") {
        if(is.null(ylim)) ylim <-  range(c(x.fitted,ci.sim.lb,ci.sim.ub))
      } else {
        if(is.null(ylim)) ylim <-  range(c(x.fitted,ci.pw.lb,ci.pw.ub,ci.bf.lb,ci.bf.ub,ci.sim.lb,ci.sim.ub))
      }
      ## First plot, then lower and upper (lines)
      plot(y.plot.eval[order(y.plot.eval)],x.fitted[order(y.plot.eval)],
           sub=sub,
           ylim=ylim,
           ylab=ylab,
           xlab=xlab,
           type=type,
           panel.first=grid(lty=1),
           ...)
      if(ci.method == "Pointwise") {
        lines(y.plot.eval[order(y.plot.eval)],ci.pw.lb[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],ci.pw.ub[order(y.plot.eval)],lty=2)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci.method == "Bonferroni") {
        lines(y.plot.eval[order(y.plot.eval)],ci.bf.lb[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],ci.bf.ub[order(y.plot.eval)],lty=2)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci & ci.method == "Simultaneous") {
        lines(y.plot.eval[order(y.plot.eval)],ci.sim.lb[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],ci.sim.ub[order(y.plot.eval)],lty=2)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci.method == "all") {
        lines(y.plot.eval[order(y.plot.eval)],ci.pw.lb[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],ci.pw.ub[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],ci.sim.lb[order(y.plot.eval)],lty=3)
        lines(y.plot.eval[order(y.plot.eval)],ci.sim.ub[order(y.plot.eval)],lty=3)
        lines(y.plot.eval[order(y.plot.eval)],ci.bf.lb[order(y.plot.eval)],lty=4)
        lines(y.plot.eval[order(y.plot.eval)],ci.bf.ub[order(y.plot.eval)],lty=4)
        legend("topright",
               legend=c("Estimated f(y|x)",
                        paste(100*(1-alpha),"% Pointwise CIs",sep=""),
                        paste(100*(1-alpha),"% Simultaneous CIs",sep=""),
                        paste(100*(1-alpha),"% Bonferroni CIs",sep="")
               ),lty=1:4,bty="n")
      }
    }
  } 
  if(plot.behavior != "plot") {
    ## Return in vector form for the user to plot as they wish (if
    ## plot.3D=TRUE convert to matrix using x.grid & y.grid)
    return(list(bias.vec=bias.vec,
                ci.bf.lb=ci.bf.lb,
                ci.bf.ub=ci.bf.ub,
                ci.cores=ci.cores,
                ci.pw.lb=ci.pw.lb,
                ci.pw.ub=ci.pw.ub,
                ci.sim.lb=ci.sim.lb,
                ci.sim.ub=ci.sim.ub,
                f=x.fitted,
                proper.cores=proper.cores,
                secs.elapsed=as.numeric(difftime(Sys.time(),secs.start,units="secs")),
                x.grid=x.grid,
                x=x.plot.eval,
                y.grid=y.grid,
                ylim=ylim,
                zlim=zlim))
  }
}

## fitted.bkcde() returns the estimated conditional density f(y|x) at the
## specified evaluation points used to estimate the density in the function
## call.

fitted.bkcde <- function(object, ...) {
  if(!inherits(object,"bkcde")) stop("object must be of class bkcde in fitted.bkcde()")
  return(object$f)
}

## predict.bkcde() returns the estimated conditional density f(y|x) at new evaluation points

predict.bkcde <- function(object, newdata, proper = NULL, ...) {
  if(!inherits(object,"bkcde")) stop("object must be of class bkcde in predict.bkcde()")
  if(!is.data.frame(newdata)) stop("newdata must be a data frame in predict.bkcde()")
  if(!all(names(newdata) %in% c("x","y"))) stop("newdata must contain columns x and y in predict.bkcde()")
  if(!is.null(proper) & !is.logical(proper)) stop("proper must be logical in predict.bkcde()")
  if(is.null(proper)) proper <- object$proper
  return(bkcde(h=object$h,
               x=object$x,
               y=object$y,
               x.eval=newdata$x,
               y.eval=newdata$y,
               y.lb=object$y.lb,
               y.ub=object$y.ub,
               x.lb=object$x.lb,
               x.ub=object$x.ub,
               proper=proper,
               degree=object$degree,
               ...)$f)
}

## summary.bkcde() provides a summary of the boundary kernel CDE object

summary.bkcde <- function(object, ...) {
  if(!inherits(object,"bkcde")) stop("object must be of class bkcde in summary.bkcde()")
  cat("Call:\n")
  cat("bkcde(h.y=",round(object$h[1],4),", h.x=",round(object$h[2],4),", x, y, x.eval, y.eval, y.lb=",object$y.lb,", y.ub=",object$y.ub,", x.lb=",object$x.lb,", x.ub=",object$x.ub,", degree=",object$degree,")\n",sep="")
  cat("\n")
  cat("Number of sample realizations: ",length(object$y),"\n",sep="")
  cat("Number of evaluation points: ",length(object$y.eval),"\n",sep="")
  cat("Bandwidths: h.y = ",object$h[1],", h.x = ",object$h[2],"\n",sep="")
  cat("Degree of local polynomial: ",object$degree,"\n",sep="")
  if(!is.na(object$f.yx.integral.pre.neg)) cat("Integral of estimate (pre any negativity correction): ",formatC(object$f.yx.integral.pre.neg,format="f",digits=12),"\n",sep="")
  if(!is.na(object$f.yx.integral)) cat("Integral of estimate (post negativity, prior to integration to 1 correction): ",formatC(object$f.yx.integral,format="f",digits=12),"\n",sep="")
  if(!is.na(object$f.yx.integral.post)) cat("Integral of estimate (post all corrections): ",formatC(object$f.yx.integral.post,format="f",digits=12),"\n",sep="")
  cat("Number of cores used for optimization in parallel processing for degree selection: ",object$optim.degree.cores,"\n",sep="")
  cat("Number of cores used for optimization in parallel processing for multistart optimization: ",object$optim.nmulti.cores,"\n",sep="")
  cat("Total number of cores used for optimization in parallel processing: ",object$ksum.cores*object$optim.degree.cores*object$optim.nmulti.cores,"\n",sep="")
  cat("Number of cores used in parallel processing for ensuring proper density: ",object$proper.cores,"\n",sep="")
  cat("Number of cores used in parallel processing for fitting density: ",object$fitted.cores,"\n",sep="")
  cat("Number of cores used in parallel processing for kernel sum: ",object$ksum.cores,"\n",sep="")
  cat("Elapsed time (total): ",formatC(object$secs.elapsed,format="f",digits=2)," seconds\n",sep="")
  cat("Optimization and estimation time: ",formatC(object$secs.estimate+sum(object$secs.optim.mat),format="f",digits=2)," seconds\n",sep="")
  cat("Optimization and estimation time per core: ",formatC((object$secs.estimate+sum(object$secs.optim.mat))/(object$ksum.cores*object$optim.degree.cores*object$optim.nmulti.cores),format="f",digits=2)," seconds/core\n",sep="")
  cat("Parallel efficiency: ",formatC(object$secs.elapsed/(object$secs.estimate+sum(object$secs.optim.mat)),format="f",digits=2),
      " (allow for overhead and blocking, ideal = ",formatC(1/(object$ksum.cores*object$optim.degree.cores*object$optim.nmulti.cores),format="f",digits=2),")\n",sep="")
  cat("\n")
  invisible()
}

## This function takes a subset of the (x,y) data, computes the optimal h and
## degree, then repeats 10 times and takes the median of the h vector and degree
## vector and returns those values. This is a fast way to compute the optimal
## bandwidth and polynomial degree based upon Racine, J.S. (1993), "An Efficient
## Cross-Validation Algorithm For Window Width Selection for Nonparametric
## Kernel Regression," Communications in Statistics, October, Volume 22, Issue
## 4, pages 1107-1114.

fast.optim <- function(x, y, n.sub = 1000, resamples = 10, proper=FALSE, ...) {
  if(!is.numeric(x)) stop("x must be numeric in fast.optim()")
  if(!is.numeric(y)) stop("y must be numeric in fast.optim()")
  if(length(x) != length(y)) stop("length of x must be equal to length of y in fast.optim()")
  if(!is.numeric(n.sub)) stop("n.sub must be numeric in fast.optim()")
  if(n.sub < 1) stop("n.sub must be at least 10 in fast.optim()")
  if(resamples < 1) stop("resamples must be at least 1 in fast.optim()")
  if(!is.logical(proper)) stop("proper must be logical in fast.optim()")
  n <- length(y)
  h.degree.mat <- matrix(NA,nrow=resamples,ncol=3)
  for(j in 1:resamples) {
    cat("\rResample ",j," of ",resamples,"...",sep="")
    ii <- sample(n,size=n.sub)
    bkcde.out <- bkcde(x=x[ii],y=y[ii],proper=proper,...)
    h.degree.mat[j,1] <- (bkcde.out$h[1]/EssDee(y[ii]))*n.sub^{1/6}
    h.degree.mat[j,2] <- (bkcde.out$h[2]/EssDee(x[ii]))*n.sub^{1/6}
    h.degree.mat[j,3] <- bkcde.out$degree
  }
  cat("\r                                          \r")
  ## Compute median of columns of h.degree.mat after rescaling for larger sample
  h.degree.mat[,1] <- h.degree.mat[,1]*EssDee(y)*n^{-1/6}
  h.degree.mat[,2] <- h.degree.mat[,2]*EssDee(x)*n^{-1/6}
  return(list(h=c(median(h.degree.mat[,1]),
                  median(h.degree.mat[,2])),
                  degree=round(median(h.degree.mat[,3])),
                  h.degree.mat=h.degree.mat))
}
