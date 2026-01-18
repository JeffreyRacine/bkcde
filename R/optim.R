## Optimization functions for bkcde package

## log.likelihood() returns a likelihood function that supports constant,
## smooth, and trim approaches for dealing with density estimates (delete-one)
## that may be improper and, in particular, negative.

log.likelihood <- function(delete.one.values,
                           cv.penalty.method=c("smooth","constant","trim","nonneg"),
                           cv.penalty.cutoff=.Machine$double.xmin,
                           verbose=FALSE,
                           degree=degree,
                           h=h) {
  cv.penalty.method <- match.arg(cv.penalty.method)
  if(cv.penalty.cutoff <= 0) stop("cv.penalty.cutoff must be positive in log.likelihood()")
  likelihood.vec <- numeric(length(delete.one.values))
  cutoff.val <- cv.penalty.cutoff
  log.cutoff <- log(cutoff.val)
  
  if(cv.penalty.method=="constant") {
    likelihood.vec[delete.one.values > cutoff.val] <- log(delete.one.values[delete.one.values > cutoff.val])
    likelihood.vec[delete.one.values <= cutoff.val] <- log.cutoff
    if(verbose && any(0 < delete.one.values & delete.one.values < cutoff.val)) warning("delete-one density lies in constant cutoff zone in log.likelihood() [degree = ",
                                                                                      degree,
                                                                                      ", ",
                                                                                      length(delete.one.values[0 < delete.one.values & delete.one.values < cutoff.val]),
                                                                                      " element(s), h.y = ",
                                                                                      round(h[1],5),
                                                                                      ", h.x = ",
                                                                                      round(h[2],5),
                                                                                      "]",
                                                                                      immediate. = TRUE)
  } else if(cv.penalty.method=="smooth") {
    likelihood.vec[delete.one.values > cutoff.val] <- log(delete.one.values[delete.one.values > cutoff.val])
    likelihood.vec[delete.one.values < -cutoff.val] <- -log(abs(delete.one.values[delete.one.values < -cutoff.val]))+2*log.cutoff
    likelihood.vec[-cutoff.val < delete.one.values & delete.one.values < cutoff.val] <- log.cutoff
    if(verbose && any(-cutoff.val < delete.one.values & delete.one.values < cutoff.val)) warning("delete-one density lies in smooth cutoff zone in log.likelihood() [degree = ",
                                                                                                degree,
                                                                                                ", ",
                                                                                                length(delete.one.values[-cutoff.val < delete.one.values & delete.one.values < cutoff.val]),
                                                                                                " element(s), h.y = ",
                                                                                                round(h[1],5),
                                                                                                ", h.x = ",
                                                                                                round(h[2],5),
                                                                                                "]",
                                                                                                immediate. = TRUE)
  } else if(cv.penalty.method=="trim") {
    likelihood.vec[delete.one.values > cutoff.val] <- log(delete.one.values[delete.one.values>cutoff.val])
    likelihood.vec[delete.one.values <= cutoff.val] <- NA
    likelihood.vec <- likelihood.vec[!is.na(likelihood.vec)]
  } else if(cv.penalty.method=="nonneg") {
    if(any(delete.one.values < cv.penalty.cutoff)) {
      likelihood.vec <- log(rep(cutoff.val,length(delete.one.values)))
    } else {
      likelihood.vec <- log(delete.one.values)
    }
  }
  return(likelihood.vec)
}

## bkcde.optim.fn() is the cross-validation objective function

bkcde.optim.fn <- function(h=NULL,
                           x=NULL,
                           y=NULL,
                           x.eval=NULL,
                           y.eval=NULL,
                           y.lb=NULL,
                           y.ub=NULL,
                           x.lb=NULL,
                           x.ub=NULL,
                           poly.raw=FALSE,
                           degree=NULL,
                           n.integrate=NULL,
                           optim.ksum.cores=1,
                           cv.penalty.method=NULL,
                           cv.penalty.cutoff=NULL,
                           verbose=FALSE,
                           bwmethod=NULL,
                           proper.cv=NULL,
                           X=NULL,
                           X.eval=NULL) {
  if(y.lb>=y.ub) stop("y.lb must be less than y.ub in bkcde.optim.fn()")
  if(x.lb>=x.ub) stop("x.lb must be less than x.ub in bkcde.optim.fn()")
  if(is.null(x)) stop("must provide x in bkcde.optim.fn()")
  if(is.null(y)) stop("must provide y in bkcde.optim.fn()")
  if(is.null(degree)) stop("must provide degree in bkcde.optim.fn()")
  if(!is.logical(poly.raw)) stop("poly.raw must be logical in bkcde.optim.fn()")
  if(optim.ksum.cores < 1) stop("optim.ksum.cores must be at least 1 in bkcde.optim.fn()")
  if(is.null(cv.penalty.method)) stop("must provide cv.penalty.method in bkcde.optim.fn()")
  if(is.null(cv.penalty.cutoff)) stop("must provide cv.penalty.cutoff in bkcde.optim.fn()")
  if(is.null(bwmethod)) stop("must provide bwmethod in bkcde.optim.fn()")
  if(is.null(n.integrate)) stop("must provide n.integrate in bkcde.optim.fn()")
  if(is.null(proper.cv)) stop("must provide proper in bkcde.optim.fn()")
  if(degree < 0 || degree >= length(y)) stop("degree must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde.optim.fn()")

  denom.x <- h[2]*(if(is.infinite(x.ub)) 1 else pnorm((x.ub-x)/h[2]) - (if(is.infinite(x.lb)) 0 else pnorm((x.lb-x)/h[2])))
  denom.y <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y)/h[1])))

  ## When non-negativity penalty is requested for degree > 0 models, 
  ## we must check if the density is negative at evaluation points.
  if(cv.penalty.method=="nonneg" && degree>0) {
    if(!identical(y,y.eval) || !identical(x,x.eval))  {
      ## Local polynomial estimation at evaluation points
      if(is.null(X)) {
        if(degree > 0) {
          X.poly <- poly(x,raw=poly.raw,degree=degree)
          X <- cbind(1,X.poly)
        } else {
          X <- matrix(1,nrow=length(x),ncol=1)
        }
      }
      if(is.null(X.eval)) {
        if(degree > 0) {
          if(!exists("X.poly")) X.poly <- poly(x,raw=poly.raw,degree=degree)
          X.eval <- cbind(1,predict(X.poly,x.eval))
        } else {
          X.eval <- matrix(1,nrow=length(x.eval),ncol=1)
        }
      }
      denom.x.eval <- h[2]*(if(is.infinite(x.ub)) 1 else pnorm((x.ub-x.eval)/h[2]) - (if(is.infinite(x.lb)) 0 else pnorm((x.lb-x.eval)/h[2])))
      denom.y.eval <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y.eval)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y.eval)/h[1])))
      f <- as.numeric(mcmapply(function(i){
        w <- NZD_pos(sqrt(pdf.kernel.bk(x.eval[i],x,h[2],x.lb,x.ub, denom=denom.x.eval[i])))
        beta.hat <- .lm.fit(X*w,pdf.kernel.bk(y.eval[i],y,h[1],y.lb,y.ub, denom=denom.y.eval[i])*w)$coefficients
        beta.hat%*%t(X.eval[i,,drop=FALSE])
      },seq_along(y.eval),mc.cores=optim.ksum.cores))
      if(any(f < 0)) return(-sqrt(.Machine$double.xmax))
    }
    if(is.null(X)) X <- if(degree>0) cbind(1,poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=length(x),ncol=1)
    f <- as.numeric(mcmapply(function(i){
      w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x,h[2],x.lb,x.ub, denom=denom.x[i])))
      beta.hat <- .lm.fit(X*w,pdf.kernel.bk(y[i],y,h[1],y.lb,y.ub, denom=denom.y[i])*w)$coefficients
      beta.hat%*%t(X[i,,drop=FALSE])
    },seq_along(y),mc.cores=optim.ksum.cores))
    if(any(f < 0)) return(-sqrt(.Machine$double.xmax))
  }

  if(bwmethod == "cv.ml") {
    if(degree==0) {
      f.loo <- as.numeric(mcmapply(function(i){
        pdf.kernel.bk.x<-pdf.kernel.bk(x[i],x[-i],h[2],x.lb,x.ub, denom=denom.x[i]);
        mean(pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i])*pdf.kernel.bk.x)/NZD_pos(mean(pdf.kernel.bk.x))
      },seq_along(y),mc.cores=optim.ksum.cores))
    } else {
      if(proper.cv) {
        if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
        if(is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(y.lb,extendrange(y,f=2)[2],length=n.integrate)
        if(!is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(extendrange(y,f=2)[1],y.ub,length=n.integrate)
        if(!is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(extendrange(y,f=2)[1],extendrange(y,f=2)[2],length=n.integrate)
        denom.y.seq <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y.seq)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y.seq)/h[1])))
        Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y, h[1], y.lb, y.ub, denom=denom.y.seq[i]),seq_along(y.seq))
        if(is.null(X)) X <- if(degree>0) cbind(1,poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=length(x),ncol=1)
        f.loo <- as.numeric(mcmapply(function(i){
          w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x[-i],h[2],x.lb,x.ub, denom=denom.x[i])))
          beta.hat <- .lm.fit(X[-i,,drop=FALSE]*w,cbind(pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i]),Y.seq.mat[-i,,drop=FALSE])*w)$coefficients
          f.loo <- X[i,,drop=FALSE]%*%beta.hat[,1,drop=FALSE]
          f.seq <- as.numeric(X[i,,drop=FALSE]%*%beta.hat[,2:dim(beta.hat)[2],drop=FALSE])
          f.seq[f.seq<0] <- 0
          f.loo[f.loo<0] <- 0
          f.loo/integrate.trapezoidal(y.seq,f.seq)[n.integrate]
        },seq_along(y),mc.cores=optim.ksum.cores))
      } else {
        if(is.null(X)) X <- if(degree>0) cbind(1,poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=length(x),ncol=1)
        f.loo <- as.numeric(mcmapply(function(i){
          w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x[-i],h[2],x.lb,x.ub, denom=denom.x[i])))
          beta.hat <- .lm.fit(X[-i,,drop=FALSE]*w,pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i])*w)$coefficients
          beta.hat%*%t(X[i,,drop=FALSE])
        },seq_along(y),mc.cores=optim.ksum.cores))
      }
    }
    return(sum(log.likelihood(f.loo,cv.penalty.method=cv.penalty.method,cv.penalty.cutoff=cv.penalty.cutoff,verbose=verbose,degree=degree,h=h)))
  } else {
    if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
    if(is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(y.lb,extendrange(y,f=2)[2],length=n.integrate)
    if(!is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(extendrange(y,f=2)[1],y.ub,length=n.integrate)
    if(!is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(extendrange(y,f=2)[1],extendrange(y,f=2)[2],length=n.integrate)
    denom.y.seq <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y.seq)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y.seq)/h[1])))
    Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y, h[1], y.lb, y.ub, denom=denom.y.seq[i]),seq_len(n.integrate))
    if(degree==0) {
      int.f.sq <- mcmapply(function(j){
        pdf.kernel.bk.x <- pdf.kernel.bk(x[j],x,h[2],x.lb,x.ub, denom=denom.x[j]);
        integrate.trapezoidal(y.seq,colMeans(Y.seq.mat*pdf.kernel.bk.x/NZD_pos(mean(pdf.kernel.bk.x)))^2)[n.integrate]
      },seq_along(y),mc.cores = optim.ksum.cores)
      f.loo <- as.numeric(mcmapply(function(i){
        pdf.kernel.bk.x<-pdf.kernel.bk(x[i],x[-i],h[2],x.lb,x.ub, denom=denom.x[i]);
        mean(pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i])*pdf.kernel.bk.x)/NZD_pos(mean(pdf.kernel.bk.x))
      },seq_along(y),mc.cores=optim.ksum.cores))
    } else {
      if(proper.cv) {
        if(is.null(X)) X <- if(degree>0) cbind(1,poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=length(x),ncol=1)
        foo <- (mcmapply(function(j){
          w <- NZD_pos(sqrt(pdf.kernel.bk(x[j],x,h[2],x.lb,x.ub, denom=denom.x[j])))
          beta.hat <- .lm.fit(X*w,cbind(pdf.kernel.bk(y[j],y[-j],h[1],y.lb,y.ub, denom=denom.y[j]),Y.seq.mat)*w)$coefficients
          f.loo <- X[j,,drop=FALSE]%*%beta.hat[,1,drop=FALSE]
          f.seq <- X[j,,drop=FALSE]%*%beta.hat[,2:(n.integrate+1)]
          f.seq[f.seq<0] <- 0
          f.seq <- f.seq/integrate.trapezoidal(y.seq,f.seq)[n.integrate]
          f.loo[f.loo<0] <- 0
          f.loo <- f.loo/integrate.trapezoidal(y.seq,f.seq)[n.integrate]
          list(int.f.sq=integrate.trapezoidal(y.seq,f.seq^2)[n.integrate], f.loo=f.loo)
        },seq_along(y),mc.cores = optim.ksum.cores))        
        int.f.sq <- sapply(foo, function(x) x$int.f.sq)
        f.loo <- sapply(foo, function(x) x$f.loo)
      } else {
        if(is.null(X)) X <- if(degree>0) cbind(1,poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=length(x),ncol=1)
        int.f.sq <- mcmapply(function(j){
          w <- NZD_pos(sqrt(pdf.kernel.bk(x[j],x,h[2],x.lb,x.ub, denom=denom.x[j])))
          beta.hat <- .lm.fit(X*w,Y.seq.mat*w)$coefficients
          integrate.trapezoidal(y.seq,(X[j,,drop=FALSE]%*%beta.hat)^2)[n.integrate]
        },seq_along(y),mc.cores = optim.ksum.cores)
        f.loo <- as.numeric(mcmapply(function(i){
          w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x[-i],h[2],x.lb,x.ub, denom=denom.x[i])))
          beta.hat <- .lm.fit(X[-i,,drop=FALSE]*w,pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i])*w)$coefficients
          beta.hat%*%t(X[i,,drop=FALSE])
        },seq_along(y),mc.cores=optim.ksum.cores))
      }
    }
    return(-(mean(int.f.sq)-2*mean(f.loo)))
  }
}

## bkcde.optim() manages the optimization process over multiple starts and degrees

bkcde.optim <- function(x=x,
                        y=y,
                        x.eval=x.eval,
                        y.eval=y.eval,
                        y.lb=y.lb,
                        y.ub=y.ub,
                        x.lb=x.lb,
                        x.ub=x.ub,
                        bwmethod=bwmethod,
                        cv.penalty.cutoff=cv.penalty.cutoff,
                        cv.penalty.method=cv.penalty.method,
                        degree.max=degree.max,
                        degree.min=degree.min,
                        n.integrate=n.integrate,
                        nmulti=nmulti,
                        optim.degree.cores=optim.degree.cores,
                        optim.ksum.cores=optim.ksum.cores,
                        optim.nmulti.cores=optim.nmulti.cores,
                        optim.sf.y.lb=optim.sf.y.lb,
                        optim.sf.x.lb=optim.sf.x.lb,
                        poly.raw=poly.raw,
                        proper.cv=proper.cv,
                        verbose=verbose,
                        ...) {
  if(degree.min < 0 || degree.max >= length(y)) stop("degree.min must lie in [0,1,...,",
                                                    length(y)-1,
                                                    "] (i.e., [0,1,dots, n-1]) in bkcde.optim()")
  if(degree.max < 0 || degree.max >= length(y)) stop("degree.max must lie in [0,1,...,",
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

  n <- length(y)
  lower <- c(optim.sf.y.lb*EssDee(y),optim.sf.x.lb*EssDee(x))*n^(-1/6)
  upper <- 10^(5)*EssDee(cbind(y,x))

  par.init <- matrix(NA,nmulti,2)
  par.init[1,] <- EssDee(cbind(y,x))*n^(-1/6)
  if(nmulti>1) {
    par.init[2:nmulti,] <- cbind(EssDee(y)*runif(nmulti-1,optim.sf.y.lb,10+optim.sf.y.lb),
                                 EssDee(x)*runif(nmulti-1,optim.sf.x.lb,10+optim.sf.x.lb))*n^(-1/6)
  }

  ## Search over polynomial degrees. Parallelized if optim.degree.cores > 1.
  ## Pre-calculating X and X.eval once per degree saves substantial computation 
  ## during the inner optimization multi-starts.
  degree.return <- mclapply(degree.min:degree.max, function(p) {
    if(p > 0) {
      poly.obj <- poly(x, raw=poly.raw, degree=p)
      X.p <- cbind(1, poly.obj)
      X.eval.p <- if(!identical(y,y.eval) || !identical(x,x.eval)) {
        cbind(1, predict(poly.obj, x.eval))
      } else NULL
    } else {
      ## degree=0 handles are constant, poly() would error. X is just a vector of ones.
      X.p <- matrix(1, nrow=length(x), ncol=1)
      X.eval.p <- if(!identical(y,y.eval) || !identical(x,x.eval)) {
        matrix(1, nrow=length(y.eval), ncol=1)
      } else NULL
    }
    
    nmulti.return <- mclapply(1:nmulti, function(i) {
      st <- system.time(optim.return <- optim(par=par.init[i,],
                                              fn=bkcde.optim.fn,
                                              x=x,
                                              y=y,
                                              x.eval=x.eval,
                                              y.eval=y.eval,
                                              y.lb=y.lb,
                                              y.ub=y.ub,
                                              x.lb=x.lb,
                                              x.ub=x.ub,
                                              poly.raw=poly.raw,
                                              bwmethod=bwmethod,
                                              cv.penalty.method=cv.penalty.method,
                                              cv.penalty.cutoff=cv.penalty.cutoff,
                                              degree=p,
                                              n.integrate=n.integrate,
                                              optim.ksum.cores=optim.ksum.cores,
                                              proper.cv=proper.cv,
                                              verbose=verbose,
                                              lower=lower,
                                              upper=upper,
                                              method="L-BFGS-B",
                                              control=list(fnscale = -1),
                                              X=X.p,
                                              X.eval=X.eval.p))
      optim.return$secs.optim <- st["elapsed"]
      optim.return$degree <- p
      optim.return$optim.par.init <- par.init[i,]
      optim.return
    },mc.cores = optim.nmulti.cores)
    optim.out <- nmulti.return[[which.max(sapply(nmulti.return, function(x) x$value))]]
    optim.out$value.vec <- sapply(nmulti.return, function(x) x$value)
    optim.out$degree.vec <- sapply(nmulti.return, function(x) x$degree)
    optim.out$optim.y.init.vec <- sapply(nmulti.return, function(x) x$optim.par.init[1])
    optim.out$optim.x.init.vec <- sapply(nmulti.return, function(x) x$optim.par.init[2])
    optim.out$convergence.vec <- sapply(nmulti.return, function(x) x$convergence)
    optim.out$secs.optim.vec <- sapply(nmulti.return, function(x) x$secs.optim)
    optim.out
  },mc.cores = optim.degree.cores)

  output.return <- degree.return[[which.max(sapply(degree.return, function(x) x$value))]]
  output.return$par.mat <- t(sapply(degree.return, function(x) x$par))
  output.return$value.vec <- sapply(degree.return, function(x) x$value)
  output.return$value.mat <- t(sapply(degree.return, function(x) x$value.vec))
  output.return$convergence.vec <- sapply(degree.return, function(x) x$convergence)
  output.return$convergence.mat <- t(sapply(degree.return, function(x) x$convergence.vec))
  output.return$degree.mat <- t(sapply(degree.return, function(x) x$degree.vec))
  output.return$secs.optim <- sapply(degree.return, function(x) x$secs.optim)
  output.return$secs.optim.mat <- t(sapply(degree.return, function(x) x$secs.optim.vec))
  output.return$optim.y.init.mat <- t(sapply(degree.return, function(x) x$optim.y.init.vec))
  output.return$optim.x.init.mat <- t(sapply(degree.return, function(x) x$optim.y.init.vec))
  return(output.return)
}

## sub.cv() implements sub-sampling cross-validation for large datasets

sub.cv <- function(x, y, 
                   n.sub = 300, 
                   progress = FALSE,
                   replace = FALSE,
                   resamples = 10,
                   resample.cores = parallel::detectCores(),
                   ...) {
  if(!is.numeric(x)) stop("x must be numeric in sub.cv()")
  if(!is.numeric(y)) stop("y must be numeric in sub.cv()")
  if(length(x) != length(y)) stop("length of x must be equal to length of y in sub.cv()")
  if(!is.numeric(n.sub)) stop("n.sub must be numeric in sub.cv()")
  if(n.sub < 100 || n.sub > length(y)) stop("n.sub must be at least 100 and less than the length of y in sub.cv()")
  if(!is.logical(progress)) stop("progress must be logical in sub.cv()")
  if(!is.logical(replace)) stop("replace must be logical in sub.cv()")
  if(replace == TRUE && resamples < 2) stop("resamples must be at least 2 when replace=TRUE in sub.cv()")
  if(n.sub==length(y) && replace==FALSE && resamples > 1) stop("taking resamples with replace=FALSE when n.sub=n results in identical samples")
  if(resamples < 1) stop("resamples must be at least 1 in sub.cv()")
  if(!is.numeric(resample.cores) || resample.cores < 1) stop("resample.cores must be a positive integer in sub.cv()")
  resample.cores <- as.integer(resample.cores)
  n <- length(y)
  sf.mat <- matrix(NA,nrow=resamples,ncol=2)
  degree.vec <- numeric()
  cv.vec <- numeric()
  
  ## Use L'Ecuyer RNG for reproducible parallel RNG streams and explicitly
  ## set child seed generation via mc.set.seed=TRUE. Save/restore RNG kind.
  .old.RNG <- RNGkind()
  on.exit(do.call(RNGkind, as.list(.old.RNG)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  
  if(resample.cores == 1) {
    ## Serial resampling (reliable, identical appearance to v1.30 progress bar)
    if(progress) pbb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta",
                                                   clear = TRUE,
                                                   force = TRUE,
                                                   total = resamples)
    if(progress) pbb$tick(0)
    sub.results <- vector("list", resamples)
    for(j in seq_len(resamples)) {
      ii <- sample(n, size=n.sub, replace=replace)
      bkcde.out <- bkcde(x=x[ii], y=y[ii], proper=FALSE, cv.only=TRUE, ...)
      sf <- bkcde.out$h/(EssDee(cbind(y[ii], x[ii])) * n.sub^(-1/6))
      sub.results[[j]] <- list(sf=sf, degree=bkcde.out$degree, cv=bkcde.out$value)
      if(progress) pbb$tick()
    }
  } else {
    if(progress) cat("Resampling: 0%\n")
    if(progress) {
      ## Use pbmclapply for parallel resampling (fast) and use the 'txt' style
      ## (txtProgressBar layout) so the progress output is more stable and similar
      ## to the original pbb appearance while preserving speed.
      sub.results <- mclapply.progress(seq_len(resamples), function(j) {
        ii <- sample(n, size=n.sub, replace=replace)
        bkcde.out <- bkcde(x=x[ii], y=y[ii], proper=FALSE, cv.only=TRUE, ...)
        sf <- bkcde.out$h/(EssDee(cbind(y[ii], x[ii])) * n.sub^(-1/6))
        res <- list(sf=sf, degree=bkcde.out$degree, cv=bkcde.out$value)
        return(res)
      }, mc.cores = resample.cores, mc.set.seed = TRUE, mc.style = "txt", mc.substyle = 3, ignore.interactive = TRUE, progress = TRUE)
    } else {
      ## No progress requested: use fast pbmclapply without rendering progress
      sub.results <- mclapply.progress(seq_len(resamples), function(j) {
        ii <- sample(n, size=n.sub, replace=replace)
        bkcde.out <- bkcde(x=x[ii], y=y[ii], proper=FALSE, cv.only=TRUE, ...)
        sf <- bkcde.out$h/(EssDee(cbind(y[ii], x[ii])) * n.sub^(-1/6))
        res <- list(sf=sf, degree=bkcde.out$degree, cv=bkcde.out$value)
        return(res)
      }, mc.cores = resample.cores, mc.set.seed = TRUE, progress = FALSE)
    }
  }
  
  for(j in seq_len(resamples)) {
    sf.mat[j,] <- sub.results[[j]]$sf
    degree.vec[j] <- sub.results[[j]]$degree
    cv.vec[j] <- sub.results[[j]]$cv
  }
  h.mat <- sweep(sf.mat,2,EssDee(cbind(y,x))*n^(-1/6),"*")
  degree <- min(find_mode(degree.vec))
  h.median <- apply(h.mat[degree.vec==degree,,drop=FALSE],2,median)
  h.mean <- apply(h.mat[degree.vec==degree,,drop=FALSE],2,mean)
  if(length(degree.vec[degree.vec==degree]) < 4) {
    h.covMcd <- h.mean
  } else {
    h.covMcd <- robustbase::covMcd(h.mat[degree.vec==degree,,drop=FALSE])$center
  }
  h.ml <- (h.mat[degree.vec==degree,,drop=FALSE])[which.max(cv.vec[degree.vec==degree]),,drop=FALSE]
  return(list(cv.vec=cv.vec,
              degree=degree,
              degree.modal.length=length(degree.vec[degree.vec==degree]),
              degree.vec=degree.vec,
              h.covMcd=h.covMcd,
              h.mat=h.mat,
              h.mean=h.mean,
              h.median=h.median,
              h.ml=h.ml,
              sf.mat=sf.mat))
}
