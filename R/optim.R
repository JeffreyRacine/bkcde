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

## bkcde.optim.fn() is the combined unbinned and linear binned cross-validation objective function

bkcde.optim.fn <- function(h=NULL, x=NULL, y=NULL, x.eval=NULL, y.eval=NULL,
                           y.lb=NULL, y.ub=NULL, x.lb=NULL, x.ub=NULL,
                           poly.raw=FALSE, degree=NULL, n.integrate=NULL,
                           optim.ksum.cores=1, cv.penalty.method=NULL,
                           cv.penalty.cutoff=NULL, verbose=FALSE,
                           bwmethod=NULL, proper.cv=NULL, X=NULL, X.eval=NULL,
                           X.act=NULL, y.seq=NULL,
                           cv.binned = FALSE,
                           n.binned = 100,
                           binned.data = NULL,
                           x.ub.finite = TRUE,
                           x.lb.finite = TRUE,
                           y.ub.finite = TRUE,
                           y.lb.finite = TRUE) {

  # --- 1. Validation ---
  if(y.lb >= y.ub) stop("y.lb must be less than y.ub in bkcde.optim.fn()")
  if(x.lb >= x.ub) stop("x.lb must be less than x.ub in bkcde.optim.fn()")
  if(!cv.binned && (is.null(x) || is.null(y))) stop("must provide x and y in bkcde.optim.fn()")
  if(cv.binned && is.null(binned.data)) stop("must provide binned.data in bkcde.optim.fn()")
  if(is.null(degree)) stop("must provide degree in bkcde.optim.fn()")

  n.obs <- if(!cv.binned) length(y) else binned.data$n.obs
  
  if(!cv.binned) {
    denom.x <- h[2]*(if(!x.ub.finite) 1 else pnorm((x.ub-x)/h[2]) - (if(x.lb.finite) pnorm((x.lb-x)/h[2]) else 0))
    denom.y <- h[1]*(if(!y.ub.finite) 1 else pnorm((y.ub-y)/h[1]) - (if(y.lb.finite) pnorm((y.lb-y)/h[1]) else 0))
  }

  # --- 2. Non-Negativity Penalty ---
  if(cv.penalty.method=="nonneg" && degree > 0 && !cv.binned) {
    if(is.null(X)) X <- cbind(1, poly(x,raw=poly.raw,degree=degree))
    f_check_orig <- as.numeric(mcmapply(function(i){
      w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x,h[2],x.lb,x.ub, denom=denom.x[i])))
      beta.hat <- .lm.fit(X*w,pdf.kernel.bk(y[i],y,h[1],y.lb,y.ub, denom=denom.y[i])*w)$coefficients
      beta.hat%*%t(X[i,,drop=FALSE])
    },seq_along(y),mc.cores=optim.ksum.cores))
    if(any(f_check_orig < 0)) return(-sqrt(.Machine$double.xmax))
  }

  # --- 3. PATH: UNBINNED ---
  if (!cv.binned) {
    if(is.null(X)) X <- if(degree>0) cbind(1,poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=n.obs,ncol=1)

    if(bwmethod == "cv.ml") {
      if(proper.cv) {
        ## Compute y.seq and Y.seq.mat only when needed for proper.cv
        if(is.null(y.seq)) {
          if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
          else y.seq <- seq(extendrange(y,f=2)[1],extendrange(y,f=2)[2],length=n.integrate)
        }
        denom.y.seq <- h[1]*(if(!y.ub.finite) 1 else pnorm((y.ub-y.seq)/h[1]) - (if(y.lb.finite) pnorm((y.lb-y.seq)/h[1]) else 0))
        Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y, h[1], y.lb, y.ub, denom=denom.y.seq[i]),seq_along(y.seq))
        
        f.loo <- as.numeric(mcmapply(function(i){
          w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x[-i],h[2],x.lb,x.ub, denom=denom.x[i])))
          beta.hat <- .lm.fit(X[-i,,drop=FALSE]*w,cbind(pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i]),Y.seq.mat[-i,,drop=FALSE])*w)$coefficients
          f.loo.val <- max(0, X[i,,drop=FALSE]%*%beta.hat[,1])
          f.seq <- pmax(0, as.numeric(X[i,,drop=FALSE]%*%beta.hat[,2:ncol(beta.hat)]))
          f.loo.val/NZD_pos(integrate.trapezoidal(y.seq,f.seq)[n.integrate])
        },seq_along(y),mc.cores=optim.ksum.cores))
      } else {
        f.loo <- as.numeric(mcmapply(function(i){
          w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x[-i],h[2],x.lb,x.ub, denom=denom.x[i])))
          beta.hat <- .lm.fit(X[-i,,drop=FALSE]*w,pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i])*w)$coefficients
          beta.hat%*%t(X[i,,drop=FALSE])
        },seq_along(y),mc.cores=optim.ksum.cores))
      }
      val <- sum(log.likelihood(f.loo,cv.penalty.method,cv.penalty.cutoff,verbose,degree,h))
    } else {
      ## For cv.ls, always need y.seq and Y.seq.mat
      if(is.null(y.seq)) {
        if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
        else y.seq <- seq(extendrange(y,f=2)[1],extendrange(y,f=2)[2],length=n.integrate)
      }
      denom.y.seq <- h[1]*(if(!y.ub.finite) 1 else pnorm((y.ub-y.seq)/h[1]) - (if(y.lb.finite) pnorm((y.lb-y.seq)/h[1]) else 0))
      Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y, h[1], y.lb, y.ub, denom=denom.y.seq[i]),seq_along(y.seq))
      
      foo <- mcmapply(function(j){
        w <- NZD_pos(sqrt(pdf.kernel.bk(x[j],x,h[2],x.lb,x.ub, denom=denom.x[j])))
        beta.hat <- .lm.fit(X*w,Y.seq.mat*w)$coefficients
        int.f.sq <- integrate.trapezoidal(y.seq,(X[j,,drop=FALSE]%*%beta.hat)^2)[n.integrate]
        w_loo <- NZD_pos(sqrt(pdf.kernel.bk(x[j],x[-j],h[2],x.lb,x.ub, denom=denom.x[j])))
        beta.loo <- .lm.fit(X[-j,,drop=FALSE]*w_loo,pdf.kernel.bk(y[j],y[-j],h[1],y.lb,y.ub, denom=denom.y[j])*w_loo)$coefficients
        list(int.f.sq=int.f.sq, f.loo=beta.loo%*%t(X[j,,drop=FALSE]))
      },seq_along(y),mc.cores = optim.ksum.cores)
      val <- -(mean(unlist(foo[1,])) - 2*mean(unlist(foo[2,])))
    }
  }

  # --- 4. PATH: LINEAR BINNING ---
  else {
    x.act <- binned.data$x.act
    y.act <- binned.data$y.act
    w.act <- binned.data$w.act
    
    ## Use the pre-computed matrix passed from bkcde.optim
    if(is.null(X.act)) {
       X.act <- if(degree > 0) cbind(1, poly(x.act, degree=degree, raw=poly.raw)) else matrix(1, length(x.act), 1)
    }

    ## Vectorized denominator calculations for speed
    denom.x.act <- h[2] * (if(!x.ub.finite) 1 else pnorm((x.ub - x.act)/h[2]) - (if(x.lb.finite) pnorm((x.lb - x.act)/h[2]) else 0))
    denom.y.act <- h[1] * (if(!y.ub.finite) 1 else pnorm((y.ub - y.act)/h[1]) - (if(y.lb.finite) pnorm((y.lb - y.act)/h[1]) else 0))

    ## Use pre-computed y.seq if provided, otherwise compute it
    if(is.null(y.seq)) {
      y.seq <- if(is.finite(y.lb) && is.finite(y.ub)) seq(y.lb,y.ub,length=n.integrate) else 
               seq(extendrange(y.act,f=2)[1],extendrange(y.act,f=2)[2],length=n.integrate)
    }
    denom.y.seq <- h[1]*(if(!y.ub.finite) 1 else pnorm((y.ub-y.seq)/h[1]) - (if(y.lb.finite) pnorm((y.lb-y.seq)/h[1]) else 0))
    Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y.act, h[1], y.lb, y.ub, denom=denom.y.seq[i]), seq_along(y.seq))

    results <- mcmapply(function(i) {
      k.weights <- pdf.kernel.bk(x.act[i], x.act, h[2], x.lb, x.ub, denom=denom.x.act[i])
      w.loo <- w.act; w.loo[i] <- pmax(1e-7, w.loo[i] - 1)
      W_vec <- sqrt(k.weights * w.loo)
      
      if (bwmethod == "cv.ml") {
        target_y <- pdf.kernel.bk(y.act[i], y.act, h[1], y.lb, y.ub, denom=denom.y.act[i])
        if (proper.cv) {
          beta <- .lm.fit(X.act * W_vec, cbind(target_y, Y.seq.mat) * (W_vec * w.loo))$coefficients
          f_loo <- max(0, sum(X.act[i,] * beta[,1]))
          f_seq <- pmax(0, as.numeric(X.act[i,] %*% beta[, 2:ncol(beta)]))
          return(f_loo / NZD_pos(integrate.trapezoidal(y.seq, f_seq)[n.integrate]))
        } else {
          return(sum(.lm.fit(X.act * W_vec, (target_y * w.loo) * W_vec)$coefficients * X.act[i,]))
        }
      } else {
        # Batch: Compute target_y early and combine both RHS for single .lm.fit call
        target_y <- pdf.kernel.bk(y.act[i], y.act, h[1], y.lb, y.ub, denom=denom.y.act[i])
        W_scaled <- W_vec * w.loo
        
        # Combine both RHS and solve once instead of twice
        rhs_both <- cbind((Y.seq.mat * w.loo) * W_vec, (target_y * w.loo) * W_vec)
        beta_both <- .lm.fit(X.act * W_vec, rhs_both)$coefficients
        
        # Extract individual results
        n.seq.cols <- ncol(Y.seq.mat)
        beta_seq <- beta_both[, 1:n.seq.cols, drop=FALSE]
        beta_loo <- beta_both[, n.seq.cols + 1]
        
        f_seq <- if(degree==0) colMeans(Y.seq.mat * k.weights)/NZD_pos(mean(k.weights)) else X.act[i,] %*% beta_seq
        f_loo <- sum(X.act[i,] * beta_loo)
        return(c(integrate.trapezoidal(y.seq, f_seq^2)[n.integrate], f_loo))
      }
    }, seq_along(x.act), mc.cores = optim.ksum.cores, SIMPLIFY = FALSE)

    if (bwmethod == "cv.ml") {
      val <- sum(w.act * log.likelihood(unlist(results), cv.penalty.method, cv.penalty.cutoff, verbose, degree, h))
    } else {
      res.mat <- do.call(rbind, results)
      val <- -(sum(res.mat[,1] * w.act)/n.obs - 2 * sum(res.mat[,2] * w.act)/n.obs)
    }
  }

  if (verbose) cat("Type:", ifelse(cv.binned, "binned", "unbinned"), "h:", h, "Objective:", val, "\n")
  return(val)
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
                        cv.binned=FALSE,
                        n.binned=100,
                        y.seq=NULL,
                        verbose=verbose,
                        seed=NULL,
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
  if(!is.numeric(n.binned) || length(n.binned) != 1 || n.binned < 2) stop("n.binned must be numeric and at least 2 in bkcde.optim()")

  n.binned <- as.integer(n.binned)

  n <- length(y)
  lower <- c(optim.sf.y.lb*EssDee(y),optim.sf.x.lb*EssDee(x))*n^(-1/6)
  upper <- 10^(5)*EssDee(cbind(y,x))
  
  ## Pre-compute bound finiteness flags once to avoid repeated is.infinite() checks
  x.ub.finite <- is.finite(x.ub)
  x.lb.finite <- is.finite(x.lb)
  y.ub.finite <- is.finite(y.ub)
  y.lb.finite <- is.finite(y.lb)
  
  ## Pre-bin data once if using linear-binning
  if(cv.binned) {
    binned.data <- bkcde.bin.data(x, y, x.lb, x.ub, y.lb, y.ub, n.binned)
  } else {
    binned.data <- NULL
  }

  ## Pre-compute y.seq once per optimization call when provided or derivable
  if(is.null(y.seq)) {
    y.seq <- if(is.finite(y.lb) && is.finite(y.ub)) {
      seq(y.lb, y.ub, length = n.integrate)
    } else {
      seq(extendrange(y, f = 2)[1], extendrange(y, f = 2)[2], length = n.integrate)
    }
  }

  par.init <- matrix(NA,nmulti,2)
  seeds.used <- rep(NA, nmulti)
  if (!is.null(seed)) {
    set.seed(seed + 1)
    seeds.used[1] <- seed + 1
  }
  par.init[1,] <- EssDee(cbind(y,x))*n^(-1/6)
  if(nmulti>1) {
    for (i in 2:nmulti) {
      if (!is.null(seed)) {
        set.seed(seed + i)
        seeds.used[i] <- seed + i
      }
      par.init[i,] <- c(EssDee(y)*runif(1,optim.sf.y.lb,10+optim.sf.y.lb),
                        EssDee(x)*runif(1,optim.sf.x.lb,10+optim.sf.x.lb))*n^(-1/6)
    }
  }

  ## Search over polynomial degrees.
  degree.return <- mclapply(degree.min:degree.max, function(p) {
    ## Pre-calculate design matrices for this specific degree
    if(p > 0) {
      poly.obj <- poly(x, raw=poly.raw, degree=p)
      X.p <- cbind(1, poly.obj)
      X.eval.p <- if(!identical(y,y.eval) || !identical(x,x.eval)) {
        cbind(1, predict(poly.obj, x.eval))
      } else NULL
      
      ## Pre-calculate the Active Design Matrix for binned data
      X.act.p <- if(!is.null(binned.data)) {
        cbind(1, predict(poly.obj, binned.data$x.act))
      } else NULL
    } else {
      X.p <- matrix(1, nrow=length(x), ncol=1)
      X.eval.p <- if(!identical(y,y.eval) || !identical(x,x.eval)) {
        matrix(1, nrow=length(y.eval), ncol=1)
      } else NULL
      X.act.p <- if(!is.null(binned.data)) {
        matrix(1, nrow=length(binned.data$x.act), ncol=1)
      } else NULL
    }
    
    nmulti.return <- mclapply(1:nmulti, function(i) {
      if (!is.null(seed)) set.seed(seeds.used[i])
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
                                 cv.binned=cv.binned,
                                 n.binned=n.binned,
                                 binned.data=binned.data,
                                 y.seq=y.seq,
                                 x.ub.finite=x.ub.finite,
                                 x.lb.finite=x.lb.finite,
                                 y.ub.finite=y.ub.finite,
                                 y.lb.finite=y.lb.finite,
                                 lower=lower,
                                 upper=upper,
                                 method="L-BFGS-B",
                                 control=list(fnscale = -1),
                                 X=X.p,
                                 X.eval=X.eval.p,
                                 X.act=X.act.p)) # Pass the pre-computed active matrix
      optim.return$secs.optim <- st["elapsed"]
      optim.return$degree <- p
      optim.return$optim.par.init <- par.init[i, ]
      optim.return$seed <- seeds.used[i]
      optim.return
    }, mc.cores = optim.nmulti.cores)
    
    optim.out <- nmulti.return[[which.max(sapply(nmulti.return, function(x) x$value))]]
    optim.out$value.vec <- sapply(nmulti.return, function(x) x$value)
    optim.out$degree.vec <- sapply(nmulti.return, function(x) x$degree)
    optim.out$optim.y.init.vec <- sapply(nmulti.return, function(x) x$optim.par.init[1])
    optim.out$optim.x.init.vec <- sapply(nmulti.return, function(x) x$optim.par.init[2])
    optim.out$convergence.vec <- sapply(nmulti.return, function(x) x$convergence)
    optim.out$secs.optim.vec <- sapply(nmulti.return, function(x) x$secs.optim)
    optim.out$seed.vec <- sapply(nmulti.return, function(x) x$seed)
    optim.out
  }, mc.cores = optim.degree.cores)

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
  output.return$optim.x.init.mat <- t(sapply(degree.return, function(x) x$optim.x.init.vec))
  return(output.return)
}

## sub.cv() implements sub-sampling cross-validation for large datasets

sub.cv <- function(x, y, 
                   n.sub = 300, 
                   progress = FALSE,
                   replace = FALSE,
                   resamples = 10, 
                   ...) {
  if(!is.numeric(x)) stop("x must be numeric in sub.cv()")
  if(!is.numeric(y)) stop("y must be numeric in sub.cv()")
  if(length(x) != length(y)) stop("length of x must be equal to length of y in sub.cv()")
  if(!is.numeric(n.sub)) stop("n.sub must be numeric in sub.cv()")
  if(n.sub < 100 | n.sub > length(y)) stop("n.sub must be at least 100 and less than the length of y in sub.cv()")
  if(!is.logical(progress)) stop("progress must be logical in sub.cv()")
  if(!is.logical(replace)) stop("replace must be logical in sub.cv()")
  ## If only 1 resample is specified it ought to be the original sample
  ## returned, so check
  if(replace == TRUE & resamples < 2) stop("resamples must be at least 2 when replace=TRUE in sub.cv()")
  if(n.sub==length(y) & replace==FALSE & resamples > 1) stop("taking resamples with replace=FALSE when n.sub=n results in identical samples")
  if(resamples < 1) stop("resamples must be at least 1 in sub.cv()")
  n <- length(y)
  sf.mat <- matrix(NA,nrow=resamples,ncol=2)
  degree.vec <- numeric()
  cv.vec <- numeric()
  if(progress) pbb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta",
                                                 clear = TRUE,
                                                 force = TRUE,
                                                 total = resamples)
  if(progress) cat("\rResampling: 0%")
  for(j in 1:resamples) {
    ii <- sample(n,size=n.sub,replace=replace)
    ## Since cross-validation in bkcde() appropriately deals with improper
    ## densities, and since we are only using cross-validation in this call, we
    ## set proper=FALSE. We retrieve the "scale factors" after removing scale
    ## and sample size factors.
    bkcde.out <- bkcde(x=x[ii],y=y[ii],proper=FALSE,cv.only=TRUE,...)
    sf.mat[j,] <- bkcde.out$h/(EssDee(cbind(y[ii],x[ii]))*n.sub^(-1/6))
    degree.vec[j] <- bkcde.out$degree
    cv.vec[j] <- bkcde.out$value
    if(progress) pbb$tick()
  }
  ## Compute "typical" column elements of h.mat after rescaling for larger
  ## sample size and scale of data.
  h.mat <- sweep(sf.mat,2,EssDee(cbind(y,x))*n^(-1/6),"*")
  ## We use robust "typical" measures of location for h and degree since,
  ## importantly, bandwidth properties differ with degree of polynomial (rates
  ## and values) and so it is not sensible to unconditionally return e.g. the
  ## median across all degrees and bandwidths. We first select the "typical"
  ## polynomial order (smallest degree mode) then take a robust measure of the
  ## "typical" vector of bandwidths corresponding to the typical polynomial
  ## order providing n > p+1 (min required by MCD). Note that the modal degrees,
  ## when > 1 exist, may not be contiguous hence taking the mean degree may be
  ## is ill-advised.
  degree <- min(find_mode(degree.vec))
  h.median <- apply(h.mat[degree.vec==degree,,drop=FALSE],2,median)
  h.mean <- apply(h.mat[degree.vec==degree,,drop=FALSE],2,mean)
  if(length(degree.vec[degree.vec==degree]) < 4) {
    h.covMcd <- h.mean
  } else {
    h.covMcd <- covMcd(h.mat[degree.vec==degree,,drop=FALSE])$center
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

