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

## bkcde.optim.fn() is the combined unbinned and linear binned cross-validation objective function Jan 21 Gemini enhanced

bkcde.optim.fn <- function(h=NULL, x=NULL, y=NULL, x.eval=NULL, y.eval=NULL,
                           y.lb=NULL, y.ub=NULL, x.lb=NULL, x.ub=NULL,
                           poly.raw=FALSE, degree=NULL, n.integrate=NULL,
                           optim.ksum.cores=1, cv.penalty.method=NULL,
                           cv.penalty.cutoff=NULL, verbose=FALSE,
                           bwmethod=NULL, proper.cv=NULL, X=NULL, X.eval=NULL,
                           optim.fn.type = c("linear-binning", "unbinned"),
                           grid.n = 100) {

  optim.fn.type <- match.arg(optim.fn.type)

  # --- 1. Strict Input Validation (Original) ---
  if(y.lb >= y.ub) stop("y.lb must be less than y.ub in bkcde.optim.fn()")
  if(x.lb >= x.ub) stop("x.lb must be less than x.ub in bkcde.optim.fn()")
  if(is.null(x) || is.null(y)) stop("must provide x and y in bkcde.optim.fn()")
  if(is.null(degree)) stop("must provide degree in bkcde.optim.fn()")
  if(optim.ksum.cores < 1) stop("optim.ksum.cores must be at least 1")
  if(degree < 0 || degree >= length(y)) stop("degree error")

  n.obs <- length(y)
  # Original Denominators for Boundary Kernels
  denom.x <- h[2]*(if(is.infinite(x.ub)) 1 else pnorm((x.ub-x)/h[2]) - (if(is.infinite(x.lb)) 0 else pnorm((x.lb-x)/h[2])))
  denom.y <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y)/h[1])))

  # --- 2. Non-Negativity Penalty (Original) ---
  if(cv.penalty.method=="nonneg" && degree > 0) {
    if(is.null(X)) X <- if(degree > 0) cbind(1, poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=n.obs,ncol=1)
    f_check_orig <- as.numeric(mcmapply(function(i){
      w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x,h[2],x.lb,x.ub, denom=denom.x[i])))
      beta.hat <- .lm.fit(X*w,pdf.kernel.bk(y[i],y,h[1],y.lb,y.ub, denom=denom.y[i])*w)$coefficients
      beta.hat%*%t(X[i,,drop=FALSE])
    },seq_along(y),mc.cores=optim.ksum.cores))
    if(any(f_check_orig < 0)) return(-sqrt(.Machine$double.xmax))
  }

  # --- 3. PATH: UNBINNED (Original Logic) ---
  if (optim.fn.type == "unbinned") {
    if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
    else y.seq <- seq(extendrange(y,f=2)[1],extendrange(y,f=2)[2],length=n.integrate)
    
    denom.y.seq <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y.seq)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y.seq)/h[1])))
    Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y, h[1], y.lb, y.ub, denom=denom.y.seq[i]),seq_along(y.seq))
    if(is.null(X)) X <- if(degree>0) cbind(1,poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=n.obs,ncol=1)

    if(bwmethod == "cv.ml") {
      f.loo <- as.numeric(mcmapply(function(i){
        w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x[-i],h[2],x.lb,x.ub, denom=denom.x[i])))
        if(proper.cv) {
          beta.hat <- .lm.fit(X[-i,,drop=FALSE]*w,cbind(pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i]),Y.seq.mat[-i,,drop=FALSE])*w)$coefficients
          f.loo.val <- max(0, X[i,,drop=FALSE]%*%beta.hat[,1])
          f.seq <- pmax(0, as.numeric(X[i,,drop=FALSE]%*%beta.hat[,2:ncol(beta.hat)]))
          return(f.loo.val/NZD_pos(integrate.trapezoidal(y.seq,f.seq)[n.integrate]))
        } else {
          beta.hat <- .lm.fit(X[-i,,drop=FALSE]*w,pdf.kernel.bk(y[i],y[-i],h[1],y.lb,y.ub, denom=denom.y[i])*w)$coefficients
          return(beta.hat%*%t(X[i,,drop=FALSE]))
        }
      },seq_along(y),mc.cores=optim.ksum.cores))
      val <- sum(log.likelihood(f.loo,cv.penalty.method,cv.penalty.cutoff,verbose,degree,h))
    } else {
      # LSCV Unbinned
      foo <- mcmapply(function(j){
        w <- NZD_pos(sqrt(pdf.kernel.bk(x[j],x,h[2],x.lb,x.ub, denom=denom.x[j])))
        beta.hat <- .lm.fit(X*w,Y.seq.mat*w)$coefficients
        int.f.sq <- integrate.trapezoidal(y.seq,(X[j,,drop=FALSE]%*%beta.hat)^2)[n.integrate]
        w_loo <- NZD_pos(sqrt(pdf.kernel.bk(x[j],x[-j],h[2],x.lb,x.ub, denom=denom.x[j])))
        beta.loo <- .lm.fit(X[-j,,drop=FALSE]*w_loo,pdf.kernel.bk(y[j],y[-j],h[1],y.lb,y.ub, denom=denom.y[j])*w_loo)$coefficients
        return(list(int.f.sq=int.f.sq, f.loo=beta.loo%*%t(X[j,,drop=FALSE])))
      },seq_along(y),mc.cores = optim.ksum.cores)
      val <- -(mean(unlist(foo[1,])) - 2*mean(unlist(foo[2,])))
    }
  }

  # --- 4. PATH: LINEAR BINNING (High Speed) ---
  else {
    x.grid <- seq(x.lb, x.ub, length.out = grid.n)
    y.grid <- seq(y.lb, y.ub, length.out = grid.n)
    delta.x <- x.grid[2] - x.grid[1]; delta.y <- y.grid[2] - y.grid[1]
    
    x.idx <- pmin(pmax(ceiling((x - x.lb) / delta.x), 1), grid.n)
    y.idx <- pmin(pmax(ceiling((y - y.lb) / delta.y), 1), grid.n)
    counts <- as.matrix(table(factor(x.idx, levels=1:grid.n), factor(y.idx, levels=1:grid.n)))
    
    active <- which(counts > 0, arr.ind = TRUE)
    x.act <- x.grid[active[,1]]; y.act <- y.grid[active[,2]]; w.act <- counts[active]
    X.act <- if(degree > 0) cbind(1, poly(x.act, degree=degree, raw=poly.raw)) else matrix(1, length(x.act), 1)

    y.seq <- if(is.finite(y.lb) && is.finite(y.ub)) seq(y.lb,y.ub,length=n.integrate) else 
             seq(extendrange(y,f=2)[1],extendrange(y,f=2)[2],length=n.integrate)
    denom.y.seq <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y.seq)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y.seq)/h[1])))
    Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y.act, h[1], y.lb, y.ub, denom=denom.y.seq[i]), seq_along(y.seq))

    results <- mcmapply(function(i) {
      denom.xi <- h[2] * (pnorm((x.ub - x.act[i])/h[2]) - pnorm((x.lb - x.act[i])/h[2]))
      k.weights <- pdf.kernel.bk(x.act[i], x.act, h[2], x.lb, x.ub, denom=denom.xi)
      w.loo <- w.act; w.loo[i] <- pmax(1e-7, w.loo[i] - 1)
      W_vec <- sqrt(k.weights * w.loo)
      
      if (bwmethod == "cv.ml") {
        denom.yi <- h[1] * (pnorm((y.ub - y.act[i])/h[1]) - pnorm((y.lb - y.act[i])/h[1]))
        target_y <- pdf.kernel.bk(y.act[i], y.act, h[1], y.lb, y.ub, denom=denom.yi)
        if (proper.cv) {
          beta <- .lm.fit(X.act * W_vec, cbind(target_y, Y.seq.mat) * (W_vec * w.loo))$coefficients
          f_loo <- max(0, sum(X.act[i,] * beta[,1]))
          f_seq <- pmax(0, as.numeric(X.act[i,] %*% beta[, 2:ncol(beta)]))
          return(f_loo / NZD_pos(integrate.trapezoidal(y.seq, f_seq)[n.integrate]))
        } else {
          return(sum(.lm.fit(X.act * W_vec, (target_y * w.loo) * W_vec)$coefficients * X.act[i,]))
        }
      } else {
        # LSCV
        beta_seq <- .lm.fit(X.act * W_vec, (Y.seq.mat * w.loo) * W_vec)$coefficients
        f_seq <- if(degree==0) colMeans(Y.seq.mat * k.weights)/NZD_pos(mean(k.weights)) else X.act[i,] %*% beta_seq
        denom.yi <- h[1] * (pnorm((y.ub - y.act[i])/h[1]) - pnorm((y.lb - y.act[i])/h[1]))
        target_y <- pdf.kernel.bk(y.act[i], y.act, h[1], y.lb, y.ub, denom=denom.yi)
        f_loo <- sum(.lm.fit(X.act * W_vec, (target_y * w.loo) * W_vec)$coefficients * X.act[i,])
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

  if (verbose) cat("Type:", optim.fn.type, "h:", h, "Objective:", val, "\n")
  return(val)
}

## bkcde.optim.fn() is the linear binned cross-validation objective function Jan 21 Gemini enhanced

bkcde.optim.fn.linear.binned <- function(h=NULL, x=NULL, y=NULL, x.eval=NULL, y.eval=NULL,
                           y.lb=NULL, y.ub=NULL, x.lb=NULL, x.ub=NULL,
                           poly.raw=FALSE, degree=NULL, n.integrate=NULL,
                           optim.ksum.cores=1, cv.penalty.method=NULL,
                           cv.penalty.cutoff=NULL, verbose=FALSE,
                           bwmethod=NULL, proper.cv=NULL, X=NULL, X.eval=NULL,
                           grid.n = 100) {

  # --- 1. Original Validations & Stops ---
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
  if(degree < 0 || degree >= length(y)) stop(paste0("degree must lie in [0,1,...,",length(y)-1,"] in bkcde.optim.fn()"))

  n.obs <- length(y)
  denom.x <- h[2]*(if(is.infinite(x.ub)) 1 else pnorm((x.ub-x)/h[2]) - (if(is.infinite(x.lb)) 0 else pnorm((x.lb-x)/h[2])))
  denom.y <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y)/h[1])))

  # --- 2. Original Non-Negativity Check ---
  if(cv.penalty.method=="nonneg" && degree > 0) {
    if(!identical(y,y.eval) || !identical(x,x.eval))  {
      if(is.null(X)) X <- if(degree > 0) cbind(1, poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=n.obs,ncol=1)
      if(is.null(X.eval)) {
        X.poly.temp <- poly(x,raw=poly.raw,degree=degree)
        X.eval <- if(degree > 0) cbind(1,predict(X.poly.temp,x.eval)) else matrix(1,nrow=length(x.eval),ncol=1)
      }
      denom.x.eval <- h[2]*(if(is.infinite(x.ub)) 1 else pnorm((x.ub-x.eval)/h[2]) - (if(is.infinite(x.lb)) 0 else pnorm((x.lb-x.eval)/h[2])))
      denom.y.eval <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y.eval)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y.eval)/h[1])))
      f_check <- as.numeric(mcmapply(function(i){
        w <- NZD_pos(sqrt(pdf.kernel.bk(x.eval[i],x,h[2],x.lb,x.ub, denom=denom.x.eval[i])))
        beta.hat <- .lm.fit(X*w,pdf.kernel.bk(y.eval[i],y,h[1],y.lb,y.ub, denom=denom.y.eval[i])*w)$coefficients
        beta.hat%*%t(X.eval[i,,drop=FALSE])
      },seq_along(y.eval),mc.cores=optim.ksum.cores))
      if(any(f_check < 0)) return(-sqrt(.Machine$double.xmax))
    }
    if(is.null(X)) X <- if(degree>0) cbind(1,poly(x,raw=poly.raw,degree=degree)) else matrix(1,nrow=n.obs,ncol=1)
    f_check_orig <- as.numeric(mcmapply(function(i){
      w <- NZD_pos(sqrt(pdf.kernel.bk(x[i],x,h[2],x.lb,x.ub, denom=denom.x[i])))
      beta.hat <- .lm.fit(X*w,pdf.kernel.bk(y[i],y,h[1],y.lb,y.ub, denom=denom.y[i])*w)$coefficients
      beta.hat%*%t(X[i,,drop=FALSE])
    },seq_along(y),mc.cores=optim.ksum.cores))
    if(any(f_check_orig < 0)) return(-sqrt(.Machine$double.xmax))
  }

  # --- 3. Setup Grids & Binning Logic ---
  x.grid <- seq(x.lb, x.ub, length.out = grid.n)
  y.grid <- seq(y.lb, y.ub, length.out = grid.n)
  delta.x <- x.grid[2] - x.grid[1]
  delta.y <- y.grid[2] - y.grid[1]

  x.idx <- pmin(pmax(ceiling((x - x.lb) / delta.x), 1), grid.n)
  y.idx <- pmin(pmax(ceiling((y - y.lb) / delta.y), 1), grid.n)
  counts <- as.matrix(table(factor(x.idx, levels=1:grid.n), factor(y.idx, levels=1:grid.n)))
  
  active <- which(counts > 0, arr.ind = TRUE)
  x.act <- x.grid[active[,1]]; y.act <- y.grid[active[,2]]; w.act <- counts[active]
  num.act <- length(x.act)
  X.act <- if(degree > 0) cbind(1, poly(x.act, degree=degree, raw=poly.raw)) else matrix(1, num.act, 1)

  # --- 4. Integration Setup (Handling Infinite Boundaries) ---
  if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
  else if(is.finite(y.lb)) y.seq <- seq(y.lb,extendrange(y,f=2)[2],length=n.integrate)
  else if(is.finite(y.ub)) y.seq <- seq(extendrange(y,f=2)[1],y.ub,length=n.integrate)
  else y.seq <- seq(extendrange(y,f=2)[1],extendrange(y,f=2)[2],length=n.integrate)
  
  denom.y.seq <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y.seq)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y.seq)/h[1])))
  Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y.act, h[1], y.lb, y.ub, denom=denom.y.seq[i]), seq_along(y.seq))

  # --- 5. Core Loop (Active Bins) ---
  results <- mcmapply(function(i) {
    denom.xi <- h[2] * (pnorm((x.ub - x.act[i])/h[2]) - pnorm((x.lb - x.act[i])/h[2]))
    k.weights <- pdf.kernel.bk(x.act[i], x.act, h[2], x.lb, x.ub, denom=denom.xi)
    
    w.loo <- w.act
    w.loo[i] <- pmax(1e-7, w.loo[i] - 1) # Prevent zero weights for .lm.fit
    W_vec <- sqrt(k.weights * w.loo)
    
    if (bwmethod == "cv.ml") {
      denom.yi <- h[1] * (pnorm((y.ub - y.act[i])/h[1]) - pnorm((y.lb - y.act[i])/h[1]))
      target_y <- pdf.kernel.bk(y.act[i], y.act, h[1], y.lb, y.ub, denom=denom.yi)
      
      if (proper.cv) {
        targets <- cbind(target_y, Y.seq.mat)
        beta <- .lm.fit(X.act * W_vec, targets * (W_vec * w.loo))$coefficients
        f_loo_val <- max(0, sum(X.act[i,] * beta[,1]))
        f_seq <- pmax(0, as.numeric(X.act[i,] %*% beta[, 2:ncol(beta)]))
        return(f_loo_val / NZD_pos(integrate.trapezoidal(y.seq, f_seq)[n.integrate]))
      } else {
        beta <- .lm.fit(X.act * W_vec, (target_y * w.loo) * W_vec)$coefficients
        return(sum(beta * X.act[i,]))
      }
    } else {
      # LSCV Logic
      beta_seq <- .lm.fit(X.act * W_vec, (Y.seq.mat * w.loo) * W_vec)$coefficients
      f_seq <- if(degree==0) colMeans(Y.seq.mat * k.weights)/NZD_pos(mean(k.weights)) else X.act[i,] %*% beta_seq
      int_f_sq <- integrate.trapezoidal(y.seq, f_seq^2)[n.integrate]
      
      denom.yi <- h[1] * (pnorm((y.ub - y.act[i])/h[1]) - pnorm((y.lb - y.act[i])/h[1]))
      target_y <- pdf.kernel.bk(y.act[i], y.act, h[1], y.lb, y.ub, denom=denom.yi)
      beta_loo <- .lm.fit(X.act * W_vec, (target_y * w.loo) * W_vec)$coefficients
      f_loo_val <- sum(beta_loo * X.act[i,])
      return(c(int_f_sq, f_loo_val))
    }
  }, seq_along(x.act), mc.cores = optim.ksum.cores, SIMPLIFY = FALSE)

  # --- 6. Final Aggregation ---
  if (bwmethod == "cv.ml") {
    f.loo.vec <- unlist(results)
    # Mirroring original log-likelihood replication
    val <- sum(log.likelihood(rep(f.loo.vec, w.act), cv.penalty.method, cv.penalty.cutoff, verbose, degree, h))
  } else {
    res.mat <- do.call(rbind, results)
    avg_int_f_sq <- sum(res.mat[,1] * w.act) / n.obs
    avg_f_loo <- sum(res.mat[,2] * w.act) / n.obs
    val <- -(avg_int_f_sq - 2 * avg_f_loo)
  }

  if (verbose) cat("Bandwidth h:", h, "Objective:", val, "\n")
  return(val)
}

bkcde.optim.fn.unbinned <- function(h=NULL,
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
                           optim.ksum.cores=1, # Recommended: set to 1
                           cv.penalty.method=NULL,
                           cv.penalty.cutoff=NULL,
                           verbose=FALSE,
                           bwmethod=NULL,
                           proper.cv=NULL,
                           X=NULL,
                           X.eval=NULL) {
  
  ## 1. Fast boundary checks
  if(y.lb >= y.ub || x.lb >= x.ub) stop("Boundaries are inconsistent.")
  
  n <- length(y)
  
  ## 2. Pre-calculate denominators (Vectorized)
  ## These depend on h, so they must stay inside the optim function
  denom.x <- h[2] * (if(is.infinite(x.ub)) 1 else pnorm((x.ub - x)/h[2]) - 
                       (if(is.infinite(x.lb)) 0 else pnorm((x.lb - x)/h[2])))
  denom.y <- h[1] * (if(is.infinite(y.ub)) 1 else pnorm((y.ub - y)/h[1]) - 
                       (if(is.infinite(y.lb)) 0 else pnorm((y.lb - y)/h[1])))

  ## 3. Non-negativity Check (Optimization)
  if(cv.penalty.method == "nonneg" && degree > 0) {
    # Helper to avoid code duplication
    check_neg <- function(eval_x, eval_y, eval_X, d_x, d_y) {
      f_vals <- sapply(seq_along(eval_y), function(i) {
        w <- NZD_pos(sqrt(pdf.kernel.bk(eval_x[i], x, h[2], x.lb, x.ub, denom = d_x[i])))
        # Use .lm.fit for speed, but avoid nested parallelization
        beta.hat <- .lm.fit(X * w, pdf.kernel.bk(eval_y[i], y, h[1], y.lb, y.ub, denom = d_y[i]) * w)$coefficients
        sum(beta.hat * eval_X[i, ])
      })
      return(any(f_vals < 0))
    }

    if(!identical(y, y.eval) || !identical(x, x.eval)) {
      denom.x.ev <- h[2]*(if(is.infinite(x.ub)) 1 else pnorm((x.ub-x.eval)/h[2]) - (if(is.infinite(x.lb)) 0 else pnorm((x.lb-x.eval)/h[2])))
      denom.y.ev <- h[1]*(if(is.infinite(y.ub)) 1 else pnorm((y.ub-y.eval)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb-y.eval)/h[1])))
      if(check_neg(x.eval, y.eval, X.eval, denom.x.ev, denom.y.ev)) return(-sqrt(.Machine$double.xmax))
    }
    if(check_neg(x, y, X, denom.x, denom.y)) return(-sqrt(.Machine$double.xmax))
  }

  ## 4. Cross-Validation Logic
  if(bwmethod == "cv.ml") {
    
    # Handle the y.seq integration grid outside the loop
    if(proper.cv) {
      y.seq <- if(is.finite(y.lb) && is.finite(y.ub)) seq(y.lb, y.ub, length = n.integrate) else {
        rng <- extendrange(y, f = 2)
        seq(if(is.finite(y.lb)) y.lb else rng[1], if(is.finite(y.ub)) y.ub else rng[2], length = n.integrate)
      }
      denom.y.seq <- h[1] * (if(is.infinite(y.ub)) 1 else pnorm((y.ub - y.seq)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb - y.seq)/h[1])))
      Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y, h[1], y.lb, y.ub, denom = denom.y.seq[i]), seq_along(y.seq))
    }

    # LOO Loop (Optimized with lapply/sapply instead of mcmapply)
    f.loo <- sapply(1:n, function(i) {
      w <- NZD_pos(sqrt(pdf.kernel.bk(x[i], x[-i], h[2], x.lb, x.ub, denom = denom.x[i])))
      
      if (proper.cv) {
        # Joint fit for efficiency
        target <- cbind(pdf.kernel.bk(y[i], y[-i], h[1], y.lb, y.ub, denom = denom.y[i]), Y.seq.mat[-i, , drop = FALSE])
        beta.hat <- .lm.fit(X[-i, , drop = FALSE] * w, target * w)$coefficients
        
        f_val <- sum(X[i, ] * beta.hat[, 1])
        f_seq <- as.numeric(X[i, , drop = FALSE] %*% beta.hat[, 2:ncol(beta.hat), drop = FALSE])
        
        f_val <- max(0, f_val)
        f_seq[f_seq < 0] <- 0
        return(f_val / integrate.trapezoidal(y.seq, f_seq)[n.integrate])
      } else {
        target <- pdf.kernel.bk(y[i], y[-i], h[1], y.lb, y.ub, denom = denom.y[i])
        beta.hat <- .lm.fit(X[-i, , drop = FALSE] * w, target * w)$coefficients
        return(sum(beta.hat * X[i, ]))
      }
    })
    
    return(sum(log.likelihood(f.loo, cv.penalty.method, cv.penalty.cutoff, verbose, degree, h)))

  } else {
    ## 5. LSCV Logic (Non-ML CV)
    # y.seq setup for integration
    y.seq <- if(is.finite(y.lb) && is.finite(y.ub)) seq(y.lb, y.ub, length = n.integrate) else {
      rng <- extendrange(y, f = 2)
      seq(if(is.finite(y.lb)) y.lb else rng[1], if(is.finite(y.ub)) y.ub else rng[2], length = n.integrate)
    }
    denom.y.seq <- h[1] * (if(is.infinite(y.ub)) 1 else pnorm((y.ub - y.seq)/h[1]) - (if(is.infinite(y.lb)) 0 else pnorm((y.lb - y.seq)/h[1])))
    Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y, h[1], y.lb, y.ub, denom = denom.y.seq[i]), seq_len(n.integrate))

    res <- lapply(1:n, function(i) {
      w <- NZD_pos(sqrt(pdf.kernel.bk(x[i], x, h[2], x.lb, x.ub, denom = denom.x[i])))
      
      if(degree == 0) {
        # Simpler logic for degree 0
        mean_w <- mean(w^2) # Since w is sqrt(kernel)
        f_loo_w <- NZD_pos(sqrt(pdf.kernel.bk(x[i], x[-i], h[2], x.lb, x.ub, denom = denom.x[i])))
        f_loo <- mean(pdf.kernel.bk(y[i], y[-i], h[1], y.lb, y.ub, denom = denom.y[i]) * f_loo_w^2) / NZD_pos(mean(f_loo_w^2))
        
        # Integration of f^2
        f_seq_vals <- colMeans(Y.seq.mat * (w^2)) / NZD_pos(mean_w)
        int_f_sq <- integrate.trapezoidal(y.seq, f_seq_vals^2)[n.integrate]
        return(list(int_f_sq = int_f_sq, f.loo = f.loo))
      } else {
        # Degree > 0 Local Poly
        if(proper.cv) {
          beta.hat <- .lm.fit(X * w, cbind(pdf.kernel.bk(y[i], y[-i], h[1], y.lb, y.ub, denom = denom.y[i]), Y.seq.mat) * w)$coefficients
          # Logic for scaling/proper CV
          f_loo <- max(0, sum(X[i, ] * beta.hat[, 1]))
          f_seq <- X[i, , drop = FALSE] %*% beta.hat[, 2:ncol(beta.hat)]
          f_seq[f_seq < 0] <- 0
          denom_int <- integrate.trapezoidal(y.seq, f_seq)[n.integrate]
          f_seq <- f_seq / denom_int
          return(list(int_f_sq = integrate.trapezoidal(y.seq, f_seq^2)[n.integrate], f.loo = f_loo / denom_int))
        } else {
          beta.hat_seq <- .lm.fit(X * w, Y.seq.mat * w)$coefficients
          int_f_sq <- integrate.trapezoidal(y.seq, (X[i, , drop=FALSE] %*% beta.hat_seq)^2)[n.integrate]
          
          w_loo <- NZD_pos(sqrt(pdf.kernel.bk(x[i], x[-i], h[2], x.lb, x.ub, denom = denom.x[i])))
          beta.hat_loo <- .lm.fit(X[-i, , drop=FALSE] * w_loo, pdf.kernel.bk(y[i], y[-i], h[1], y.lb, y.ub, denom = denom.y[i]) * w_loo)$coefficients
          return(list(int_f_sq = int_f_sq, f.loo = sum(beta.hat_loo * X[i, ])))
        }
      }
    })
    
    int.f.sq <- sapply(res, `[[`, "int_f_sq")
    f.loo <- sapply(res, `[[`, "f.loo")
    return(-(mean(int.f.sq) - 2 * mean(f.loo)))
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

  n <- length(y)
  lower <- c(optim.sf.y.lb*EssDee(y),optim.sf.x.lb*EssDee(x))*n^(-1/6)
  upper <- 10^(5)*EssDee(cbind(y,x))

    par.init <- matrix(NA,nmulti,2)
    seeds.used <- rep(NA, nmulti)
    # First start
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
      # Set the same seed as used for par.init for reproducibility
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
                                 lower=lower,
                                 upper=upper,
                                 method="L-BFGS-B",
                                 control=list(fnscale = -1),
                                 X=X.p,
                                 X.eval=X.eval.p))
      optim.return$secs.optim <- st["elapsed"]
      optim.return$degree <- p
      optim.return$optim.par.init <- par.init[i,]
      optim.return$seed <- seeds.used[i]
      optim.return
    },mc.cores = optim.nmulti.cores)
    optim.out <- nmulti.return[[which.max(sapply(nmulti.return, function(x) x$value))]]
    optim.out$value.vec <- sapply(nmulti.return, function(x) x$value)
    optim.out$degree.vec <- sapply(nmulti.return, function(x) x$degree)
    optim.out$optim.y.init.vec <- sapply(nmulti.return, function(x) x$optim.par.init[1])
    optim.out$optim.x.init.vec <- sapply(nmulti.return, function(x) x$optim.par.init[2])
    optim.out$convergence.vec <- sapply(nmulti.return, function(x) x$convergence)
    optim.out$secs.optim.vec <- sapply(nmulti.return, function(x) x$secs.optim)
      optim.out$seed.vec <- sapply(nmulti.return, function(x) x$seed)
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

