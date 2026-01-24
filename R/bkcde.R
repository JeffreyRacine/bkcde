## Description: This file contains the R functions used in the simulation study
## for implementing the boundary kernel conditional density estimator with
## cross-validated bandwidth and polynomial order selection.  The functions are
## used in the simulation study in the file mc.R, licensed under
## GPL-3.0-or-later; written by Jeffrey Racine <racinej@mcmaster.ca>

## mclapply() and mcmapply() from the parallel package are used throughout
## instead of lapply() and mapply() to allow for multi-core computing, and when
## optim.ksum.cores, optim.degree.cores, or optim.nmulti.cores are set to a
## value greater than 1, the number of cores used in the parallel processing is
## set to the value of the respective argument.  But when these are set to 1,
## the number of cores used in the parallel processing is set to 1, i.e., serial
## processing occurs exactly as if lapply() and mapply() were being used. If
## verbose=TRUE is enabled, it appears warnings are most likely to appear
## immediately running things in serial mode.

## The functions are now distributed across utils.R, kernels.R, and
## optim.R. This file contains the main S3 methods.



## bckde() and bkcde.default() compute the conditional density \hat f(y|x)
## where, if no bandwidth is provided, either likelihood or least-squares
## cross-validation is used to select the bandwidths and polynomial order via
## numerical optimization with nmulti restarts by default and polynomial orders
## 0,1,...,degree.max by default. Restarting is used in an attempt to maximize
## the likelihood and (hopefully) avoid local optima. This function supports
## local polynomial orders [0,1,...,n-1] where n is the number of sample
## realizations (raw polynomials or orthogonal polynomials can be used and
## appear to provide identical results for modest degree.max). For large samples
## the function can be computationally intensive so we include a sub-sampling
## cross-validation procedure following Racine (1993) (cv="sub") which can be
## used to reduce computation time for large samples, say of the order 10^7,
## which can be handled a few minutes on a modern processor. Note that we rely
## on the mcmapply and mclapply functions which rely on "forking" that is not
## currently available on Windows.  The function returns a list of class "bkcde"
## with the following components: convergence.mat, convergence.vec, convergence,
## cv, degree.mat, degree.max, degree.min, degree, fitted.cores,
## f.yx.integral.post, f.yx.integral.pre.neg, f.yx.integral, f, h.mat, h,
## optim.ksum.cores, optim.degree.cores, optim.nmulti.cores, optimize,
## proper.cores, proper, secs.elapsed, secs.estimate, secs.optim.mat, value.mat,
## value.vec, value, x.eval, x.lb, x.ub, x, y.eval, y.lb, y.ub, y. S3 methods
## for the class "bkcde" include fitting, plotting, and predicting the
## conditional density estimate.

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
                          bwmethod=c("cv.ml","cv.ls"),
                          bwscaling = FALSE,
                          cv=c("auto","full","sub"),
                          cv.auto.threshold=5000,
                          cv.binned=FALSE,
                          cv.only=FALSE,
                          cv.penalty.cutoff=.Machine$double.xmin,
                          cv.penalty.method=c("smooth","constant","trim","nonneg"),
                          degree=NULL,
                          degree.max=3,
                          degree.min=0,
                          fitted.cores=detectCores(),
                          integrate.erf=1,
                          n.binned=100,
                          n.grid=25,
                          n.integrate=100,
                          n.sub=300,
                          nmulti=3,
                          optim.cores=c("auto","manual"),
                          optim.degree.cores=NULL,
                          optim.ksum.cores=NULL,
                          optim.ksum.auto.thresholds=c(2000,4000),
                          optim.ksum.auto.max=4,
                          optim.nmulti.cores=NULL,
                          optim.sf.x.lb=0.5,
                          optim.sf.y.lb=0.5,
                          poly.raw=FALSE,
                          progress=FALSE,
                          proper=FALSE,
                          proper.cores=detectCores(),
                          proper.cv=FALSE,
                          resamples=10,
                          seed=42,
                          verbose=FALSE,
                          warnings=TRUE,
                          x.erf=0,
                          y.erf=0,
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
  if(!is.logical(cv.only)) stop("cv.only must be logical in bkcde()")
  if(!is.null(h) & cv.only) stop("cannot provide h when cv.only=TRUE in bkcde()")
  if(!is.null(y.eval) & is.null(x.eval)) stop("must provide x.eval in bkcde() when y.eval is not NULL")
  ## We set x.eval and y.eval to short sequences if they are not provided to
  ## avoid excessive computation with large samples when only x and y are
  ## provided (essentially what would be the default for n.grid in
  ## plot.bkcde()). Of course these can be changed, and plot will override them
  ## if desired, as will predict.bkcde().
  optimize <- ifelse(is.null(h),TRUE,FALSE)
  if(is.null(y.lb)) y.lb <- min(y)
  if(is.null(y.ub)) y.ub <- max(y)
  if(is.null(x.lb)) x.lb <- min(x)
  if(is.null(x.ub)) x.ub <- max(x)
  ## Test if non-unix type OS like Windows is running and set all cores to 1
  if(.Platform$OS.type=="unix") {
  } else {
    if(warnings) warning("non-unix type OS detected, setting all cores to 1 in bkcde() (no support for forking)")
    fitted.cores <- 1
    proper.cores <- 1
    optim.degree.cores <- 1
    optim.ksum.cores <- 1
    optim.nmulti.cores <- 1
  }
  if(is.null(x.eval) & is.null(y.eval)) {
    ## When infinite lower or upper values are provided, we provide the
    ## flexibility to use min/max for the evaluation data (default, *.trim=0),
    ## otherwise we allow the user to extend (*.trim < 0) or narrow (*.trim > 0)
    ## the range of the evaluation data. Admittedly clunky and could be refined
    ## further (i.e., with finite bounds could narrow the range of the
    ## evaluation data without injury, leave for a rainy day - the purpose is to
    ## assess behaviour of the integrated mean function).
    if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.grid)
    if(is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(y.lb,extendrange(y,f=y.erf)[2],length=n.grid)
    if(!is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(extendrange(y,f=y.erf)[1],y.ub,length=n.grid)
    if(!is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(extendrange(y,f=y.erf)[1],extendrange(y,f=y.erf)[2],length=n.grid)
    if(is.finite(x.lb) && is.finite(x.ub)) x.seq <- seq(x.lb,x.ub,length=n.grid)
    if(is.finite(x.lb) && !is.finite(x.ub)) x.seq <- seq(x.lb,extendrange(x,f=x.erf)[2],length=n.grid)
    if(!is.finite(x.lb) && is.finite(x.ub)) x.seq <- seq(extendrange(x,f=x.erf)[1],x.ub,length=n.grid)
    if(!is.finite(x.lb) && !is.finite(x.ub)) x.seq <- seq(extendrange(x,f=x.erf)[1],extendrange(x,f=x.erf)[2],length=n.grid)
    data.grid <- expand.grid(x.seq,y.seq)
    x.eval <- data.grid$Var1
    y.eval <- data.grid$Var2
    ## Set logical flag for evaluation data being auto-generated hence grid
    ## structure guaranteed
  }
  ## Set logical flag for evaluation data being auto-generated hence grid
  ## structure guaranteed (well, someone could break it I am sure but
  ## internally, if you let objects be created it will be fine)
  if((length(y.eval)/length(unique(y.eval))) == n.grid && (length(x.eval)/length(unique(x.eval))) == n.grid) {
    is.grid <- TRUE
  } else {
    is.grid <- FALSE
  }
  if(any(y<y.lb) | any(y>y.ub)) stop("y must lie in [y.lb,y.ub] in bkcde()")
  if(any(y.eval<y.lb) | any(y.eval>y.ub)) stop("y.eval must lie in [y.lb,y.ub] in bkcde()")
  if(any(x<x.lb) | any(x>x.ub)) stop("x must lie in [x.lb,x.ub] in bkcde()")
  if(any(x.eval<x.lb) | any(x.eval>x.ub)) stop("x.eval must lie in [x.lb,x.ub] in bkcde()")
  if(y.lb>=y.ub) stop("y.lb must be less than y.ub in bkcde()")
  if(x.lb>=x.ub) stop("x.lb must be less than x.ub in bkcde()")
  if(!is.logical(poly.raw)) stop("poly.raw must be logical in bkcde()")
  if(!is.logical(proper)) stop("proper must be logical in bkcde()")
  if(!is.logical(verbose)) stop("verbose must be logical in bkcde()")
  if(!is.logical(warnings)) stop("warnings must be logical in bkcde()")
  if(nmulti < 1) stop("nmulti must be at least 1 in bkcde()")
  if(n.integrate < 1) stop("n.integrate must be at least 1 in bkcde()")
  if(!is.numeric(n.binned) || length(n.binned) != 1 || n.binned < 2) stop("n.binned must be numeric and at least 2 in bkcde()")
  n.binned <- as.integer(n.binned)
  if(!is.null(h) & is.null(degree)) stop("must provide degree in bkcde() when h is not NULL")
  if(!is.null(h) & length(h) != 2) stop("h must be a vector of length 2 in bkcde()")
  if(!is.null(degree) && (degree < 0 | degree >= length(y))) stop("degree must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.min < 0 | degree.min >= length(y)) stop("degree.min must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.max < 0 | degree.max >= length(y)) stop("degree.max must lie in [0,1,...,",length(y)-1,"] (i.e., [0,1,dots, n-1]) in bkcde()")
  if(degree.min > degree.max) stop("degree.min must be <= degree.max in bkcde()")
  if(proper.cores < 1) stop("proper.cores must be at least 1 in bkcde()")
  if(fitted.cores < 1) stop("fitted.cores must be at least 1 in bkcde()")
  if(!is.null(optim.degree.cores) && optim.degree.cores < 1) stop("optim.degree.cores must be at least 1 in bkcde()")
  if(!is.null(optim.nmulti.cores) && optim.nmulti.cores < 1) stop("optim.nmulti.cores must be at least 1 in bkcde()")
  if(optim.sf.y.lb <= 0) stop("optim.sf.y.lb must be positive in bkcde()")
  if(optim.sf.x.lb <= 0) stop("optim.sf.x.lb must be positive in bkcde()")
  ## This chunk of code is to determine the number of cores to use for the
  ## optimization based on the number of models and the number of multistarts and
  ## the number of cores available via detectCores() in the parallel package. If
  ## neither optim.degree.cores nor optim.nmulti.cores is set, this tries to
  ## balance the load between the two, attempting to make full use of the
  ## available cores.
  nmodels <- degree.max-degree.min+1
  n.obs <- length(y)
  optim.cores <- match.arg(optim.cores)

  ## Manual (default) core allocation mirrors previous behavior
  if(optim.cores == "manual") {
    if(is.null(optim.ksum.cores)) {
      optim.ksum.cores <- 1L
    } else {
      optim.ksum.cores <- as.integer(optim.ksum.cores)
    }
    if(nmodels==1 & nmulti==1) {
      combn.out <- 1
    } else {
      combn.out <- combn(max(nmodels,nmulti),2)
      combn.out <- combn.out[,which(apply(combn.out,2,prod)<=detectCores()),drop=FALSE]
      combn.out <- combn.out[,ncol(combn.out)]
    }
    if(is.null(optim.degree.cores)) optim.degree.cores <- ifelse(nmodels >= nmulti,max(combn.out),min(combn.out))
    if(is.null(optim.nmulti.cores)) optim.nmulti.cores <- ifelse(nmodels < nmulti,max(combn.out),min(combn.out))
  }

  ## Auto mode: scale cores using detectCores(), sample size, and degree range
  if(optim.cores == "auto" && .Platform$OS.type=="unix") {
    C <- detectCores()
    ksum_max <- max(1L, as.integer(optim.ksum.auto.max))
    ksum_thresholds <- optim.ksum.auto.thresholds
    if(is.null(ksum_thresholds)) ksum_thresholds <- integer(0)
    ksum_thresholds <- sort(unique(as.integer(ksum_thresholds[ksum_thresholds > 0])))
    total_tasks <- nmodels * nmulti
    if(is.null(optim.ksum.cores)) {
      ## If the outer grid already saturates available cores and n is below the first threshold, avoid ksum overhead
      if(total_tasks <= C && length(ksum_thresholds) > 0 && n.obs < ksum_thresholds[1]) {
        optim.ksum.cores <- 1L
      } else {
        tier_levels <- pmin(ksum_max, 2L^(0:length(ksum_thresholds)))
        tier_idx <- findInterval(n.obs, ksum_thresholds) + 1L
        tier_idx <- min(tier_idx, length(tier_levels))
        target <- tier_levels[tier_idx]
        cap <- max(1L, floor(C/2))
        if(degree.max < 3 || C < 4) target <- 1L
        optim.ksum.cores <- max(1L, min(target, cap, ksum_max))
      }
    } else {
      optim.ksum.cores <- max(1L, as.integer(optim.ksum.cores))
    }
    limit <- max(1L, floor(C / optim.ksum.cores))

    ## Derive a balanced split within the available limit
    if(nmodels==1 & nmulti==1) {
      base_deg <- base_nmulti <- 1L
    } else {
      combn.out <- combn(max(nmodels,nmulti),2)
      combn.out <- combn.out[,which(apply(combn.out,2,prod)<=limit),drop=FALSE]
      if(ncol(combn.out)==0) combn.out <- matrix(c(1L,1L),nrow=2)
      choice <- combn.out[,ncol(combn.out)]
      base_deg <- if(nmodels >= nmulti) choice[1] else choice[2]
      base_nmulti <- if(nmodels >= nmulti) choice[2] else choice[1]
    }

    if(is.null(optim.degree.cores) || optim.degree.cores == detectCores()) optim.degree.cores <- base_deg
    if(is.null(optim.nmulti.cores) || optim.nmulti.cores == detectCores()) optim.nmulti.cores <- base_nmulti

    ## Ensure we do not exceed available cores after ksum split
    if(optim.degree.cores * optim.nmulti.cores > limit) {
      while(optim.degree.cores * optim.nmulti.cores > limit && optim.nmulti.cores > 1) {
        optim.nmulti.cores <- max(1L, floor(optim.nmulti.cores/2))
      }
      while(optim.degree.cores * optim.nmulti.cores > limit && optim.degree.cores > 1) {
        optim.degree.cores <- max(1L, floor(optim.degree.cores/2))
      }
      if(optim.degree.cores * optim.nmulti.cores > limit) {
        optim.nmulti.cores <- 1L
        optim.degree.cores <- 1L
      }
    }
  }
  if(optim.ksum.cores < 1) stop("optim.ksum.cores must be at least 1 in bkcde()")
  cv.penalty.method <- match.arg(cv.penalty.method)
  if(!is.logical(cv.binned)) stop("cv.binned must be logical in bkcde()")
  bwmethod <- match.arg(bwmethod)
  cv <- match.arg(cv)
  if(cv == "auto") cv <- ifelse(length(y) > cv.auto.threshold,"sub","full")
  if(cv.penalty.cutoff <= 0) stop("cv.penalty.cutoff must be positive in bkcde()")
  if(!is.null(h) & bwscaling) h <- h*EssDee(cbind(y,x))*length(y)^(-1/6)
  if(warnings && is.null(h) && (length(y) > 1e4 && cv == "full" && !cv.binned)) warning("large sample size for full sample cross-validation, consider cv='sub' in bkcde() [n = ",length(y),"]",immediate. = TRUE)
  
  ## Pre-compute y.seq for integration if required (computed once to avoid redundant computation)
  y.seq <- NULL
  if(proper || proper.cv) {
    y.seq <- seq(y.lb, y.ub, length.out = n.integrate)
  }
  ## This ensures Monte Carlo simulations are not disrupted
  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed <- TRUE
  } else {
    exists.seed <- FALSE
  }
  
  secs.start.total <- Sys.time()
  ## If no bandwidth is provided, then either likelihood or least-squares
  ## cross-validation is used to obtain the bandwidths and polynomial order (use
  ## optim.ksum.cores, optim.degree.cores, optim.nmulti.cores)
  if(is.null(h) & cv == "full") {
    if(progress) cat("\rNested optimization running (",degree.max-degree.min+1," models with ",nmulti," multistarts per model)...",sep="")
      optim.out <- bkcde.optim(x=x,
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
                               nmulti=nmulti,
                               n.integrate=n.integrate,
                               optim.degree.cores=optim.degree.cores,
                               optim.ksum.cores=optim.ksum.cores,
                               optim.nmulti.cores=optim.nmulti.cores,
                               optim.sf.y.lb=optim.sf.y.lb,
                               optim.sf.x.lb=optim.sf.x.lb,
                               cv.binned=cv.binned,
                               n.binned=n.binned,
                               poly.raw=poly.raw,
                               proper.cv=proper.cv,
                               y.seq=y.seq,
                               verbose=verbose,
                               seed=seed,
                               ...)
    h <- optim.out$par
    h.mat <- optim.out$par.mat
    h.x.init.mat <- optim.out$optim.x.init.mat
    h.y.init.mat <- optim.out$optim.y.init.mat
    degree <- optim.out$degree
    degree.mat <- optim.out$degree.mat
    value <- optim.out$value
    value.vec <- optim.out$value.vec
    value.mat <- optim.out$value.mat
    if(identical(value,-sqrt(.Machine$double.xmax))) {
      convergence <- 99
    } else {
      convergence <- optim.out$convergence
    }
    convergence.vec <- optim.out$convergence.vec
    convergence.mat <- optim.out$convergence.mat
    secs.optim <- optim.out$secs.optim
    secs.optim.mat <- optim.out$secs.optim.mat
    if(progress) cat("\rNested optimization complete (",degree.max-degree.min+1," models with ",nmulti," multistarts) in ",round(as.numeric(difftime(Sys.time(),secs.start.total,units="secs"))), " seconds\n",sep="")
    if(convergence != 0 && warnings) warning("optimization did not converge in bkcde(), consider increasing nmulti [degree = ",
                                 degree,
                                 ", h.y = ",
                                 round(h[1],5),
                                 ", h.x = ",
                                 round(h[2],5),
                                 "]",
                                 immediate. = TRUE)
  } else if(is.null(h) & cv == "sub") { 
    ## Code recursion in R is a thing of beauty sub.cv() calls bkcde()...
    if(progress) cat("\rSub-sample optimization (",degree.max-degree.min+1," models, ",nmulti," restarts/model, n.sub = ",n.sub,", resamples = ", resamples, ")\n",sep="")
    optimal <- sub.cv(x=x,
                      y=y,
                      x.eval=x.eval,
                      y.eval=y.eval,
                      n.sub=n.sub,
                      resamples=resamples,
                      nmulti=nmulti,
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
                      optim.degree.cores=optim.degree.cores,
                      optim.ksum.cores=optim.ksum.cores,
                      optim.nmulti.cores=optim.nmulti.cores,
                      optim.sf.y.lb=optim.sf.y.lb,
                      optim.sf.x.lb=optim.sf.x.lb,
                      cv.binned=cv.binned,
                      n.binned=n.binned,
                      seed=seed,
                      poly.raw=poly.raw,
                      progress=progress,
                      ...)
    h <- optimal$h.median
    degree <- optimal$degree
    h.mat <- NULL
    h.x.init.mat <- NULL
    h.y.init.mat <- NULL
    degree.mat <- NULL
    value <- NULL
    value.vec <- NULL
    value.mat <- NULL
    convergence <- NULL
    convergence.vec <- NULL
    convergence.mat <- NULL
    secs.optim <- NULL
    secs.optim.mat <- NULL
    if(progress) cat("\rSub-sample optimization (",degree.max-degree.min+1," models, ",nmulti," restarts/model, n.sub = ",n.sub,", resamples = ", resamples, ") in ",round(as.numeric(difftime(Sys.time(),secs.start.total,units="secs"))), " seconds\n",sep="")
  } else {
    h.mat <- NULL
    h.x.init.mat <- NULL
    h.y.init.mat <- NULL
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
  if(!cv.only) {
    if(progress) cat("\rFitting conditional density estimate...",sep="")
    secs.start.estimate <- Sys.time()
    if(!is.grid) {
      ## Compute the fitted conditional density estimate (use fitted.cores). Grid
      ## structure cannot be assumed, so instead we must go brute force as
      ## y.eval and/or x.eval are provided but could be equal to the sample x
      ## and y, etc. There is no savings possible here as even if the x values
      ## were identical the y values could all be different, so we need to go
      ## brute force for all (x.eval[i],y.eval[i]) pairs which could be
      ## (x[i],y[i]) if the user chooses. But if a grid structure is present and
      ## generated by this code, indeed we branch to the grid structure code
      ## below.
      if(degree == 0) {
        ## For degree 0 don't invoke the overhead associated with .lm.fit(), just
        ## compute the estimate \hat f(y|x) as efficiently as possible
        foo <- t(mcmapply(function(i){
          ## Here we exploit the fact that we can multiply columns of a matrix by
          ## the vector of kernel weights for x and avoid unnecessary computation
          pdf.kernel.bk.x<-pdf.kernel.bk(x.eval[i],x,h[2],x.lb,x.ub);
          colMeans(sweep(cbind(pdf.kernel.bk(y.eval[i],y,h[1],y.lb,y.ub),y,cdf.kernel.bk(y.eval[i],y,h[1],y.lb,y.ub)),1,pdf.kernel.bk.x,"*"))/NZD_pos(mean(pdf.kernel.bk.x))
        },seq_along(y.eval),mc.cores=fitted.cores))
        f.yx <- foo[,1]
        E.yx <- foo[,2]
        F.yx <- foo[,3]
        ## Derivatives extracted from degree 1 polynomial fit only (below), and
        ## only if poly.raw=TRUE
        f1.yx <- NULL
        E1.yx <- NULL
        F1.yx <- NULL
      } else {
        ## Choice of raw or orthogonal polynomials
        X.poly <- poly(x,raw=poly.raw,degree=degree)
        X <- cbind(1,X.poly)
        X.eval <- cbind(1,predict(X.poly,x.eval))
        ## For degree > 0 we use, e.g., lm(y~I(x^2)) and fitted values from the
        ## regression to estimate \hat f(y|x) rather than the intercept term from
        ## lm(y-I(x[i]-X)^2), which produce identical results for raw polynomials
        foo <- t(mcmapply(function(i){
          ## Here we exploit the fact that .lm.fit() can work on multivariate Y
          ## objects, so we don't need to recompute the X weight kernel sums and
          ## can avoid unnecessary computation
          w <- NZD_pos(sqrt(pdf.kernel.bk(x.eval[i],x,h[2],x.lb,x.ub)))
          beta.hat <- .lm.fit(X*w,cbind(pdf.kernel.bk(y.eval[i],y,h[1],y.lb,y.ub),y,cdf.kernel.bk(y.eval[i],y,h[1],y.lb,y.ub))*w)$coefficients
          ## The first three rows are the conditional density, the conditional
          ## mean, and conditional distribution, respectively, while the last
          ## three rows are the derivatives of the conditional density, mean, and
          ## distribution, respectively (computed only if degree=1)
          c(beta.hat[,1]%*%t(X.eval[i,,drop=FALSE]),
            beta.hat[,2]%*%t(X.eval[i,,drop=FALSE]),
            beta.hat[,3]%*%t(X.eval[i,,drop=FALSE]),
            beta.hat[2,1],
            beta.hat[2,2],
            beta.hat[2,3])
        },seq_along(y.eval),mc.cores=fitted.cores))
        f.yx <- foo[,1]
        E.yx <- foo[,2]
        F.yx <- foo[,3]
        ## Derivatives extracted from degree 1 polynomial fit only and only if
        ## poly.raw=TRUE
        if(degree==1 & poly.raw) {
          f1.yx <- foo[,4]
          E1.yx <- foo[,5]
          F1.yx <- foo[,6]
        } else {
          f1.yx <- NULL
          E1.yx <- NULL
          F1.yx <- NULL
        }
      }
    } else {
      ## Grid structure (hopefully) guaranteed since is.grid is TRUE. Compute
      ## the fitted conditional density estimate (use fitted.cores) on only the
      ## unique x values and use the multivariate Y capabilities to efficiently
      ## compute the coefficient matrix
      x.eval.unique <- unique(x.eval)
      y.eval.unique <- unique(y.eval)
      pdf.kernel.mat <- mapply(function(i) pdf.kernel.bk(y.eval.unique[i], y, h[1], y.lb, y.ub),seq_len(n.grid))
      cdf.kernel.mat <- mapply(function(i) cdf.kernel.bk(y.eval.unique[i], y, h[1], y.lb, y.ub),seq_len(n.grid))
      if(degree == 0) {
        ## For degree 0 don't invoke the overhead associated with .lm.fit(), just
        ## compute the estimate \hat f(y|x) as efficiently as possible.
        output <- mclapply(1:n.grid,function(i) {
          ## Here we exploit the fact that we can multiply columns of a matrix by
          ## the vector of kernel weights for x and avoid unnecessary computation
          pdf.kernel.bk.x<-pdf.kernel.bk(x.eval.unique[i],x,h[2],x.lb,x.ub);
          ## Coefficient matrix in order are column 1 E.y.x, 2 - n.grid+1 f.yx,
          ## n.grid+2 - 2*n.grid+1 F.yx
          foo <- colMeans(sweep(cbind(y,pdf.kernel.mat,cdf.kernel.mat),1,pdf.kernel.bk.x,"*"))/NZD_pos(mean(pdf.kernel.bk.x))
          return(list(E.yx=foo[1],
                      f.yx=foo[2:(n.grid+1)],
                      F.yx=foo[(n.grid+2):(2*n.grid+1)]))
        },mc.cores=fitted.cores)
        f.yx <- sapply(output, function(x) x$f.yx)
        f.yx <- f.yx[order(x.eval)]
        F.yx <- sapply(output, function(x) x$F.yx)
        F.yx <- F.yx[order(x.eval)]
        E.yx <- sapply(output, function(x) x$E.yx)
        E.yx <- E.yx[match(x.eval, x.eval.unique)]
        ## Derivatives extracted from degree 1 polynomial fit only (below), and
        ## only if poly.raw=TRUE
        f1.yx <- NULL
        E1.yx <- NULL
        F1.yx <- NULL
      } else {
        ## Choice of raw or orthogonal polynomials
        X.poly <- poly(x,raw=poly.raw,degree=degree)
        X <- cbind(1,X.poly)
        X.eval.unique <- cbind(1,predict(X.poly,x.eval.unique))
        ## For degree > 0 we use, e.g., lm(y~I(x^2)) and fitted values from the
        ## regression to estimate \hat f(y|x) rather than the intercept term from
        ## lm(y-I(x[i]-X)^2), which produce identical results for raw polynomials
        output <- mclapply(1:n.grid,function(i) {
          ## Here we exploit the fact that .lm.fit() can work on multivariate Y
          ## objects, so we don't need to recompute the X weight kernel sums and
          ## can avoid unnecessary computation
          w <- NZD_pos(sqrt(pdf.kernel.bk(x.eval.unique[i],x,h[2],x.lb,x.ub)))
          beta.hat <- .lm.fit(X*w,cbind(y,pdf.kernel.mat,cdf.kernel.mat)*w)$coefficients
          return(list(E.yx=as.numeric(X.eval.unique[i,,drop=FALSE]%*%beta.hat[,1,drop=FALSE]),
                      f.yx=as.numeric(X.eval.unique[i,,drop=FALSE]%*%beta.hat[,2:(n.grid+1),drop=FALSE]),
                      F.yx=as.numeric(X.eval.unique[i,,drop=FALSE]%*%beta.hat[,(n.grid+2):(2*n.grid+1),drop=FALSE]),
                      E1.yx=as.numeric(beta.hat[2,1]),
                      f1.yx=as.numeric(beta.hat[2,2:(n.grid+1)]),
                      F1.yx=as.numeric(beta.hat[2,(n.grid+2):(2*n.grid+1)])))          
        },mc.cores=fitted.cores)
        f.yx <- sapply(output, function(x) x$f.yx)
        f.yx <- f.yx[order(x.eval)]
        F.yx <- sapply(output, function(x) x$F.yx)
        F.yx <- F.yx[order(x.eval)]
        E.yx <- sapply(output, function(x) x$E.yx)
        E.yx <- E.yx[match(x.eval, x.eval.unique)]
        ## Derivatives extracted from degree 1 polynomial fit only and only if
        ## poly.raw=TRUE
        if(degree==1 & poly.raw) {
          f1.yx <- sapply(output, function(x) x$f1.yx)
          f1.yx <- f1.yx[order(x.eval)]
          F1.yx <- sapply(output, function(x) x$F1.yx)
          F1.yx <- F1.yx[order(x.eval)]
          E1.yx <- sapply(output, function(x) x$E1.yx)
          E1.yx <- E1.yx[match(x.eval, x.eval.unique)]
        } else {
          f1.yx <- NULL
          E1.yx <- NULL
          F1.yx <- NULL
        }
      }
      ## Kill objects that can be memory intensive and no longer needed
      if(degree==0) {
        rm(pdf.kernel.mat,cdf.kernel.mat)
      } else {
        rm(pdf.kernel.mat,cdf.kernel.mat,X,X.poly,X.eval.unique)
      }
    }
    if(progress) cat("\rFitted conditional density estimate complete in ",round(as.numeric(difftime(Sys.time(),secs.start.estimate,units="secs"))), " seconds\n",sep="")
    ## Ensure the estimate is proper (use proper.cores over unique(x.eval) which
    ## could be < # proper.cores allocated)
    if(proper) {
      ## When rendering proper, we cannot presume a grid structure, however, we
      ## definitely want to exploit computing the integral over the sequence of
      ## y values for only the unique X values, which could be even 1 if we are
      ## are doing a 2D plot which will have a vector of evaluation points for x
      ## equal in length to y.eval but they are all identical.
      f.yx.unadjusted <- f.yx
      E.yx.unadjusted <- E.yx
      F.yx.unadjusted <- F.yx
      if(progress) cat("\rComputing integrals to ensure estimate is proper...\n",sep="")
      ## Create a sequence of values along an appropriate grid to compute the integral.
      if(is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(y.lb,y.ub,length=n.integrate)
      if(is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(y.lb,extendrange(y,f=integrate.erf)[2],length=n.integrate)
      if(!is.finite(y.lb) && is.finite(y.ub)) y.seq <- seq(extendrange(y,f=integrate.erf)[1],y.ub,length=n.integrate)
      if(!is.finite(y.lb) && !is.finite(y.ub)) y.seq <- seq(extendrange(y,f=integrate.erf)[1],extendrange(y,f=integrate.erf)[2],length=n.integrate)
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
      int.f.seq.post.mat <- matrix(NA,nrow=length(x.eval.unique),ncol=n.integrate)
      E.yx <- numeric()
      ## We test for only 1 unique value of x.eval to avoid parallel processing in
      ## the outer mcmapply call and invoke fitting the mcmapply sequence of
      ## f(y|x) values with proper.cores
      Y.seq.mat <- mapply(function(i) pdf.kernel.bk(y.seq[i], y, h[1], y.lb, y.ub),seq_len(n.integrate))
      if(degree > 0) {
        ## Not function of j, so we can compute this once
        X.poly <- poly(x,raw=poly.raw,degree=degree)
        X <- cbind(1,X.poly)
        X.eval.unique <- cbind(1,predict(X.poly,x.eval.unique))
      }
      ## Precompute kernels between unique x.eval values and sample x to avoid
      ## re-evaluating pdf.kernel.bk() repeatedly (preserves numerical identity)
      pdf.kernel.bk.x.mat <- sapply(x.eval.unique, function(xx) pdf.kernel.bk(xx, x, h[2], x.lb, x.ub))

      proper.out <- mclapply.progress(seq_along(x.eval.unique),function(j) {
        if(degree == 0) {
          pdf.kernel.bk.x <- pdf.kernel.bk.x.mat[, j]
          f.seq <- colMeans(Y.seq.mat * pdf.kernel.bk.x / NZD_pos(mean(pdf.kernel.bk.x)))
        } else {
          w <- NZD_pos(sqrt(pdf.kernel.bk.x.mat[, j]))
          beta.hat <- .lm.fit(X * w, Y.seq.mat * w)$coefficients
          f.seq <- as.numeric(X.eval.unique[j, , drop = FALSE] %*% beta.hat)
        }
        ## Compute integral of f.seq including any possible negative values
        int.f.seq.pre.neg[j] <- integrate.trapezoidal(y.seq, f.seq)[n.integrate]
        ## Set any possible negative f.seq values to 0
        f.seq[f.seq < 0] <- 0
        ## Compute integral of f.seq after setting any possible negative values to 0
        int.f.seq[j] <- integrate.trapezoidal(y.seq, f.seq)[n.integrate]
        ## Compute integral of f.seq after setting any possible negative values to 0
        ## and correcting to ensure final estimate integrates to 1
        int.f.seq.post.mat[j, ] <- integrate.trapezoidal(y.seq, f.seq / int.f.seq[j])
        int.f.seq.post[j] <- int.f.seq.post.mat[j, n.integrate]
        E.yx[j] <- integrate.trapezoidal(y.seq, y.seq * f.seq / int.f.seq[j])[n.integrate]
        return(list(int.f.seq.pre.neg = int.f.seq.pre.neg[j],
                    int.f.seq = int.f.seq[j],
                    int.f.seq.post = int.f.seq.post[j],
                    int.f.seq.post.mat = int.f.seq.post.mat[j, ],
                    E.yx = E.yx[j]))
      },mc.cores = ifelse(length(x.eval.unique)>1,proper.cores,1),progress=progress)
      ## Now gather the results, correct for negative entries then divide elements
      ## of f.xy by the corresponding integral (one for each x.eval.unique) to
      ## ensure the estimate is proper
      if(verbose & any(f.yx < 0)) if(warnings) warning("negative density estimate reset to 0 via option proper=TRUE in bkcde() [degree = ",
                                          degree,
                                          ", j = ",
                                          length(f.yx[f.yx < 0]),
                                          " element(s), min = ",
                                          formatC(min(f.yx[f.yx < 0]),format="f",digits=6), 
                                          ", h.y = ",
                                          round(h[1],5),
                                          ", h.x = ",
                                          round(h[2],5),
                                          "]",
                                          immediate. = TRUE)
      f.yx[f.yx < 0] <- 0
      int.f.seq.pre.neg <- sapply(proper.out, function(x) x$int.f.seq.pre.neg)
      int.f.seq <- sapply(proper.out, function(x) x$int.f.seq)
      int.f.seq.post <- sapply(proper.out, function(x) x$int.f.seq.post)
      int.f.seq.post.mat <- t(sapply(proper.out, function(x) x$int.f.seq.post.mat))
      int.f.seq.post.mat <- int.f.seq.post.mat[match(x.eval, x.eval.unique),]
      F.yx <- sapply(seq_along(x.eval), function(j) {max(int.f.seq.post.mat[j, y.seq <= y.eval[j]])})
      F.yx <- (F.yx-min(F.yx))/(max(F.yx)-min(F.yx))
      E.yx <- sapply(proper.out, function(x) x$E.yx)
      E.yx <- E.yx[match(x.eval, x.eval.unique)]
      f.yx <- f.yx/int.f.seq[match(x.eval, x.eval.unique)]
      ## As a summary measure report the mean of the integrals
      int.f.seq.pre.neg <- mean(int.f.seq.pre.neg)
      int.f.seq <- mean(int.f.seq)
      int.f.seq.post <- mean(int.f.seq.post)
      ## Kill objects that can be memory intensive and no longer needed
      if(degree==0) {
        rm(Y.seq.mat)
      } else {
        rm(Y.seq.mat,X,X.poly,X.eval.unique)
      }
      if(progress) cat("\rComputed integrals to ensure estimate is proper complete in ",round(as.numeric(difftime(Sys.time(),secs.start.estimate,units="secs"))), " seconds\n",sep="")
    } else {
      int.f.seq.pre.neg <- NULL
      int.f.seq <- NULL
      int.f.seq.post <- NULL
      f.yx.unadjusted <- NULL
      E.yx.unadjusted <- NULL
      F.yx.unadjusted <- NULL
      if(any(f.yx < 0)) if(warnings) warning("negative density estimate encountered, consider option proper=TRUE in bkcde() [degree = ",
                                degree,
                                ", ", 
                                length(f.yx[f.yx < 0]),
                                " element(s), min = ",
                                formatC(min(f.yx),format="f",digits=6), 
                                ", h.y = ",
                                round(h[1],5),
                                ", h.x = ",
                                round(h[2],5),
                                "]",
                                immediate. = TRUE)
    }
  } else {
    f.yx <- NULL
    f.yx.unadjusted <- NULL
    F.yx <- NULL
    F.yx.unadjusted <- NULL
    E.yx <- NULL
    E.yx.unadjusted <- NULL
    f1.yx <- NULL
    E1.yx <- NULL
    F1.yx <- NULL
    int.f.seq.pre.neg <- NULL
    int.f.seq <- NULL
    int.f.seq.post <- NULL
  }
  return.list <- list(convergence.mat=convergence.mat,
                      convergence.vec=convergence.vec,
                      convergence=convergence,
                      optim.cores=optim.cores,
                      cv=cv,
                      cv.binned=cv.binned,
                      cv.only=cv.only,
                      bwmethod=bwmethod,
                      degree.mat=degree.mat,
                      degree.max=degree.max,
                      degree.min=degree.min,
                      degree=degree,
                      F=F.yx,
                      F1=F1.yx,
                      F.unadjusted=F.yx.unadjusted,
                      fitted.cores=fitted.cores,
                      f.yx.integral.post=int.f.seq.post,
                      f.yx.integral.pre.neg=int.f.seq.pre.neg,
                      f.yx.integral=int.f.seq,
                      f=f.yx,
                      f1=f1.yx,
                      f.unadjusted=f.yx.unadjusted,
                      g=E.yx,
                      g1=E1.yx,
                      g.unadjusted=E.yx.unadjusted,
                      h.mat=h.mat,
                      h=h,
                      h.x.init.mat=h.x.init.mat,
                      h.y.init.mat=h.y.init.mat,
                      h.sf=h/(EssDee(cbind(y,x))*length(y)^(-1/6)),
                      n.binned=n.binned,
                      n.grid=n.grid,
                      n.sub=n.sub,
                      n.integrate=n.integrate,
                      optim.degree.cores=optim.degree.cores,
                      optim.ksum.cores=optim.ksum.cores,
                      optim.nmulti.cores=optim.nmulti.cores,
                      optimize=optimize,
                      proper.cores=proper.cores,
                      proper=proper,
                      resamples=resamples,
                      secs.elapsed=as.numeric(difftime(Sys.time(),secs.start.total,units="secs")),
                      secs.estimate=ifelse(cv.only,NA,as.numeric(difftime(Sys.time(),secs.start.estimate,units="secs"))),
                      secs.optim.mat=secs.optim.mat,
                      value.mat=value.mat,
                      value.vec=value.vec,
                      value=value,
                      warnings=warnings,
                      x.eval=x.eval,
                      x.lb=x.lb,
                      x.ub=x.ub,
                      x=x,
                      y.eval=y.eval,
                      y.lb=y.lb,
                      y.ub=y.ub,
                      y=y)
  class(return.list) <- "bkcde"
  
  ## Restore seed before returning (if it was saved at the start)
  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
  
  return(return.list)
}


## persp() does not allow ylim and zlim to be set to NULL, so this function
## allows for these options to be passable thereby avoiding unnecessary
## duplication of code in plot.bkcde()




## Allow for progress to be displayed while running in parallel, use pbmcapply()
## instead of mcmapply() and pbmclapply() instead of mclapply() to display a
## progress bar. The progress bar is only displayed if progress=TRUE in the
## function calls.




## plot.bkcde() is used to plot the results of the boundary kernel estimate of
## f(y|x) along with bootstrap confidence intervals generated as either
## pointwise, Bonferroni, or simultaneous intervals. Both 2D and 3D plots can be
## produced (persp=TRUE or FALSE). A simple bootstrap mean correction is
## applied (there are likely better ways of doing this outside of a bootstrap
## procedure, leave this for future investigation). A handful of options are
## available, including returning the confidence intervals (pointwise,
## Bonferroni and simultaneous) and estimates (plot.behavior=...). Note that a
## large number of bootstrap replications should be used, particularly if
## Bonferroni corrections are requested, and in such cases the number of
## bootstrap replications should increase with the grid size (i.e., the number
## of evaluation points). The function uses the SCSrank() function from the
## MCPAN package to compute simultaneous confidence intervals.

plot.bkcde <- function(x,
                       B = 3999,
                       alpha = 0.05,
                       ci = FALSE,
                       ci.bias.correct = TRUE,
                       ci.cores = NULL,
                       ci.method = c("Pointwise","Bonferroni","Simultaneous","all"),
                       ci.preplot = TRUE,
                       fitted.cores = NULL,
                       n.grid = NULL,
                       persp = TRUE,
                       phi = NULL,
                       plot.behavior = c("plot","plot-data","data"),
                       plot.unadjusted = FALSE,
                       progress = FALSE,
                       proper = NULL,
                       proper.cores = NULL,
                       sub = NULL,
                       theta = NULL,
                       type = NULL,
                       x.eval = NULL,
                       xlab = NULL,
                       y.eval = NULL,
                       ylab = NULL,
                       ylim = NULL,
                       zlab = NULL,
                       zlim = NULL,
                       ...) {
  if(!inherits(x,"bkcde")) stop("x must be of class bkcde in plot.bkcde()")
  if(x$cv.only) stop("x must be the output of bkcde() called with cv.only=FALSE in plot.bkcde()")
  if(!is.logical(ci)) stop("ci must be logical in plot.bkcde()")
  if(!is.logical(ci.bias.correct)) stop("ci.bias.correct must be logical in plot.bkcde()")
  if(!is.logical(ci.preplot)) stop("ci.preplot must be logical in plot.bkcde()")
  if(!is.logical(plot.unadjusted)) stop("plot.unadjusted must be logical in plot.bkcde()")
  if(!is.logical(progress)) stop("progress must be logical in plot.bkcde()")
  if(!is.logical(ci.preplot)) stop("ci.preplot must be logical in plot.bkcde()")
  ## Note that the Bonferroni method places a restriction on the smallest number
  ## of bootstrap replications such that (alpha/(2*n.grid))*(B+1) or
  ## (alpha/(2*n.grid^2))*(B+1) is a positive integer - this is on the
  ## user to ensure. The simultaneous method will almost certainly require fewer
  ## bootstrap replications than the Bonferroni method. The pointwise method
  ## requires the fewest, alpha*(B+1) must be a positive integer. So, for the
  ## Bonferroni bounds the smallest number of bootstrap replications using the
  ## default (alpha=0.05, n.grid=100 or alpha=0.05, n.grid=10 is
  ## 3999. For non-default values the user should adjust B accordingly.
  ci.method <- match.arg(ci.method)
  plot.behavior <- match.arg(plot.behavior)
  if(alpha <= 0 | alpha >= 1) stop("alpha must lie in (0,1) in plot.bkcde()")
  if(B < 1) stop("B must be at least 1 in plot.bkcde()")
  if(.Platform$OS.type=="unix") {
  } else {
    warning("non-unix type OS detected, setting all cores to 1 in plot.bkcde() (no support for forking)")
    ci.cores <- 1
    fitted.cores <-1
    proper.cores <- 1
  }
  if(!is.null(ci.cores) && ci.cores < 1) stop("ci.cores must be at least 1 in plot.bkcde()")
  f.ci.pw.lb <- f.ci.pw.ub <- f.ci.bf.lb <- f.ci.bf.ub <- f.ci.sim.lb <- f.ci.sim.ub <- bias.vec <- NULL
  g.ci.pw.lb <- g.ci.pw.ub <- g.ci.bf.lb <- g.ci.bf.ub <- g.ci.sim.lb <- g.ci.sim.ub <- NULL
  F.ci.pw.lb <- F.ci.pw.ub <- F.ci.bf.lb <- F.ci.bf.ub <- F.ci.sim.lb <- F.ci.sim.ub <- NULL
  f1.ci.pw.lb <- f1.ci.pw.ub <- f1.ci.bf.lb <- f1.ci.bf.ub <- f1.ci.sim.lb <- f1.ci.sim.ub <- NULL
  g1.ci.pw.lb <- g1.ci.pw.ub <- g1.ci.bf.lb <- g1.ci.bf.ub <- g1.ci.sim.lb <- g1.ci.sim.ub <- NULL
  F1.ci.pw.lb <- F1.ci.pw.ub <- F1.ci.bf.lb <- F1.ci.bf.ub <- F1.ci.sim.lb <- F1.ci.sim.ub <- NULL
  if(!is.null(proper) & !is.logical(proper)) stop("proper must be logical in plot.bkcde()")
  if(persp & !is.null(x.eval) & !is.null(y.eval) & length(x.eval) != length(y.eval)) stop("length of x.eval must be equal to length of y.eval in plot.bkcde()")
  if(!is.null(n.grid) && n.grid < 2) stop("n.grid must be at least 2 in plot.bkcde()")
  if(persp & is.null(n.grid)) {
    n.grid <- x$n.grid
  } else if(is.null(n.grid)) {
    ## Default for 2D grid if not specified is x$n.grid*4 (100)
    n.grid <- x$n.grid*4
  }
  if(ci.preplot==FALSE & ci==FALSE & ci.method != "data") stop("ci.preplot must be TRUE when ci is TRUE and ci.method is not 'data' in plot.bkcde()")
  if(is.null(proper)) {
    plot.proper <- NULL
    proper <- x$proper
  } else {
    plot.proper <- proper
  }
  if(!proper & plot.unadjusted) stop("plot.unadjusted=TRUE requires proper=TRUE in bkcde() or in plot.bkcde() (i.e., plot())")
  if(!persp & is.null(x.eval)) x.eval <- median(x$x.eval)
  if(!persp & !is.null(x.eval) & length(x.eval) >1) stop("x.eval must be a scalar in plot.bkcde() when persp=FALSE")
  ## proper.cores and fitted.cores are used in the predict() function and only
  ## in the bootstrap if ci=TRUE and ci.cores>1
  if(is.null(proper.cores)) proper.cores <- detectCores()
  if(is.null(fitted.cores)) fitted.cores <- detectCores()
  secs.start <- Sys.time()
  ## For the user, whether ci=TRUE or not get the estimate plotted asap
  ## otherwise they are faced with a blank screen
  if(persp) {
    ## Plot 3D
    if(is.null(x.eval)) {
      x.eval <- seq(min(x$x.eval),max(x$x.eval),length=n.grid)
      y.eval <- seq(min(x$y.eval),max(x$y.eval),length=n.grid)
    } else {
      x.eval <- x.eval
      y.eval <- y.eval
      n.grid <- length(x.eval)
    }
    if(length(unique(x.eval))==1) stop("only one unique x.eval value, cannot call persp() in plot.bkcde() (perhaps call bkcde() with non-unique x.eval OR call plot() with persp=FALSE OR call plot() and provide x.eval?)")
    if(length(unique(y.eval))==1) stop("only one unique y.eval value, cannot call persp() in plot.bkcde() (perhaps call bkcde() with non-unique y.eval OR call plot() and provide y.eval?)")
    data.grid <- expand.grid(x.eval,y.eval)
    x.plot.eval <- data.grid$Var1
    y.plot.eval <- data.grid$Var2
    rm(data.grid)
    if((n.grid != x$n.grid) || (!is.null(plot.proper) && plot.proper!=x$proper) || !identical(x.plot.eval,x$x.eval) || !identical(y.plot.eval,x$y.eval)) {
      f.yx.plot <- bkcde(h=x$h,
                         x=x$x,
                         y=x$y,
                         x.eval=x.plot.eval,
                         y.eval=y.plot.eval,
                         y.lb=x$y.lb,
                         y.ub=x$y.ub,
                         x.lb=x$x.lb,
                         x.ub=x$x.ub,
                         n.grid=n.grid,
                         proper=proper,
                         degree=x$degree,
                         fitted.cores=fitted.cores,
                         proper.cores=proper.cores,
                         progress=progress,
                         ...)
      f.fitted <- f.yx.plot$f
      f.fitted.unadjusted <- f.yx.plot$f.unadjusted
      g.fitted <- f.yx.plot$g
      F.fitted <- f.yx.plot$F
      f1.fitted <- f.yx.plot$f1
      g1.fitted <- f.yx.plot$g1
      F1.fitted <- f.yx.plot$F1
    } else {
      f.fitted <- x$f
      f.fitted.unadjusted <- x$f.unadjusted
      g.fitted <- x$g
      F.fitted <- x$F
      f1.fitted <- x$f1
      g1.fitted <- x$g1
      F1.fitted <- x$F1
    }
    predict.mat <- matrix(f.fitted,n.grid,n.grid)
    if(plot.unadjusted) predict.mat.unadjusted <- matrix(f.fitted.unadjusted,n.grid,n.grid)
    if(is.null(theta)) theta <- 120
    if(is.null(phi)) phi <- 45
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- "y"
    if(is.null(zlab)) zlab <- "f(y|x)"
    if(ci.preplot & plot.behavior != "data") {
      if(!plot.unadjusted) {
        persp.lim(x=x.eval,y=y.eval,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,shade=.25,...)    
      } else {
        if(is.null(zlim)) zlim <- range(predict.mat,predict.mat.unadjusted)
        persp.lim(x=x.eval,y=y.eval,z=predict.mat.unadjusted,xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="red",col=NA,lty=2,lwd=2,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",border="black",col=NA,ylim=ylim,zlim=zlim,...)
        legend("topleft",legend=c("Adjusted","Unadjusted"),lty=c(1,2),col=c(1,2),lwd=c(1,2),bty="n")
      }
    }
  } else if(!persp) {
    ## Plot 2D
    predict.mat <- NULL
    if(is.null(y.eval)) {
      y.plot.eval <- y.eval <- seq(min(x$y.eval),max(x$y.eval),length=n.grid)
    } else {
      y.plot.eval <- y.eval
      n.grid <- length(y.eval)
    }
    x.plot.eval <- x.eval <- rep(x.eval,length(y.plot.eval))
    if((n.grid != x$n.grid) || (!is.null(plot.proper) && plot.proper!=x$proper) || length(unique(x$x.eval)) > 1) {
      f.yx.plot <- bkcde(h=x$h,
                         x=x$x,
                         y=x$y,
                         x.eval=x.plot.eval,
                         y.eval=y.plot.eval,
                         y.lb=x$y.lb,
                         y.ub=x$y.ub,
                         x.lb=x$x.lb,
                         x.ub=x$x.ub,
                         n.grid=n.grid,
                         proper=proper,
                         degree=x$degree,
                         fitted.cores=fitted.cores,
                         proper.cores=proper.cores,
                         progress=progress,
                         ...)
      f.fitted <- f.yx.plot$f
      f.fitted.unadjusted <- f.yx.plot$f.unadjusted
      g.fitted <- f.yx.plot$g
      F.fitted <- f.yx.plot$F
      f1.fitted <- f.yx.plot$f1
      g1.fitted <- f.yx.plot$g1
      F1.fitted <- f.yx.plot$F1
    } else {
      f.fitted <- x$f
      f.fitted.unadjusted <- x$f.unadjusted
      g.fitted <- x$g
      F.fitted <- x$F
      f1.fitted <- x$f1
      g1.fitted <- x$g1
      F1.fitted <- x$F1
    }
    if(is.null(sub)) sub <- paste("(degree = ",x$degree,", h.y = ",round(x$h[1],3), ", h.x = ",round(x$h[2],3),", n = ",length(x$y),")",sep="")
    if(is.null(ylab)) ylab <- "f(y|x)"
    if(is.null(xlab)) xlab <- paste("y|x=",round(x.eval[1],digits=2),sep="")
    if(is.null(type)) type <- "l"
    if(ci.preplot & plot.behavior != "data") {
      if(!plot.unadjusted) {
        plot(y.plot.eval[order(y.plot.eval)],f.fitted[order(y.plot.eval)],
             sub=sub,
             ylim=ylim,
             ylab=ylab,
             xlab=xlab,
             type=type,
             panel.first=grid(lty=1),
             ...)
      } else {
        if(is.null(ylim)) ylim <- range(f.fitted,f.fitted.unadjusted)
        plot(y.plot.eval[order(y.plot.eval)],f.fitted[order(y.plot.eval)],
             sub=sub,
             ylim=ylim,
             ylab=ylab,
             xlab=xlab,
             type=type,
             panel.first=grid(lty=1),
             ...)        
        lines(y.plot.eval[order(y.plot.eval)],f.fitted.unadjusted[order(y.plot.eval)],lty=2,col=2,lwd=2)
        legend("topleft",legend=c("Adjusted","Unadjusted"),lty=c(1,2),col=c(1,2),lwd=c(1,2),bty="n")
      }
    }
  }
  if(ci) {
    if(progress) cat("\rComputing bootstrap confidence intervals...\n",sep="") 
    ## All processing goes into computing the matrix of bootstrap estimates, so
    ## once this is done it makes sense to then generate all three types of
    ## confidence intervals
    if(is.null(ci.cores)) ci.cores <- detectCores()
    boot.return <- mclapply.progress(seq_len(B),function(b){
      ii <- sample(seq_along(x$y),replace=TRUE)
      bkcde.out <- bkcde(h=x$h,
                         x=x$x[ii],
                         y=x$y[ii],
                         x.eval=x.plot.eval,
                         y.eval=y.plot.eval,
                         y.lb=x$y.lb,
                         y.ub=x$y.ub,
                         x.lb=x$x.lb,
                         x.ub=x$x.ub,
                         n.grid=n.grid,
                         fitted.cores=ifelse(ci.cores>1,1,fitted.cores),
                         proper=proper,
                         proper.cores=ifelse(ci.cores>1,1,proper.cores),
                         degree=x$degree)
      return(list(f=bkcde.out$f,
                  g=bkcde.out$g,
                  F=bkcde.out$F,
                  f1=bkcde.out$f1,
                  g1=bkcde.out$g1,
                  F1=bkcde.out$F1))
    },mc.cores=ci.cores,progress=progress)
    f.boot.mat <- t(sapply(boot.return, function(x) x$f))
    F.boot.mat <- t(sapply(boot.return, function(x) x$F))
    g.boot.mat <- t(sapply(boot.return, function(x) x$g))
    if(!is.null(x$f1) & !is.null(x$F1) & !is.null(x$g1)) {
      f1.boot.mat <- t(sapply(boot.return, function(x) x$f1))
      F1.boot.mat <- t(sapply(boot.return, function(x) x$F1))
      g1.boot.mat <- t(sapply(boot.return, function(x) x$g1))
    } else {
      f1.boot.mat <- NULL
      F1.boot.mat <- NULL
      g1.boot.mat <- NULL
    }
    if(ci.bias.correct) {
      bias.vec <- colMeans(f.boot.mat) - f.fitted
      f.boot.mat <- sweep(f.boot.mat,2,bias.vec,"-")
      bias.vec <- colMeans(F.boot.mat) - F.fitted
      F.boot.mat <- sweep(F.boot.mat,2,bias.vec,"-")
      bias.vec <- colMeans(g.boot.mat) - g.fitted
      g.boot.mat <- sweep(g.boot.mat,2,bias.vec,"-")
      if(!is.null(f1.fitted) & !is.null(F1.fitted) & !is.null(g1.fitted)) {
        bias.vec <- colMeans(f1.boot.mat) - f1.fitted
        f1.boot.mat <- sweep(f1.boot.mat,2,bias.vec,"-")
        bias.vec <- colMeans(F1.boot.mat) - F1.fitted
        F1.boot.mat <- sweep(F1.boot.mat,2,bias.vec,"-")
        bias.vec <- colMeans(g1.boot.mat) - g1.fitted
        g1.boot.mat <- sweep(g1.boot.mat,2,bias.vec,"-")
      }
    }
    f.ci.pw.lb <- apply(f.boot.mat, 2, quantile, probs = alpha / 2)
    f.ci.pw.ub <- apply(f.boot.mat, 2, quantile, probs = 1 - alpha / 2)
    f.ci.bf.lb <- apply(f.boot.mat, 2, quantile, probs = alpha / (2 * length(y.plot.eval)))
    f.ci.bf.ub <- apply(f.boot.mat, 2, quantile, probs = 1 - alpha / (2 * length(y.plot.eval)))
    f.ci.SCS <- SCSrank(f.boot.mat, conf.level=1-alpha)$conf.int
    f.ci.sim.lb <- f.ci.SCS[,1]
    f.ci.sim.ub <- f.ci.SCS[,2]
    F.ci.pw.lb <- apply(F.boot.mat, 2, quantile, probs = alpha / 2)
    F.ci.pw.ub <- apply(F.boot.mat, 2, quantile, probs = 1 - alpha / 2)
    F.ci.bf.lb <- apply(F.boot.mat, 2, quantile, probs = alpha / (2 * length(y.plot.eval)))
    F.ci.bf.ub <- apply(F.boot.mat, 2, quantile, probs = 1 - alpha / (2 * length(y.plot.eval)))
    F.ci.SCS <- SCSrank(F.boot.mat, conf.level=1-alpha)$conf.int
    F.ci.sim.lb <- F.ci.SCS[,1]
    F.ci.sim.ub <- F.ci.SCS[,2]
    if(persp) {
      g.ci.pw.lb <- apply(g.boot.mat, 2, quantile, probs = alpha / 2)
      g.ci.pw.ub <- apply(g.boot.mat, 2, quantile, probs = 1 - alpha / 2)
      g.ci.bf.lb <- apply(g.boot.mat, 2, quantile, probs = alpha / (2 * length(y.plot.eval)))
      g.ci.bf.ub <- apply(g.boot.mat, 2, quantile, probs = 1 - alpha / (2 * length(y.plot.eval)))
      g.ci.SCS <- SCSrank(g.boot.mat, conf.level=1-alpha)$conf.int
      g.ci.sim.lb <- g.ci.SCS[,1]
      g.ci.sim.ub <- g.ci.SCS[,2]
    }
    if(!is.null(x$f1) & !is.null(x$F1) & !is.null(x$g1)){
      f1.ci.pw.lb <- apply(f1.boot.mat, 2, quantile, probs = alpha / 2)
      f1.ci.pw.ub <- apply(f1.boot.mat, 2, quantile, probs = 1 - alpha / 2)
      f1.ci.bf.lb <- apply(f1.boot.mat, 2, quantile, probs = alpha / (2 * length(y.plot.eval)))
      f1.ci.bf.ub <- apply(f1.boot.mat, 2, quantile, probs = 1 - alpha / (2 * length(y.plot.eval)))
      f1.ci.SCS <- SCSrank(f1.boot.mat, conf.level=1-alpha)$conf.int
      f1.ci.sim.lb <- f1.ci.SCS[,1]
      f1.ci.sim.ub <- f1.ci.SCS[,2]
      F1.ci.pw.lb <- apply(F1.boot.mat, 2, quantile, probs = alpha / 2)
      F1.ci.pw.ub <- apply(F1.boot.mat, 2, quantile, probs = 1 - alpha / 2)
      F1.ci.bf.lb <- apply(F1.boot.mat, 2, quantile, probs = alpha / (2 * length(y.plot.eval)))
      F1.ci.bf.ub <- apply(F1.boot.mat, 2, quantile, probs = 1 - alpha / (2 * length(y.plot.eval)))
      F1.ci.SCS <- SCSrank(F1.boot.mat, conf.level=1-alpha)$conf.int
      F1.ci.sim.lb <- F1.ci.SCS[,1]
      F1.ci.sim.ub <- F1.ci.SCS[,2]
      if(persp) {
        g1.ci.pw.lb <- apply(g1.boot.mat, 2, quantile, probs = alpha / 2)
        g1.ci.pw.ub <- apply(g1.boot.mat, 2, quantile, probs = 1 - alpha / 2)
        g1.ci.bf.lb <- apply(g1.boot.mat, 2, quantile, probs = alpha / (2 * length(y.plot.eval)))
        g1.ci.bf.ub <- apply(g1.boot.mat, 2, quantile, probs = 1 - alpha / (2 * length(y.plot.eval)))
      }
    }
    if(progress) cat("\rComputed bootstrap confidence intervals in ",round(as.numeric(difftime(Sys.time(),secs.start,units="secs"))), " seconds\n",sep="")
  } else {
    if(is.null(ylim)) ylim <- NULL
  }
  if(plot.behavior != "data") {
    if(persp & ci) {
      ## Plot 3D again with zlim set for the confidence intervals (not checking for if(is.null(zlim)) yet)
      if(ci.method == "Pointwise") {
        if(is.null(zlim)) zlim <-  range(c(f.fitted,f.ci.pw.lb,f.ci.pw.ub))
      } else if(ci.method == "Bonferroni") {
        if(is.null(zlim)) zlim <-  range(c(f.fitted,f.ci.bf.lb,f.ci.bf.ub))
      } else if(ci.method == "Simultaneous") {
        if(is.null(zlim)) zlim <-  range(c(f.fitted,f.ci.sim.lb,f.ci.sim.ub))
      } else {
        if(is.null(zlim)) zlim <-  range(c(f.fitted,f.ci.pw.lb,f.ci.pw.ub,f.ci.bf.lb,f.ci.bf.ub,f.ci.sim.lb,f.ci.sim.ub))
      }
      ## Unlike plot() persp() does accept a null ylim argument so we need to check...
      if(ci.method == "Pointwise") {
        ## First lower, then plot, then upper (surfaces)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.pw.lb,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.pw.ub,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci.method == "Bonferroni") {
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.bf.lb,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.bf.ub,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci.method == "Simultaneous") {
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.sim.lb,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.sim.ub,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci & ci.method == "all") {
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.pw.lb,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.pw.ub,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=2,ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.sim.lb,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=3,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=predict.mat,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,ticktype="detailed",ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.sim.ub,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=3,ylim=ylim,zlim=zlim,...)
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.bf.lb,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=4,ylim=ylim,zlim=zlim,...) 
        par(new = TRUE)
        persp.lim(x=x.eval,y=y.eval,z=matrix(f.ci.bf.ub,n.grid,n.grid),xlab="",ylab="",zlab="",theta=theta,phi=phi,ticktype="detailed",border="grey",col=NA,lty=4,ylim=ylim,zlim=zlim,...)
        legend("topright",
               legend=c("Estimated f(y|x)",
                        paste(100*(1-alpha),"% Pointwise CIs",sep=""),
                        paste(100*(1-alpha),"% Simultaneous CIs",sep=""),
                        paste(100*(1-alpha),"% Bonferroni CIs",sep="")
               ),lty=1:4,col=c(1,rep("grey",3)),bty="n")
      }
    } else if(ci) {
      ## Plot 2D again with ylim set for the confidence intervals
      if(ci.method == "Pointwise") {
        if(is.null(ylim)) ylim <-  range(c(f.fitted,f.ci.pw.lb,f.ci.pw.ub))
      } else if(ci.method == "Bonferroni") {
        if(is.null(ylim)) ylim <-  range(c(f.fitted,f.ci.bf.lb,f.ci.bf.ub))
      } else if(ci.method == "Simultaneous") {
        if(is.null(ylim)) ylim <-  range(c(f.fitted,f.ci.sim.lb,f.ci.sim.ub))
      } else {
        if(is.null(ylim)) ylim <-  range(c(f.fitted,f.ci.pw.lb,f.ci.pw.ub,f.ci.bf.lb,f.ci.bf.ub,f.ci.sim.lb,f.ci.sim.ub))
      }
      ## First plot, then lower and upper (lines)
      plot(y.plot.eval[order(y.plot.eval)],f.fitted[order(y.plot.eval)],
           sub=sub,
           ylim=ylim,
           ylab=ylab,
           xlab=xlab,
           type=type,
           panel.first=grid(lty=1),
           ...)
      if(ci.method == "Pointwise") {
        lines(y.plot.eval[order(y.plot.eval)],f.ci.pw.lb[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],f.ci.pw.ub[order(y.plot.eval)],lty=2)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci.method == "Bonferroni") {
        lines(y.plot.eval[order(y.plot.eval)],f.ci.bf.lb[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],f.ci.bf.ub[order(y.plot.eval)],lty=2)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci & ci.method == "Simultaneous") {
        lines(y.plot.eval[order(y.plot.eval)],f.ci.sim.lb[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],f.ci.sim.ub[order(y.plot.eval)],lty=2)
        legend("topright",legend=c("Estimated f(y|x)",paste(100*(1-alpha),"% ",ci.method, " CIs",sep="")),lty=c(1,2),bty="n")
      } else if(ci.method == "all") {
        lines(y.plot.eval[order(y.plot.eval)],f.ci.pw.lb[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],f.ci.pw.ub[order(y.plot.eval)],lty=2)
        lines(y.plot.eval[order(y.plot.eval)],f.ci.sim.lb[order(y.plot.eval)],lty=3)
        lines(y.plot.eval[order(y.plot.eval)],f.ci.sim.ub[order(y.plot.eval)],lty=3)
        lines(y.plot.eval[order(y.plot.eval)],f.ci.bf.lb[order(y.plot.eval)],lty=4)
        lines(y.plot.eval[order(y.plot.eval)],f.ci.bf.ub[order(y.plot.eval)],lty=4)
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
    ## persp=TRUE convert to matrix using x.eval & y.eval)
    return(list(bias.vec=bias.vec,
                ci.cores=ci.cores,
                f.ci.bf.lb=f.ci.bf.lb,
                f.ci.bf.ub=f.ci.bf.ub,
                f.ci.pw.lb=f.ci.pw.lb,
                f.ci.pw.ub=f.ci.pw.ub,
                f.ci.sim.lb=f.ci.sim.lb,
                f.ci.sim.ub=f.ci.sim.ub,
                F.ci.bf.lb=F.ci.bf.lb,
                F.ci.bf.ub=F.ci.bf.ub,
                F.ci.pw.lb=F.ci.pw.lb,
                F.ci.pw.ub=F.ci.pw.ub,
                F.ci.sim.lb=F.ci.sim.lb,
                F.ci.sim.ub=F.ci.sim.ub,
                g.ci.bf.lb=g.ci.bf.lb,
                g.ci.bf.ub=g.ci.bf.ub,
                g.ci.pw.lb=g.ci.pw.lb,
                g.ci.pw.ub=g.ci.pw.ub,
                g.ci.sim.lb=g.ci.sim.lb,
                g.ci.sim.ub=g.ci.sim.ub,
                f1.ci.bf.lb=f1.ci.bf.lb,
                f1.ci.bf.ub=f1.ci.bf.ub,
                f1.ci.pw.lb=f1.ci.pw.lb,
                f1.ci.pw.ub=f1.ci.pw.ub,
                f1.ci.sim.lb=f1.ci.sim.lb,
                f1.ci.sim.ub=f1.ci.sim.ub,
                F1.ci.bf.lb=F1.ci.bf.lb,
                F1.ci.bf.ub=F1.ci.bf.ub,
                F1.ci.pw.lb=F1.ci.pw.lb,
                F1.ci.pw.ub=F1.ci.pw.ub,
                F1.ci.sim.lb=F1.ci.sim.lb,
                F1.ci.sim.ub=F1.ci.sim.ub,
                g1.ci.bf.lb=g1.ci.bf.lb,
                g1.ci.bf.ub=g1.ci.bf.ub,
                g1.ci.pw.lb=g1.ci.pw.lb,
                g1.ci.pw.ub=g1.ci.pw.ub,
                g1.ci.sim.lb=g1.ci.sim.lb,
                g1.ci.sim.ub=g1.ci.sim.ub,                
                f=f.fitted,
                F=F.fitted,
                g=g.fitted,
                f1=f1.fitted,
                F1=F1.fitted,
                g1=g1.fitted,
                n.grid=n.grid,
                proper.cores=proper.cores,
                secs.elapsed=as.numeric(difftime(Sys.time(),secs.start,units="secs")),
                x.eval=x.eval,
                x=x.plot.eval,
                y.eval=y.eval,
                y=y.plot.eval,
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

## summary.bkcde() provides a summary of a bkcde object

summary.bkcde <- function(object, ...) {
  if(!inherits(object,"bkcde")) stop("object must be of class bkcde in summary.bkcde()")
  cat("Call:\n")
  cat("bkcde(h.y=",round(object$h[1],4),", h.x=",round(object$h[2],4),", x, y, x.eval, y.eval, y.lb=",object$y.lb,", y.ub=",object$y.ub,", x.lb=",object$x.lb,", x.ub=",object$x.ub,", degree=",object$degree,")\n",sep="")
  cat("\n")
  cat("Number of sample realizations: ",format(length(object$y), big.mark=",", scientific=FALSE),"\n",sep="")
  cat("Number of evaluation points: ",format(length(object$y.eval), big.mark=",", scientific=FALSE),"\n",sep="")
  cat("Bandwidths: h.y = ",object$h[1],", h.x = ",object$h[2],"\n",sep="")
  cat("Bandwidth scale factors: sf.y = ",object$h.sf[1],", sf.x = ",object$h.sf[2],"\n",sep="")
  cat("Degree of local polynomial: ",object$degree,"\n",sep="")
  cat("Warnings enabled: ",object$warnings,"\n",sep="")
  if(!is.null(object$f.yx.integral.pre.neg)) cat("Integral of estimate (pre any negativity correction): ",formatC(object$f.yx.integral.pre.neg,format="f",digits=12),"\n",sep="")
  if(!is.null(object$f.yx.integral)) cat("Integral of estimate (post negativity, prior to integration to 1 correction): ",formatC(object$f.yx.integral,format="f",digits=12),"\n",sep="")
  if(!is.null(object$f.yx.integral.post)) cat("Integral of estimate (post all corrections): ",formatC(object$f.yx.integral.post,format="f",digits=12),"\n",sep="")
  if(object$optimize) {
    cat("Bandwidth selection criterion: ",object$bwmethod,"\n",sep="")
    cat("Optimization cross-validation method: ",object$cv,"\n",sep="")
    cat("Optimization cross-validation binned: ",object$cv.binned,"\n",sep="")
    if(object$cv.binned) cat("Optimization cross-validation bin count per dimension: ",object$n.binned,"\n",sep="")
    cat("Convergence code: ",object$convergence,"\n",sep="")
    if(object$cv=="sub") {
      cat("Number of sub-cv resamples: ",object$resamples,"\n",sep="")
      cat("Sample size of sub-cv resamples: ",format(object$n.sub, big.mark=",", scientific=FALSE),"\n",sep="")
    }
    cat("Core allocation mode: ",object$optim.cores,"\n",sep="")
    cat("Number of cores used for optimization in parallel processing for degree selection: ",object$optim.degree.cores,"\n",sep="")
    cat("Number of cores used for optimization in parallel processing for multistart optimization: ",object$optim.nmulti.cores,"\n",sep="")
    cat("Total number of cores used for optimization in parallel processing: ",object$optim.ksum.cores*object$optim.degree.cores*object$optim.nmulti.cores,"\n",sep="")
  }
  if(!object$cv.only) {
    if(object$proper) cat("Number of cores used in parallel processing for ensuring proper density: ",object$proper.cores,"\n",sep="")
    cat("Number of cores used in parallel processing for fitting density: ",object$fitted.cores,"\n",sep="")
    cat("Number of cores used in parallel processing for kernel sum: ",object$optim.ksum.cores,"\n",sep="")
  }
  cat("Elapsed time (total): ",formatC(object$secs.elapsed,format="f",digits=2)," seconds\n",sep="")
  if(object$optimize & !object$cv.only & object$cv != "sub") {
    optim_time <- sum(object$secs.optim.mat)
    fit_time <- object$secs.estimate
    opt_cores <- max(1, object$optim.ksum.cores*object$optim.degree.cores*object$optim.nmulti.cores)
    fit_cores <- max(1, object$fitted.cores)
    ideal_elapsed <- optim_time/opt_cores + fit_time/fit_cores
    eff_overall <- ideal_elapsed / max(1e-9, object$secs.elapsed)
    eff_opt <- (optim_time/opt_cores) / max(1e-9, object$secs.elapsed)
    eff_fit <- ifelse(fit_time > 0, (fit_time/fit_cores) / max(1e-9, object$secs.elapsed), 0)

    cat("Optimization time: ",formatC(optim_time,format="f",digits=2)," seconds (cores = ",opt_cores,", per-core = ",formatC(optim_time/opt_cores,format="f",digits=3),")\n",sep="")
    cat("Fitting time: ",formatC(fit_time,format="f",digits=2)," seconds (cores = ",fit_cores,", per-core = ",formatC(fit_time/fit_cores,format="f",digits=3),")\n",sep="")
    cat("Stage efficiencies (opt / fit, ideal = 1): ",formatC(eff_opt,format="f",digits=2)," / ",formatC(eff_fit,format="f",digits=2),"\n",sep="")
    cat("Overall parallel efficiency (ideal = 1): ",formatC(eff_overall,format="f",digits=2),"\n",sep="")
  } else if(object$optimize & object$cv.only & object$cv != "sub") {
    optim_time <- sum(object$secs.optim.mat)
    opt_cores <- max(1, object$optim.ksum.cores*object$optim.degree.cores*object$optim.nmulti.cores)
    ideal_elapsed <- optim_time/opt_cores
    eff_overall <- ideal_elapsed / max(1e-9, object$secs.elapsed)

    cat("Optimization time: ",formatC(optim_time,format="f",digits=2)," seconds (cores = ",opt_cores,", per-core = ",formatC(ideal_elapsed,format="f",digits=3),")\n",sep="")
    cat("Parallel efficiency (optimization only, ideal = 1): ",formatC(eff_overall,format="f",digits=2),"\n",sep="")
  }
  cat("\n")
  invisible()
}

