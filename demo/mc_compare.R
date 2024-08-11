## This Monte Carlo script can be used to run simulations and also pick up where
## one left off if you wished to interrupt, or if you ran out of allocated
## runtime on a cluster. You need to modify the DGP below to suit your needs.

library(bkcde)
library(lpcde)
library(progress)
suppressPackageStartupMessages(library(parallel))

## MC parameters

n <- scan("n.dat",quiet=TRUE)
n.grid <- scan("ngrid.dat",quiet=TRUE)
M <- scan("M.dat",quiet=TRUE)
trim <- scan("trim.dat",quiet=TRUE)
nmulti <- scan("nmulti.dat",quiet=TRUE)
plot.during <- scan("plot_during.dat",logical(),quiet=TRUE)
degree.min <- scan("degree_min.dat",quiet=TRUE)
degree.max <- scan("degree_max.dat",quiet=TRUE)
plot.pdf <- scan("plot_pdf.dat",logical(),quiet=TRUE)
ksum.cores <- scan("ksum_cores.dat",quiet=TRUE)
optim.degree.cores <- scan("degree_cores.dat",quiet=TRUE)
optim.nmulti.cores <- scan("nmulti_cores.dat",quiet=TRUE)
progress <- TRUE
true.dgp <- match.arg(scan("dgp.dat",character(),quiet=TRUE),c("normal","beta"))
bounds <- match.arg(scan("bounds.dat",character(),quiet=TRUE),c("empirical","known"))

## Need these declared regardless of whether starting from scratch or picking up
## where it left off (if runtime exceeded, reboot of system, etc.)

rmse.models.mat <- matrix(NA,nrow=M,ncol=4)
cv.models.mat <- matrix(NA,nrow=M,ncol=3)
rmse.candidates.mat <- matrix(NA,nrow=M,ncol=degree.max-degree.min+1)
cv.candidates.mat <- matrix(NA,nrow=M,ncol=degree.max-degree.min+1)
degree.mat <- matrix(NA,nrow=M,ncol=3)
times.mat <- matrix(NA,nrow=M,ncol=2)

## If the output files exist and have data in them, read them in and continue

if(file.exists("rmse_models.out") && nrow(read.table("rmse_models.out",header=TRUE)) > 0) {
  m <- nrow(read.table("rmse_models.out",header=TRUE))
  rmse.models.mat[1:m,] <- data.matrix(read.table("rmse_models.out",header=TRUE))
  cv.models.mat[1:m,] <- data.matrix(read.table("cv_models.out",header=TRUE))
  rmse.candidates.mat[1:m,] <- data.matrix(read.table("rmse_candidates.out",header=TRUE))
  cv.candidates.mat[1:m,] <- data.matrix(read.table("cv_candidates.out",header=TRUE))
  degree.mat[1:m,] <- data.matrix(read.table("degree.out",header=TRUE))
  times.mat[1:m,] <- data.matrix(read.table("times.out",header=TRUE))
  m <- m+1
  aborts <- scan("aborts.out")
  if(m > M) {
    stop("All replications completed.\n")
  }
} else {
  write(c("BKPA","HRL","FYT","CCJM"),file="rmse_models.out",ncolumns=4)
  write(c("BKPA","HRL","FYT"),file="cv_models.out",ncolumns=3)
  write(paste("p=",degree.min:degree.max,sep=""),file="rmse_candidates.out",ncolumns=degree.max-degree.min+1)
  write(paste("p=",degree.min:degree.max,sep=""),file="cv_candidates.out",ncolumns=degree.max-degree.min+1)
  write(c("BKPA","HRL","FYT"),file="degree.out",ncolumns=3)
  write(c("CCJM","BKPA"),file="times.out",ncolumns=2)
  m <- 1
  aborts <- 0
}

## Main loop

set.seed(m)

if(progress) pbb <- progress_bar$new(format = "  Monte Carlo Simulation [:bar] :percent eta: :eta",
                                     clear = TRUE,
                                     force = TRUE,
                                     width = 60,
                                     total = M-m+1)

for(i in m:M) {
  
  lpcde.success <- FALSE
  
  ## Time the lpcde and data step
  
  st.lpcde <- Sys.time()
  
  while(!lpcde.success) {
    
    ## The lpcde method of Catteano et al 2024 aborts on some DGPs, so we
    ## "forgive" it and draw another "friendlier" resample for this implementation
    ## (0.1.4 as of this writing)
    
    if(true.dgp=="normal") {
      x <- runif(n,-0.5,0.5)
      y <- rnorm(n,mean=x)
      x.eval <- 0
      x.eval.grid <- rep(x.eval,n.grid)
      y.grid <- seq(quantile(y,trim),quantile(y,1-trim),length=n.grid)
      dgp <- dnorm(y.grid,mean=x.eval)
      if(bounds=="empirical") {
        y.lb <- min(y)
        y.ub <- max(y)
      } else if(bounds=="known") {
        y.lb <- -Inf
        y.ub <- Inf
      }
    } else if(true.dgp=="beta") {
      x <- runif(n,-.25,.25)
      s1 <- 1
      s2 <- 1.5
      y <- rbeta(n,s1+x,s2+x)
      x.eval <- 0
      x.eval.grid <- rep(x.eval,n.grid)
      y.grid <- seq(quantile(y,trim),quantile(y,1-trim),length=n.grid)
      dgp <- dbeta(y.grid,s1+x.eval,s2+x.eval)
      if(bounds=="empirical") {
        y.lb <- min(y)
        y.ub <- max(y)
      } else if(bounds=="known") {
        y.lb <- 0
        y.ub <- 1
      }
    }
    
    f.yx.lpcde <- tryCatch(lpcde(x_data=x, 
                                 y_data=y, 
                                 y_grid=y.grid, 
                                 x=x.eval, 
                                 bw_type="mse-rot", 
                                 nonneg=TRUE, 
                                 normalize=TRUE),
                           error = function(e){return(FALSE)})
    
    if(class(f.yx.lpcde)=="lpcde") {
      lpcde.success <- TRUE
    } else {
      aborts <- aborts+1
    }
    
  }
  
  st.lpcde <- as.numeric(difftime(Sys.time(),st.lpcde,units="secs"))
  
  ## Fit models from degree degree.min to degree.max using empirical support bounds
  
  st.bkcde <- Sys.time()
  
  output <- bkcde(h=NULL,
                  x=x,
                  y=y,
                  x.eval=x.eval.grid,
                  y.eval=y.grid,
                  degree.min=degree.min,
                  degree.max=degree.max,
                  nmulti=nmulti,
                  ksum.cores=ksum.cores,
                  optim.degree.cores=optim.degree.cores,
                  optim.nmulti.cores=optim.nmulti.cores)
  
  options(scipen=9)
  print(output$value.mat)
  foo <-data.frame(cbind(t(apply(output$value.mat,1,range)),
                         apply(output$value.mat,1,IQR),
                         apply(output$value.mat,1,max)-apply(output$value.mat,1,min)))
  colnames(foo) <- c("min","max","iqr","range")
  rownames(foo) <- paste("p=",degree.min:degree.max,sep="")
  print(foo)
  
  ## Fit the Hall et al model (p=0) and Fan et al model (p=1) with an infinite
  ## support kernel
  
  output.inf <- mclapply(0:1,function(p) {
    bkcde(h=NULL,
          x=x,
          y=y,
          x.eval=x.eval.grid,
          y.eval=y.grid,
          y.lb=-Inf,
          y.ub=Inf,
          x.lb=-Inf,
          x.ub=Inf,
          degree.min=p,
          degree.max=p,
          nmulti=nmulti,
          ksum.cores=ksum.cores,
          optim.degree.cores=1,
          optim.nmulti.cores=optim.nmulti.cores)
  },mc.cores=2)
  
  output.hrl <- output.inf[[1]]
  output.fyt <- output.inf[[2]]
  
  st.bkcde <- as.numeric(difftime(Sys.time(),st.bkcde,units="secs"))
  
  rmse.models.mat[i,] <- c(sqrt(mean((output$f-dgp)^2)),
                           sqrt(mean((output.hrl$f-dgp)^2)),
                           sqrt(mean((output.fyt$f-dgp)^2)),
                           sqrt(mean((f.yx.lpcde$Estimate[,"est_RBC"]-dgp)^2)))
  
  ## Write fitted values and rmse values from each degree degree.min,..., degree.max
  ## (these are used to construct the final bkcde() model)
  
  for(p in degree.min:degree.max) {
    f.p <- bkcde(h=output$h.mat[p-degree.min+1,],
                 x=x,
                 y=y,
                 x.eval=x.eval.grid,
                 y.eval=y.grid,
                 degree=p)$f
    rmse.candidates.mat[i,p-degree.min+1] <- sqrt(mean((f.p-dgp)^2))
    write(f.p,file=paste("f_p_",p,".out",sep=""),append=TRUE, ncolumns=n.grid)
  }
  ## Write fitted values from the bkcde and infinite support models
  write(output$f,file="f_f.yx",append=TRUE,ncolumns=n.grid)
  write(output.hrl$f,file="f_hrl.out",append=TRUE,ncolumns=n.grid)
  write(output.fyt$f,file="f_fyt.out",append=TRUE,ncolumns=n.grid)
  write(f.yx.lpcde$Estimate[,"est_RBC"],file="f_lpcde.out",append=TRUE,ncolumns=n.grid)
  ## Write data, grid, and dgp so that we can recreate the results
  write(x,file="x.out",append=TRUE,ncolumns=n)
  write(y,file="y.out",append=TRUE,ncolumns=n)
  write(dgp,file="dgp.out",append=TRUE,ncolumns=n.grid)
  write(y.grid,file="ygrid.out",append=TRUE,ncolumn=n.grid)
  ## Write the RMSE and CV values for each model and candidate, along with
  ## summary statistics
  write(rmse.models.mat[i,],file="rmse_models.out",ncolumns=4,append=TRUE)
  write(apply(rmse.models.mat[1:i,,drop=FALSE],2,mean),file="mean_rmse_models.out",ncolumns=4)
  write(apply(rmse.models.mat[1:i,,drop=FALSE],2,median),file="median_rmse_models.out",ncolumns=4)    
  degree.mat[i,] <- c(output$degree,output.hrl$degree,output.fyt$degree)
  write(degree.mat[i,],file="degree.out",ncolumns=3,append=TRUE)
  cv.models.mat[i,] <- c(output$value,output.hrl$value,output.fyt$value)
  write(cv.models.mat[i,],file="cv_models.out",ncolumns=3,append=TRUE)
  cv.candidates.mat[i,] <- output$value.vec
  write(cv.candidates.mat[i,],file="cv_candidates.out",ncolumns=degree.max-degree.min+1,append=TRUE)
  write(apply(cv.candidates.mat[1:i,,drop=FALSE],2,mean),file="mean_cv_candidates.out",ncolumns=degree.max-degree.min+1)
  write(apply(cv.candidates.mat[1:i,,drop=FALSE],2,median),file="median_cv_candidates.out",ncolumns=degree.max-degree.min+1)
  write(rmse.candidates.mat[i,],file="rmse_candidates.out",ncolumns=degree.max-degree.min+1,append=TRUE)
  write(apply(rmse.candidates.mat[1:i,,drop=FALSE],2,mean),file="mean_rmse_candidates.out",ncolumns=degree.max-degree.min+1)
  write(apply(rmse.candidates.mat[1:i,,drop=FALSE],2,median),file="median_rmse_candidates.out",ncolumns=degree.max-degree.min+1)
  write(aborts,file="aborts.out")
  write(c(st.lpcde,st.bkcde),file="times.out",ncolumns=2,append=TRUE)
  write(output$h,file="h.out",append=TRUE,ncolumns=2)
  write(output.hrl$h,file="h_hrl.out",append=TRUE,ncolumns=2)
  write(output.fyt$h,file="h_fyt.out",append=TRUE,ncolumns=2)
  write(c(output$f.yx.integral.pre.neg,output$f.yx.integral,output$f.yx.integral.post),file="integral.out",append=TRUE,ncolumns=3)
  write(c(output.fyt$f.yx.integral.pre.neg,output.fyt$f.yx.integral,output.fyt$f.yx.integral.post),file="integral_fyt.out",append=TRUE,ncolumns=3)
  
  if(plot.during) {
    if(plot.during & !plot.pdf) par(mfrow=c(2,2),cex=.5)
    if(plot.pdf) pdf(file="rmse.pdf",pointsize=7)
    ## RMSE
    df <- data.frame(rmse.candidates.mat[1:i,,drop=FALSE],rmse.models.mat[1:i,,drop=FALSE])
    names(df) <- c(paste("p=",degree.min:degree.max,sep=""),"BKPA","HRL","FYT","CCJM")
    boxplot(df,
            main=paste("RMSE Boxplots (", nrow(rmse.models.mat[1:i,,drop=FALSE])," of ",M," replications completed)\nMedian: ",paste(formatC(apply(df,2,median),format="f",digits=3),collapse=", "),sep=""),
            sub=paste("Relative Median: ", paste(formatC(apply(df,2,median)/apply(df,2,median)[7],format="f",digits=2),collapse=", "),sep=""),
            ylab="RMSE",
            outline=FALSE,
            notch=TRUE)
    abline(h=median(rmse.models.mat[1:i,1,drop=FALSE]),lty=2)
    if(plot.pdf) dev.off()
    if(plot.pdf) pdf(file="cv.pdf",pointsize=7)
    ## Cross-validation
    df <- data.frame(cv.candidates.mat[1:i,,drop=FALSE],cv.models.mat[1:i,,drop=FALSE])
    names(df) <- c(paste("p=",degree.min:degree.max,sep=""),"BKPA","HRL","FYT")
    boxplot(df,
            main=paste("CV Boxplots (", nrow(cv.models.mat[1:i,,drop=FALSE])," of ",M," replications completed)\nMedian: ",paste(formatC(apply(df,2,median),format="f",digits=3),collapse=", "),sep=""),
            sub=paste("Relative Median: ", paste(formatC(apply(df,2,median)/apply(df,2,median)[4],format="f",digits=2),collapse=", "),sep=""),
            ylab="CV",
            outline=FALSE,
            notch=TRUE)
    abline(h=median(cv.models.mat[1:i,1,drop=FALSE]),lty=2)
    if(plot.pdf) dev.off()
    if(plot.pdf) pdf(file="degree.pdf",pointsize=7)
    ## create a barchart for the counts
    barplot(prop.table(table(ordered(degree.mat[1:i,1,drop=FALSE],levels=degree.min:degree.max))), xlab="Polynomial Order", ylab="Proportion", col="lightblue", border="black")
    if(plot.pdf) dev.off()
  }
  
  if(progress) { 
    pbb$tick()
  } else {
    cat("\rm = ",i," of ",M," (",aborts," lpcde aborts)",sep="")
  }
  
}

warnings()
