## This Monte Carlo script can be used to run simulations and also pick up where
## one left off if you wished to interrupt, or if you ran out of allocated
## runtime on a cluster. You need to modify the DGP below to suit your needs.

library(bkcde)
library(lpcde)
library(progress)
suppressPackageStartupMessages(library(parallel))

## MC parameters

n <- 100 # scan("n.dat",quiet=TRUE) 
n.grid <- 25 # scan("ngrid.dat",quiet=TRUE)
M <- 1000 # scan("M.dat",quiet=TRUE)
trim <- 0.00 # scan("trim.dat",quiet=TRUE)
nmulti <- 2 # scan("nmulti.dat",quiet=TRUE)
plot.during <- TRUE # scan("plot.dat",logical(),quiet=TRUE)
p.min <- 0 # scan("pmin.dat",quiet=TRUE)
p.max <- 3 # scan("pmax.dat",quiet=TRUE)
plot.pdf <- FALSE # scan("pdf.dat",logical(),quiet=TRUE)
ksum.cores <- 1 # scan("ksum_cores.dat",quiet=TRUE)
degree.cores <- 4 # scan("degree_cores.dat",quiet=TRUE)
nmulti.cores <- 2 # scan("nmulti_cores.dat",quiet=TRUE)
progress <- TRUE

## Need these declared regardless of whether starting from scratch or picking up
## where it left off (if runtime exceeded, reboot of system, etc.)

rmse.mat <- matrix(NA,nrow=M,ncol=4)
cv.mat <- matrix(NA,nrow=M,ncol=3)
degree.mat <- matrix(NA,nrow=M,ncol=3)
times.mat <- matrix(NA,nrow=M,ncol=2)

## If the output files exist and have data in them, read them in and continue

if(file.exists("rmse.out") && nrow(read.table("rmse.out",header=TRUE)) > 0) {
  m <- nrow(read.table("rmse.out",header=TRUE))
  rmse.mat[1:m,] <- data.matrix(read.table("rmse.out",header=TRUE))
  cv.mat[1:m,] <- data.matrix(read.table("cv.out",header=TRUE))
  degree.mat[1:m,] <- data.matrix(read.table("degree.out",header=TRUE))
  times.mat[1:m,] <- data.matrix(read.table("times.out",header=TRUE))
  m <- m+1
  aborts <- scan("aborts.out")
} else {
  write(c("BKPA","HRL","FYT","CCJM"),file="rmse.out",ncolumns=4)
  write(c("BKPA","HRL","FYT"),file="cv.out",ncolumns=3)
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

    x.min <- -.25
    x.max <- .25
    x <- runif(n,x.min,x.max)
    x.eval <- 0
    x.eval.grid <- rep(x.eval,n.grid)
    s1 <- 1
    s2 <- 1.5
    y <- rbeta(n,s1+x,s2+x)
    y.grid <- seq(quantile(y,trim),quantile(y,1-trim),length=n.grid)
    dgp <- dbeta(y.grid,s1+x.eval,s2+x.eval)      
    
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
  
  ## Fit models from degree p.min to p.max using empirical support bounds
  
  st.bkcde <- Sys.time()
  
  output <- bkcde(h=NULL,
                  x=x,
                  y=y,
                  x.eval=x.eval.grid,
                  y.eval=y.grid,
                  degree.min=p.min,
                  degree.max=p.max,
                  nmulti=nmulti,
                  ksum.cores=ksum.cores,
                  degree.cores=degree.cores,
                  nmulti.cores=nmulti.cores)
  
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
          degree.cores=1,
          nmulti.cores=nmulti.cores)
  },mc.cores=2)
  
  output.hrl <- output.inf[[1]]
  output.fyt <- output.inf[[2]]
  
  st.bkcde <- as.numeric(difftime(Sys.time(),st.bkcde,units="secs"))
  
  rmse.mat[i,] <- c(sqrt(mean((output$f-dgp)^2)),
                    sqrt(mean((output.hrl$f-dgp)^2)),
                    sqrt(mean((output.fyt$f-dgp)^2)),
                    sqrt(mean((f.yx.lpcde$Estimate[,"est_RBC"]-dgp)^2)))
  write(rmse.mat[i,],file="rmse.out",ncolumns=4,append=TRUE)
  write(apply(rmse.mat[1:i,,drop=FALSE],2,mean),file="mean_rmse.out",ncolumns=4)
  write(apply(rmse.mat[1:i,,drop=FALSE],2,median),file="median_rmse.out",ncolumns=4)    
  degree.mat[i,] <- c(output$degree,output.hrl$degree,output.fyt$degree)
  write(degree.mat[i,],file="degree.out",ncolumns=3,append=TRUE)
  cv.mat[i,] <- c(output$value,output.hrl$value,output.fyt$value)
  write(cv.mat[i,],file="cv.out",ncolumns=3,append=TRUE)
  write(aborts,file="aborts.out")
  write(c(st.lpcde,st.bkcde),file="times.out",ncolumns=2,append=TRUE)
  
  ## Write fitted values and rmse values from each degree p.min,..., p.max
  ## (these are used to construct the final bkcde() model)
  
  for(p in p.min:p.max) {
    f.p <- bkcde(h=output$h.mat[p-p.min+1,],
                 x=x,
                 y=y,
                 x.eval=x.eval.grid,
                 y.eval=y.grid,
                 degree.min=p,
                 degree.max=p)$f
    rmse.p <- sqrt(mean((f.p-dgp)^2))
    write(f.p,file=paste("f_p_",p,".out",sep=""),append=TRUE, ncolumns=n.grid)
    write(rmse.p,file=paste("rmse_p_",p,".out",sep=""),append=TRUE)
  }
  ## Write fitted values from the bkcde and infinite support models
  write(output$f,file="f_bkcde.out",append=TRUE,ncolumns=n.grid)
  write(output.hrl$f,file="f_hrl.out",append=TRUE,ncolumns=n.grid)
  write(output.fyt$f,file="f_fyt.out",append=TRUE,ncolumns=n.grid)
  write(f.yx.lpcde$Estimate[,"est_RBC"],file="f_lpcde.out",append=TRUE,ncolumns=n.grid)
  ## Write data, grid, and dgp so that we can recreate the results
  write(x,file="x.out",append=TRUE,ncolumns=n)
  write(y,file="y.out",append=TRUE,ncolumns=n)
  write(dgp,file="dgp.out",append=TRUE,ncolumns=n.grid)
  write(y.grid,file="ygrid.out",append=TRUE,ncolumn=n.grid)
  
  if(plot.during) {
    if(plot.pdf) pdf()
    boxplot(rmse.mat[1:i,,drop=FALSE],
            names=c("BKPA","HRL","FYT","CCJM"),
            col=c("red","green","blue","purple"),
            main=c("Mean RMSE",paste(formatC(apply(rmse.mat[1:i,,drop=FALSE],2,mean),format="f",digits=3),collapse=", ")),
            sub=paste(formatC(apply(rmse.mat[1:i,,drop=FALSE],2,median),format="f",digits=3),collapse=", "),
            ylab="RMSE",
            outline=FALSE,
            notch=TRUE)
    abline(h=median(rmse.mat[1:i,1,drop=FALSE]),lty=2)
    if(plot.pdf) dev.off()
  }
  
  if(progress) { 
    pbb$tick()
  } else {
    cat("\rm = ",i," of ",M," (",aborts," lpcde aborts)",sep="")
  }
  
}
