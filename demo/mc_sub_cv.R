## Monte Carlo R code for sub-sampling cross-validation to compare the full and
## sub-sampling cross-validation methods for bandwidth selection in bivariate
## kernel conditional density estimation.

library(bkcde)
library(progress)

n <- scan("n.dat",quiet=TRUE)
n.sub <- seq(100,400,by=100)
if(n <= n.sub[1]) stop("n must be greater than ",n.sub[1])
n.sub <- n.sub[n.sub<n]
M <- scan("M.dat",quiet=TRUE)
resamples <- scan("resamples.dat",quiet=TRUE)
nmulti <- scan("nmulti.dat",quiet=TRUE)
degree.min <- scan("degree_min.dat",quiet=TRUE)
degree.max <- scan("degree_max.dat",quiet=TRUE)
n.grid <- scan("n_grid.dat",quiet=TRUE)
true.dgp <- match.arg(scan("dgp.dat",character(),quiet=TRUE),c("normal","beta"))
bounds <- match.arg(scan("bounds.dat",character(),quiet=TRUE),c("empirical","known"))
plot.pdf <- scan("plot_pdf.dat",logical(),quiet=TRUE)
progress.bar <- scan("progress_bar.dat",logical(),quiet=TRUE)

rmse.mat <- matrix(numeric(),nrow=M,ncol=length(n.sub)+1)
degree.mat <- matrix(numeric(),nrow=M,ncol=length(n.sub)+1)
h.x.mat <- matrix(numeric(),nrow=M,ncol=length(n.sub)+1)
h.y.mat <- matrix(numeric(),nrow=M,ncol=length(n.sub)+1)
sf.x.mat <- matrix(numeric(),nrow=M,ncol=length(n.sub)+1)
sf.y.mat <- matrix(numeric(),nrow=M,ncol=length(n.sub)+1)
secs.mat <- matrix(numeric(),nrow=M,ncol=length(n.sub)+1)

if(file.exists("rmse.out") && nrow(read.table("rmse.out",header=TRUE)) > 0) {
  m <- nrow(read.table("rmse.out",header=FALSE))
  rmse.mat[1:m,] <- data.matrix(read.table("rmse.out",header=FALSE))
  degree.mat[1:m,] <- data.matrix(read.table("degree.out",header=FALSE))
  h.x.mat[1:m,] <- data.matrix(read.table("h.x.out",header=FALSE))
  h.y.mat[1:m,] <- data.matrix(read.table("h.y.out",header=FALSE))
  sf.x.mat[1:m,] <- data.matrix(read.table("sf.x.out",header=FALSE))
  sf.y.mat[1:m,] <- data.matrix(read.table("sf.y.out",header=FALSE))
  secs.mat[1:m,] <- data.matrix(read.table("secs.out",header=FALSE))
  m.start <- m+1
} else {
  m.start <- 1
}
if(m.start > M) {
  stop("Monte Carlo completed.\n")
}

if(progress.bar) pbb <- progress_bar$new(format = "  Monte Carlo Simulation [:bar] :percent eta: :eta",
                                         clear = TRUE,
                                         force = TRUE,
                                         width = 60,
                                         total = M)

for(m in m.start:M) {
  set.seed(m)
  if(!progress.bar) cat("\r",m," of ",M,"     ")
  
  if(true.dgp=="normal") {
    x <- runif(n,-0.5,0.5)
    y <- rnorm(n,mean=x)
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
    if(bounds=="empirical") {
      y.lb <- min(y)
      y.ub <- max(y)
    } else if(bounds=="known") {
      y.lb <- 0
      y.ub <- 1
    }
  }
  
  f.yx.full <- bkcde(x=x, 
                     y=y, 
                     y.lb=y.lb,
                     y.ub=y.ub,
                     degree.min=degree.min,
                     degree.max=degree.max,
                     nmulti=nmulti,
                     n.grid=n.grid,
                     cv="full")
  
  if(true.dgp=="normal") {
    dgp <- dnorm(f.yx.full$y.eval,mean=f.yx.full$x.eval)
  } else if(true.dgp=="beta") {
    dgp <- dbeta(f.yx.full$y.eval,s1+f.yx.full$x.eval,s2+f.yx.full$x.eval)
  }
  
  rmse.mat[m,1] <- sqrt(mean((f.yx.full$f-dgp)^2))
  degree.mat[m,1] <- f.yx.full$degree
  h.x.mat[m,1] <- f.yx.full$h[2]
  h.y.mat[m,1] <- f.yx.full$h[1]
  sf.x.mat[m,1] <- f.yx.full$h.sf[2]
  sf.y.mat[m,1] <- f.yx.full$h.sf[1]
  secs.mat[m,1] <- f.yx.full$secs.elapsed
  
  for(j in 1:length(n.sub)) {
    if(!progress.bar) cat("\r",m," of ",M," (",j," of ",length(n.sub),")",sep="")
    f.yx.sub <- bkcde(x=x, 
                      y=y, 
                      y.lb=y.lb,
                      y.ub=y.ub,
                      cv="sub", 
                      degree.min=degree.min,
                      degree.max=degree.max,
                      nmulti=nmulti,
                      n.grid=n.grid,
                      n.sub=n.sub[j], 
                      resamples=resamples)
    rmse.mat[m,j+1] <- sqrt(mean((f.yx.sub$f-dgp)^2))
    degree.mat[m,j+1] <- f.yx.sub$degree
    h.x.mat[m,j+1] <- f.yx.sub$h[2]
    h.y.mat[m,j+1] <- f.yx.sub$h[1]
    sf.x.mat[m,j+1] <- f.yx.sub$h.sf[2]
    sf.y.mat[m,j+1] <- f.yx.sub$h.sf[1]
    secs.mat[m,j+1] <- f.yx.sub$secs.elapsed
  }
  
  write(rmse.mat[m,],file="rmse.out",ncolumns=length(n.sub)+1,append=TRUE)
  write(degree.mat[m,],file="degree.out",ncolumns=length(n.sub)+1,append=TRUE)
  write(h.x.mat[m,],file="h.x.out",ncolumns=length(n.sub)+1,append=TRUE)
  write(h.y.mat[m,],file="h.y.out",ncolumns=length(n.sub)+1,append=TRUE)
  write(sf.x.mat[m,],file="sf.x.out",ncolumns=length(n.sub)+1,append=TRUE)
  write(sf.y.mat[m,],file="sf.y.out",ncolumns=length(n.sub)+1,append=TRUE)
  write(secs.mat[m,],file="secs.out",ncolumns=length(n.sub)+1,append=TRUE)
  
  write(apply(rmse.mat[1:m,,drop=FALSE],2,mean),ncolumns=length(n.sub)+1,file="mean_rmse.out")
  write(apply(rmse.mat[1:m,,drop=FALSE],2,median),ncolumns=length(n.sub)+1,file="median_rmse.out")
  
  write(apply(degree.mat[1:m,,drop=FALSE],2,mean),ncolumns=length(n.sub)+1,file="mean_degree.out")
  write(apply(degree.mat[1:m,,drop=FALSE],2,median),ncolumns=length(n.sub)+1,file="median_degree.out")
  
  write(apply(h.x.mat[1:m,,drop=FALSE],2,mean),ncolumns=length(n.sub)+1,file="mean_h.x.out")
  write(apply(h.x.mat[1:m,,drop=FALSE],2,median),ncolumns=length(n.sub)+1,file="median_h.x.out")
  
  write(apply(sf.x.mat[1:m,,drop=FALSE],2,mean),ncolumns=length(n.sub)+1,file="mean_sf.x.out")
  write(apply(sf.x.mat[1:m,,drop=FALSE],2,median),ncolumns=length(n.sub)+1,file="median_sf.x.out")
  
  write(apply(h.y.mat[1:m,,drop=FALSE],2,mean),ncolumns=length(n.sub)+1,file="mean_h.y.out")
  write(apply(h.y.mat[1:m,,drop=FALSE],2,median),ncolumns=length(n.sub)+1,file="median_h.y.out")
  
  write(apply(sf.y.mat[1:m,,drop=FALSE],2,mean),ncolumns=length(n.sub)+1,file="mean_sf.y.out")
  write(apply(sf.y.mat[1:m,,drop=FALSE],2,median),ncolumns=length(n.sub)+1,file="median_sf.y.out")
  
  write(apply(secs.mat[1:m,,drop=FALSE],2,mean),ncolumns=length(n.sub)+1,file="mean_secs.out")
  write(apply(secs.mat[1:m,,drop=FALSE],2,median),ncolumns=length(n.sub)+1,file="median_secs.out")
  
  colnames(rmse.mat) <- c(n,n.sub)
  if(plot.pdf) pdf(file="rmse.pdf",pointsize = 12)
  par(mfrow=c(1,1),cex=.8)
  boxplot(rmse.mat[1:m,,drop=FALSE],
          outline=FALSE,
          notch=TRUE,
          ylab="RMSE",
          xlab="Sample Size (n on left then n.sub)",
          main=c(paste(true.dgp," (",bounds,", m = ",m," of ",M,", res = ",resamples,", nmul = ",nmulti,", d.max = ",degree.max,")",sep=""),
                 paste(round(apply(rmse.mat[1:m,,drop=FALSE],2,mean),4),collapse=", ")),
          sub=paste(round(apply(rmse.mat[1:m,,drop=FALSE],2,median),4),collapse=", "),
          names=c(n,n.sub))  
  if(plot.pdf) dev.off()
  
  if(progress.bar) pbb$tick()
  
}
