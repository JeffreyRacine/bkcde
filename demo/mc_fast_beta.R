## This code runs a Monte Carlo experiment to compare the performance of the
## fast beta kernel conditional density estimator (fast.optim) with the full
## sample estimator (bkcde) for a beta distribution. The fast.optim function is
## used to select the bandwidth and degree of the kernel conditional density
## estimator (kde) using the covMcd and non-covMcd methods. The fast.optim
## function is called within a loop that generates a random sample from a beta
## distribution, estimates the bandwidth and degree of the kde using the
## fast.optim function, and then estimates the kde using the bkcde function. The
## RMSE of the kde is calculated for the full sample estimator and the
## fast.optim estimator using the covMcd and non-covMcd methods. The RMSE,
## bandwidth, and degree of the kde are saved to a file for each replication.
## The RMSE, bandwidth, and degree of the kde are then plotted for each
## replication.

library(bkcde)

set.seed(42)
n.sub <- c(100,200)
M <- 1000
n.grid <- 25

n <- 800 # scan("n.dat")
resamples <- 20 # scan("resamples.dat",quiet=TRUE)
nmulti <- 3 # scan("nmulti.dat",quiet=TRUE)
degree.min <- 0 # scan("degree_min.dat",quiet=TRUE)
degree.max <- 3 # scan("degree_max.dat",quiet=TRUE)
proper <- TRUE # scan("proper.dat",logical(),quiet=TRUE)

rmse.n <- numeric()
degree.n <- numeric()
h.x.n <- numeric()
h.y.n <- numeric()

rmse.covMcd.mat <- matrix(0,nrow=M,ncol=length(n.sub))
degree.covMcd.mat <- matrix(0,nrow=M,ncol=length(n.sub))
h.x.covMcd.mat <- matrix(0,nrow=M,ncol=length(n.sub))
h.y.covMcd.mat <- matrix(0,nrow=M,ncol=length(n.sub))

rmse.non.covMcd.mat <- matrix(0,nrow=M,ncol=length(n.sub))
degree.non.covMcd.mat <- matrix(0,nrow=M,ncol=length(n.sub))
h.x.non.covMcd.mat <- matrix(0,nrow=M,ncol=length(n.sub))
h.y.non.covMcd.mat <- matrix(0,nrow=M,ncol=length(n.sub))

for(m in 1:M) {
  
  cat("\rReplication ",m," of ",M,sep="")
  
  x <- runif(n,-.25,.25)
  s1 <- s2 <- 1.25
  y <- rbeta(n,s1-x,s2+x)
  
  ## Full sample estimator
  
  bkcde.n <- bkcde(x=x,
                   y=y,
                   n.grid=n.grid,
                   nmulti=nmulti,
                   proper=proper,
                   degree.min=degree.min,
                   degree.max=degree.max)
  
  ## Ensure dgp is consistent with y <- above
  
  dgp <- dbeta(bkcde.n$y.eval,s1-bkcde.n$x.eval,s2+bkcde.n$x.eval)
  
  rmse.n[m] <- sqrt(mean((bkcde.n$f - dgp)^2))
  degree.n[m] <- bkcde.n$degree
  h.x.n[m] <- bkcde.n$h[1]
  h.y.n[m] <- bkcde.n$h[2]
  
  for(j in 1:length(n.sub)) {
    
    cat("\rReplication ",m," of ",M," n.sub = ",n.sub[j],sep="")
    
    optimal.covMcd <- fast.optim(x=x,y=y,
                                 n.sub=n.sub[j],
                                 use.covMcd=TRUE,
                                 non.covMcd="mean",
                                 resamples=resamples,
                                 nmulti=nmulti,
                                 degree.min=degree.min,
                                 degree.max=degree.max)
    
    optimal.non.covMcd <- fast.optim(x=x,y=y,
                                     n.sub=n.sub[j],
                                     use.covMcd=FALSE,
                                     non.covMcd="median",
                                     resamples=resamples,
                                     nmulti=nmulti,
                                     degree.min=degree.min,
                                     degree.max=degree.max)
    
    bkcde.covMcd <- bkcde(h=optimal.covMcd$h,
                          degree=optimal.covMcd$degree,
                          x=x,
                          y=y,
                          proper=proper,
                          n.grid=n.grid)
    
    bkcde.non.covMcd <- bkcde(h=optimal.non.covMcd$h,
                              degree=optimal.non.covMcd$degree,
                              x=x,
                              y=y,
                              proper=proper,
                              n.grid=n.grid)
    
    rmse.covMcd.mat[m,j] <- sqrt(mean((bkcde.covMcd$f - dgp)^2))
    rmse.non.covMcd.mat[m,j] <- sqrt(mean((bkcde.non.covMcd$f - dgp)^2))
    h.x.covMcd.mat[m,j] <- bkcde.covMcd$h[1]
    h.y.covMcd.mat[m,j] <- bkcde.covMcd$h[2]
    
    degree.covMcd.mat[m,j] <- bkcde.covMcd$degree
    degree.non.covMcd.mat[m,j] <- bkcde.non.covMcd$degree
    h.x.non.covMcd.mat[m,j] <- bkcde.non.covMcd$h[1]
    h.y.non.covMcd.mat[m,j] <- bkcde.non.covMcd$h[2]
    
  }
  
  write(rmse.n[m],file="rmse_n.out",append=TRUE)
  write(degree.n[m],file="degree_n.out",append=TRUE)
  write(h.x.n[m],file="h_x_n.out",append=TRUE)
  write(h.y.n[m],file="h_y_n.out",append=TRUE)
  
  write(c(rmse.covMcd.mat[m,]),file="rmse_covMcd.out",ncolumns=length(n.sub),append=TRUE)
  write(c(degree.covMcd.mat[m,]),file="degree_covMcd.out",ncolumns=length(n.sub),append=TRUE)
  write(c(h.x.covMcd.mat[m,]),file="h_x_covMcd.out",ncolumns=length(n.sub),append=TRUE)
  write(c(h.y.covMcd.mat[m,]),file="h_y_covMcd.out",ncolumns=length(n.sub),append=TRUE)
  
  write(c(rmse.non.covMcd.mat[m,]),file="rmse_non_covMcd.out",ncolumns=length(n.sub),append=TRUE)
  write(c(degree.non.covMcd.mat[m,]),file="degree_non_covMcd.out",ncolumns=length(n.sub),append=TRUE)
  write(c(h.x.non.covMcd.mat[m,]),file="h_x_non_covMcd.out",ncolumns=length(n.sub),append=TRUE)
  write(c(h.y.non.covMcd.mat[m,]),file="h_y_non_covMcd.out",ncolumns=length(n.sub),append=TRUE)
  
  d.levs <- degree.min:(degree.max-degree.min+1)
  d.counts.mat <- as.matrix(table(ordered(degree.n[1:m],levels=d.levs)))
  d.covMcd.counts.mat <- matrix(0,nrow=length(d.levs),ncol=length(n.sub))
  d.non.covMcd.counts.mat <- matrix(0,nrow=length(d.levs),ncol=length(n.sub))
  for(j in 1:length(n.sub)) {
    d.covMcd.counts.mat[,j] <- as.matrix(table(ordered(degree.covMcd.mat[1:m,j],levels=d.levs)))
    d.non.covMcd.counts.mat[,j] <- as.matrix(table(ordered(degree.non.covMcd.mat[1:m,j],levels=d.levs)))
  }
  degree.mat <- cbind(d.counts.mat,d.covMcd.counts.mat,d.non.covMcd.counts.mat)
  colnames(degree.mat) <- c("n",paste("covMcd",n.sub,sep="_"),paste("non_covMcd",n.sub,sep="_"))
  rownames(degree.mat) <- d.levs
  
  write(colnames(degree.mat),file="degree_mat.out",ncolumns=ncol(degree.mat),append=TRUE)
  write(t(degree.mat),file="degree_mat.out",append=TRUE)
  
  pdf("rmse_fast_beta.pdf",pointsize=8)
  
  par(mfrow=c(2,2),cex=0.5)
  
  boxplot(data.frame(h.x.n[1:m],h.x.covMcd.mat[1:m,,drop=FALSE],h.x.non.covMcd.mat[1:m,,drop=FALSE]),
          names=c("n",paste("covMcd",n.sub,sep="_"),paste("non_covMcd",n.sub,sep="_")),
          ylab="h.x",
          main=c("Mean h.x\n",
                 paste(c(round(mean(h.x.n[1:m]),4),
                         round(colMeans(h.x.covMcd.mat[1:m,,drop=FALSE]),4),
                         round(colMeans(h.x.non.covMcd.mat[1:m,,drop=FALSE]),4)),collapse=",")),
          sub=c("Median h.x\n",
                paste(c(round(median(h.x.n[1:m]),4),
                        round(apply(h.x.covMcd.mat[1:m,,drop=FALSE],2,median),4),
                        round(apply(h.x.non.covMcd.mat[1:m,,drop=FALSE],2,median),4)),collapse=",")),
          #col=c("blue","red","green"),
          las=1,
          outline=FALSE,
          notch=TRUE)
  
  boxplot(data.frame(h.y.n[1:m],h.y.covMcd.mat[1:m,,drop=FALSE],h.y.non.covMcd.mat[1:m,,drop=FALSE]),
          names=c("n",paste("covMcd",n.sub,sep="_"),paste("non_covMcd",n.sub,sep="_")),
          ylab="h.y",
          main=c("Mean h.y\n",
                 paste(c(round(mean(h.y.n[1:m]),4),
                         round(colMeans(h.y.covMcd.mat[1:m,,drop=FALSE]),4),
                         round(colMeans(h.y.non.covMcd.mat[1:m,,drop=FALSE]),4)),collapse=",")),
          sub=c("Median h.y\n",
                paste(c(round(median(h.y.n[1:m]),4),
                        round(apply(h.y.covMcd.mat[1:m,,drop=FALSE],2,median),4),
                        round(apply(h.y.non.covMcd.mat[1:m,,drop=FALSE],2,median),4)),collapse=",")),
          #col=c("blue","red","green"),
          las=1,
          outline=FALSE,
          notch=TRUE)
  
  
  boxplot(data.frame(rmse.n[1:m],rmse.covMcd.mat[1:m,,drop=FALSE],rmse.non.covMcd.mat[1:m,,drop=FALSE]),
          names=c("n",paste("covMcd",n.sub,sep="_"),paste("non_covMcd",n.sub,sep="_")),
          ylab="RMSE",
          main=c("Mean RMSE\n",
                 paste(c(round(mean(rmse.n[1:m]),4),
                         round(colMeans(rmse.covMcd.mat[1:m,,drop=FALSE]),4),
                         round(colMeans(rmse.non.covMcd.mat[1:m,,drop=FALSE]),4)),collapse=",")),
          sub=c("Median RMSE\n",
                paste(c(round(median(rmse.n[1:m]),4),
                        round(apply(rmse.covMcd.mat[1:m,,drop=FALSE],2,median),4),
                        round(apply(rmse.non.covMcd.mat[1:m,,drop=FALSE],2,median),4)),collapse=",")),
          #col=c("blue","red","green"),
          las=1,
          outline=FALSE,
          notch=TRUE)
  
  barplot(degree.mat,
          beside=TRUE,
          #col=c("blue","red","green"),
          las=1,
          main="Degree",
          ylab="Count",
          xlab="Degree",
          legend=rownames(degree.mat),
          args.legend=list(x="topright",bty="n"))
  
  dev.off()
}
