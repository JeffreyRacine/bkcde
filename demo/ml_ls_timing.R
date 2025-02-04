## We investigate timing the speed of least-squares cross-validation
## (cv.ls) versus likelihood cross-validation (cv.ml). We set the
## number of multistarts to 1 and the polynomial degree to 1 to ensure
## we are using lm.wfit() for the comparison.

library(bkcde)
ml.time <- numeric()
ls.time <- numeric()
time.rat <- numeric()
n.vec <- seq(100,2000,50)

for(i in 1:length(n.vec)){
  set.seed(42)
  x <- rnorm(n.vec[i])
  y <- rnorm(n.vec[i],mean=x)
  ml.time[i] <- as.numeric(system.time(f.yx <- bkcde(x=x,y=y,cv="full",degree.min=1,degree.max=1,nmulti=1,bwmethod="cv.ml",proper=TRUE))["elapsed"])
  ls.time[i] <- as.numeric(system.time(f.yx <- bkcde(x=x,y=y,cv="full",degree.min=1,degree.max=1,nmulti=1,bwmethod="cv.ls",proper=TRUE))["elapsed"])
  time.rat[i] <- ml.time[i]/ls.time[i]
  if(i>2){
    par(mfrow=c(2,1))
    matplot(n.vec[1:i],cbind(ml.time,ls.time,fitted(lm(ml.time~poly(n.vec[1:i],degree=2))),fitted(lm(ls.time~poly(n.vec[1:i],degree=2))),time.rat[1:i]),type="b",lty=1,col=c(1,2,3,3,4),ylab="Time",xlab="Sample size",main="Time comparison")
    legend("topleft",legend=c("ML","LS","ML fit","LS fit","Ratio"),col=c(1,2,3,3,4),lty=1,bty="n")
    plot(n.vec[1:i],time.rat[1:i],type="b",lty=1,col=1,ylab="Ratio",xlab="Sample size",main="Ratio between ML and LS")
    abline(h=mean(ml.time)/mean(ls.time),col=2,lty=2)
    lines(n.vec[1:i],fitted(lm(time.rat[1:i]~poly(n.vec[1:i],degree=2))),col=3,lty=3)
  }
  cat("\r avg. ratio: ",mean(ml.time)/mean(ls.time))
}
