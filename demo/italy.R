library(np)
data(Italy)
library(bkcde)

set.seed(42)
n.grid <- 100

x <- as.integer(as.character(Italy$year))
x.eval <- uocquantile(x,0.5)
x.eval.grid <- rep(x.eval,n.grid)
y <- Italy$gdp
y.grid <- seq(min(y),max(y),length=n.grid)
dgp <- dnorm(y.grid,mean=x.eval)

bkcde.out <- bkcde(x=x,
                   y=y,
                   x.eval=x.eval.grid,
                   y.eval=y.grid)

write(bkcde.out$h, file="bkcde_Italy_h.out",ncol=2)
write(bkcde.out$degree, file="bkcde_Italy_degree.out")
save(bkcde.out, file="bkcde_Italy.RData")

bkcde.out$value.mat

summary(bkcde.out)

plot(bkcde.out)

plot(bkcde.out, ci = TRUE, alpha = 0.01)

x.seq <- seq(min(x),max(x),length=200)
y.seq <- seq(min(y),max(y),length=length(x.seq))
gdp_x.mat <- matrix(NA,nrow=length(x.seq),ncol=length(x.seq))
for(i in 1:length(x.seq)){
  gdp_x.mat[,i] <- predict(bkcde.out, newdata=data.frame(y=y.seq,x=rep(x.seq[i],length(x.seq))))
}

plot(y.seq,gdp_x.mat[,1],ylim=range(gdp_x.mat),type="l",col=1)
legend("topright",legend=x.seq,col=1:length(x.seq),lty=1:length(x.seq),title="Year",cex=.5,bty="n")
for(i in 2:length(x.seq)){
  Sys.sleep(1)
  lines(y.seq,gdp_x.mat[,i],col=i,lty=i)
}
par(mfrow=c(2,1))
for(i in 1:length(x.seq)){
  plot(y.seq,gdp_x.mat[,i],ylim=range(gdp_x.mat),ylab="f(gdp|year)",xlab="GDP",type="l",col=1)
  legend("topright",legend=round(x.seq[i],0),col=1,lty=1,title="Year",bty="n")
  Sys.sleep(.1)
}

x.seq <- seq(min(x),max(x),length=100)
y.seq <- seq(min(y),max(y),length=length(x.seq))
gdp_x.mat <- matrix(NA,nrow=length(x.seq),ncol=length(x.seq))
for(i in 1:length(x.seq)){
  gdp_x.mat[,i] <- predict(bkcde.out, newdata=data.frame(y=y.seq,x=rep(x.seq[i],length(x.seq))))
}
for(i in 1:length(x.seq)){
  plot(y.seq,gdp_x.mat[,i],ylim=range(gdp_x.mat),ylab="f(gdp|year)",xlab="GDP",type="l",col=1)
  legend("topright",legend=round(x.seq[i],0),col=1,lty=1,title="Year",bty="n")
  Sys.sleep(.1)
}
persp(x=x.seq,y=y.seq,z=gdp_x.mat,xlab="Year",ylab="GDP",zlab="f(gdp|year)",theta=-35,phi=45,ticktype="detailed")
#persp(x=x.seq,y=y.seq,z=gdp_x.mat,xlab="Year",ylab="GDP",zlab="f(gdp|year)",theta=150,phi=45,ticktype="detailed")

par(mfrow=c(2,2),cex=.5)
x.plot <- c(5,15,35,45)
for(i in x.plot){
  plot(y.seq,gdp_x.mat[,i],type="l",ylim=range(gdp_x.mat),sub=paste("Year = ",x.seq[i],sep=""),col=i,lty=i)
}


