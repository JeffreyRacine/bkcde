library(pglm)
## Patent data often has a preponderance of zeros and "overdispersion". Since
## the data has a boundary at zero, a zero-inflated negative binomial model is
## often used. Instead we consider a kernel density estimation approach that
## adapts to the boundary automatically and also chooses the polynomial order in
## a data-driven fashion.
data(PatentsRD)
attach(PatentsRD)
# plot(year, patent, type="l", col="blue", lwd=2, xlab="Year", ylab="Number of patent")
barplot(table(patent), col="blue", xlab="Number of patent", ylab="Frequency")

library(bkcde)
n.grid <- 2500
x <- as.integer(as.character(year))
x.eval <- np::uocquantile(x,0.5)
x.eval.grid <- rep(x.eval,n.grid)
y <- patent
y.grid <- seq(min(y),max(y),length=n.grid)

## Boundary kernel with empirical support (default)

bkcde.out <- bkcde(x=x,
                   y=y,
                   x.eval=x.eval.grid,
                   y.eval=y.grid)

summary(bkcde.out)

par(mfrow=c(1,2))

hist(patent[year==x.eval],xlab="patent",main="",ylim=c(0,0.06),xlim=c(0,100),prob=TRUE,breaks=100)
lines(y.grid,fitted(bkcde.out),lwd=2)
lines(density(patent[year==1987]),col=2)
legend("topright",legend=c("bkcde","density"),col=c(1,2),lwd=2,bty="n")

x.seq <- seq(min(x),max(x),length=100)
y.seq <- seq(min(y),100,length=length(x.seq))
patents_x.mat <- matrix(NA,nrow=length(x.seq),ncol=length(x.seq))
for(i in 1:length(x.seq)){
  patents_x.mat[i,] <- predict(bkcde.out, newdata=data.frame(y=y.seq,x=rep(x.seq[i],length(x.seq))))
}
persp(x=x.seq,y=y.seq,z=patents_x.mat,xlab="Year",ylab="patents",zlab="f(patents|year)",theta=120,phi=25,ticktype="detailed")

## HRL

bkcde_hrl.out <- bkcde(x=x,
                       y=y,
                       x.eval=x.eval.grid,
                       y.eval=y.grid,
                       degree.min=0,
                       degree.max=0,
                       x.lb=-Inf,
                       x.ub=Inf,
                       y.lb=-Inf,
                       y.ub=Inf)

summary(bkcde_hrl.out)

par(mfrow=c(1,2))

hist(patent[year==x.eval],xlab="patent",main="",ylim=c(0,0.06),xlim=c(0,100),prob=TRUE,breaks=100)
lines(y.grid,fitted(bkcde_hrl.out),lwd=2)
lines(density(patent[year==1987]),col=2)
legend("topright",legend=c("hrl","density"),col=c(1,2),lwd=2,bty="n")

x.seq <- seq(min(x),max(x),length=100)
y.seq <- seq(min(y),100,length=length(x.seq))
patents_x.mat <- matrix(NA,nrow=length(x.seq),ncol=length(x.seq))
for(i in 1:length(x.seq)){
  patents_x.mat[i,] <- predict(bkcde_hrl.out, newdata=data.frame(y=y.seq,x=rep(x.seq[i],length(x.seq))))
}
persp(x=x.seq,y=y.seq,z=patents_x.mat,xlab="Year",ylab="patents",zlab="f(patents|year)",theta=120,phi=25,ticktype="detailed")

## FYT

bkcde_fyt.out <- bkcde(x=x,
                       y=y,
                       x.eval=x.eval.grid,
                       y.eval=y.grid,
                       degree.min=1,
                       degree.max=1,
                       x.lb=-Inf,
                       x.ub=Inf,
                       y.lb=-Inf,
                       y.ub=Inf)

summary(bkcde_fyt.out)

par(mfrow=c(1,2))

hist(patent[year==x.eval],xlab="patent",main="",ylim=c(0,0.06),xlim=c(0,100),prob=TRUE,breaks=100)
lines(y.grid,fitted(bkcde_fyt.out),lwd=2)
lines(density(patent[year==1987]),col=2)
legend("topright",legend=c("fyt","density"),col=c(1,2),lwd=2,bty="n")

x.seq <- seq(min(x),max(x),length=100)
y.seq <- seq(min(y),100,length=length(x.seq))
patents_x.mat <- matrix(NA,nrow=length(x.seq),ncol=length(x.seq))
for(i in 1:length(x.seq)){
  patents_x.mat[i,] <- predict(bkcde_fyt.out, newdata=data.frame(y=y.seq,x=rep(x.seq[i],length(x.seq))))
}
persp(x=x.seq,y=y.seq,z=patents_x.mat,xlab="Year",ylab="patents",zlab="f(patents|year)",theta=120,phi=25,ticktype="detailed")

save.image("patent_world.RData")
