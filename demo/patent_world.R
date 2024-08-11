library(bkcde)
suppressPackageStartupMessages(library(np))
options(np.messages=FALSE)
suppressPackageStartupMessages(library(pglm))
## Patent data often has a preponderance of zeros and "overdispersion". Since
## the data has a boundary at zero, a zero-inflated negative binomial model is
## often used. Instead we consider a kernel density estimation approach that
## adapts to the boundary automatically and also chooses the polynomial order in
## a data-driven fashion.
data(PatentsRD)
attach(PatentsRD)

x <- as.integer(as.character(year))
y <- patent

f.yx <- bkcde(x=x,y=y)

## Estimator of Hall, Racine & Li (polynomial order 0 conditional density estimate)
f.np <- npcdens(txdat=x,exdat=rep(1987,100),tydat=y,eydat=seq(min(y),100,length=100))

## In the following plots note that the data has a very long tail (right skewed)
## so we hone in on firms having 100 or fewer patents/year. We first take
## univariate histogram and density() estimators for the year 1987, then overlay
## the conditional kernel density estimates from npcdens() and bkcde() for the
## same year (these use data for all years). We also plot the bkcde() 2D
## estimates for the year 1987 without, then with simultaneous confidence
## intervals. Finally, we plot the bkcde() estimate in 3D with and without
## confidence intervals (all years

par(mfrow=c(2,2),cex=.7)

hist(y[x==1987],xlab="patent",main="",ylim=c(0,0.06),xlim=c(0,100),prob=TRUE,breaks=500)
lines(density(y[x==1987],from=0,to=100),lty=1,col=1)
lines(seq(min(y),100,length=100),fitted(f.np),col=2,lty=2,lwd=2)
lines(seq(min(y),100,length=100),predict(f.yx,newdata=data.frame(x=rep(1987,100),y=seq(min(y),100,length=100))),col=3,lty=3,lwd=3)
legend("topright",legend=c("density() kernel density","npcdens() conditional kernel","bckde() conditional kernel"),col=1:3,lty=1:3,lwd=c(1,2,3),bty="n")

plot(f.yx,plot.3D=FALSE,x.eval=1987,plot.2D.y.grid=seq(min(y),100,length=100),xlim=c(0,100))

plot(f.yx,plot.3D=FALSE,x.eval=1987,plot.2D.y.grid=seq(min(y),100,length=100),xlim=c(0,100),ci=TRUE,ci.method="Simultaneous",ci.preplot=FALSE,B=1000)

plot(f.yx,plot.3D=TRUE,plot.3D.x.grid=seq(min(x),max(x),length=25),plot.3D.y.grid=seq(min(y),100,length=25),proper=FALSE,ci=TRUE,ci.method="Simultaneous",ci.preplot=FALSE,B=1000)
