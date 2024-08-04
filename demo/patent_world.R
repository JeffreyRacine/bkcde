library(pglm)
## Patent data often has a preponderance of zeros and "overdispersion". Since
## the data has a boundary at zero, a zero-inflated negative binomial model is
## often used. Instead we consider a kernel density estimation approach that
## adapts to the boundary automatically and also chooses the polynomial order in
## a data-driven fashion.
data(PatentsRD)
attach(PatentsRD)

library(bkcde)
x <- as.integer(as.character(year))
y <- patent

bkcde.out <- bkcde(x=x,y=y)

summary(bkcde.out)

## The data has a very long tail so we hone in on firms having 100 or fewer
## patents/year

plot(bkcde.out,plot.3D=FALSE,x.eval=1984,plot.2D.y.grid=seq(min(y),100,length=100),xlim=c(0,100),ci=TRUE,ci.method="all",B=100)
plot(bkcde.out,plot.3D=TRUE,plot.3D.x.grid=seq(min(x),max(x),length=100),plot.3D.y.grid=seq(min(y),100,length=100),ci=TRUE,ci.method="Simultaneous",B=100)
