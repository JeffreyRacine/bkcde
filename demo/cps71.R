library(bkcde)
library(np)

## Take the cps71 data from the R package np (see ?cps71 for details)

data(cps71)
x <- cps71$age
y <- cps71$logwage

## npcdens() in the np package is a local polynomial of degree 0 estimator of
## f(y|x) that selects the bandwidth via likelihood cross-validation

np.out <- npcdens(y~x)
summary(np.out)

## bkcde() in the bkcde package selects the degree of the polynomial and
## bandwidths via likelihood cross-validation

f.yx <- bkcde(x=x,y=y,degree.max=4,nmulti=5)

## Create a 2x2 plot of the npcdens() and bkcde() figures

par(mfrow=c(2,2))

## Manually select theta and phi, no confidence intervals

plot(np.out,view="fixed",theta=-50,phi=30,main="")
plot(f.yx,proper=TRUE,plot.3D.n.grid=50,theta=-50,phi=30,col="cyan",main="")

## Now just look at the bkcde() figure, produce bootstrapped confidence
## intervals and report progress (default theta and phi, note twp plots created
## since ci.preplot = TRUE is the default)

plot(f.yx,ci=TRUE,B=100,plot.3D.n.grid=50,progress=TRUE)

summary(f.yx)