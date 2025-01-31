library(bkcde)
library(np)

## Canadian cross-section wage data consisting of a random sample taken from the
## 1971 Canadian Census Public Use Tapes for male individuals having common
## education (grade 13). There are 205 observations in total.

data(cps71)
x <- cps71$age
y <- cps71$logwage

## npcdens() in the np package is a local polynomial of degree 0 estimator of
## f(y|x) that selects the bandwidth via likelihood cross-validation

np.out <- npcdens(y~x)
summary(np.out)

## bkcde() in the bkcde package selects the degree of the polynomial and
## bandwidths via likelihood cross-validation

f.yx <- bkcde(x=x,y=y,degree.max=4,nmulti=5,proper=TRUE)

## Create a 2x2 plot of the npcdens() and bkcde() figures

par(mfrow=c(2,2))

## Manually select theta and phi, no confidence intervals

plot(np.out,view="fixed",theta=-50,phi=30,main="")
plot(f.yx,proper=TRUE,n.grid=50,theta=-50,phi=30,col="cyan",main="")

## Now just look at the bkcde() figure, produce bootstrapped confidence
## intervals and report progress (default theta and phi, note two plots created
## since ci.preplot = TRUE is the default)

plot(f.yx,ci=TRUE,theta=-50,phi=30,B=100,n.grid=50,progress=TRUE)

summary(f.yx)
