library(np)
data(Italy)
library(bkcde)

## Italian GDP growth panel for 21 regions covering the period 1951-1998
## (millions of Lire, 1990=base). There are 1008 observations in total.

x <- as.integer(as.character(Italy$year))
y <- Italy$gdp

f.yx <- bkcde(x=x,y=y,proper=TRUE)

## Create a 2D plot of the estimated conditional density for the year 1974,
## f(y|x=1974), and a 3D plot of the estimated conditional density for all
## years.

par(mfrow=c(1,2))
plot(f.yx,persp=FALSE,x.eval=1974,n.grid=500)
plot(f.yx,n.grid=50)

summary(f.yx)
