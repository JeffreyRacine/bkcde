library(np)
data(Italy)
library(bkcde)

x <- as.integer(as.character(Italy$year))
y <- Italy$gdp

bkcde.out <- bkcde(x=x,y=y)

summary(bkcde.out)

plot(bkcde.out,plot.3D=FALSE,x.eval=np::uocquantile(x,prob=0.5))
plot(bkcde.out,plot.3D=TRUE,plot.3D.n.grid=50)

