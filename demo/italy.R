library(np)
data(Italy)
library(bkcde)

x <- as.integer(as.character(Italy$year))
y <- Italy$gdp

f.yx <- bkcde(x=x,y=y)

summary(f.yx)

plot(f.yx,plot.3D=FALSE,x.eval=np::uocquantile(x,prob=0.5))
plot(f.yx,plot.3D=TRUE,plot.3D.n.grid=50)

