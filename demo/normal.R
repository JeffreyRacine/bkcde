library(bkcde)
set.seed(42)
n <- 1e+06
ii <- 1:1e+04
n.grid <- 50
x <- runif(n,0,1)
y <- rnorm(n,mean=2*sin(4*pi*x),sd=1+abs(x))

f.yx <- bkcde(x=x,y=y,n.grid=n.grid,proper=TRUE,progress=TRUE)

f.dgp.mat <- matrix(dnorm(f.yx$y.eval,mean=2*sin(4*pi*f.yx$x.eval),sd=1+abs(f.yx$x.eval)),n.grid,n.grid)
F.dgp.mat <- matrix(pnorm(f.yx$y.eval,mean=2*sin(4*pi*f.yx$x.eval),sd=1+abs(f.yx$x.eval)),n.grid,n.grid)

x.seq <- unique(f.yx$x.eval)
y.seq <- unique(f.yx$y.eval)
dgp.seq <- 2*sin(4*pi*x.seq)

par(mfrow=c(2,3),cex=.75)

persp(y=y.seq,
      x=x.seq,
      z=f.dgp.mat,
      xlab="x",
      ylab="y",
      zlab="f(y|x)",
      zlim=range(f.dgp.mat,f.yx$f),
      theta=120,
      phi=45,
      main="True Conditional Density",
      ticktype="detailed",
      expand=0.75,
      shade=.25)

persp(y=y.seq,
      x=x.seq,
      z=F.dgp.mat,
      xlab="x",
      ylab="y",
      zlab="F(y|x)",
      zlim=range(F.dgp.mat,f.yx$F),
      theta=120,
      phi=45,
      main="True Conditional Distribution",
      ticktype="detailed",
      expand=0.75,
      shade=.25)

plot(x[ii],y[ii],cex=.1,col="lightgrey",xlab="x",ylab="y",main="True Conditional Mean")
lines(x.seq,dgp.seq)

persp(y=y.seq,
      x=x.seq,
      z=matrix(f.yx$f,n.grid,n.grid),
      xlab="x",
      ylab="y",
      zlab="f(y|x)",
      zlim=range(f.dgp.mat,f.yx$f),
      theta=120,
      phi=45,
      main="Estimaed Conditional Density",
      ticktype="detailed",
      expand=0.75,
      shade=.25)

persp(y=y.seq,
      x=x.seq,
      z=matrix(f.yx$F,n.grid,n.grid),
      xlab="x",
      ylab="y",
      zlab="F(y|x)",
      zlim=range(F.dgp.mat,f.yx$F),
      theta=120,
      phi=45,
      main="Estimaed Conditional Distribution",
      ticktype="detailed",
      expand=0.75,
      shade=.25)

plot(x[ii],y[ii],cex=.1,col="lightgrey",xlab="x",ylab="y",main="Estimated Conditional Mean")
lines(x.seq,f.yx$E.yx)

summary(f.yx)
