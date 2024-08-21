library(bkcde)
set.seed(42)
n <- 1e+06
ii <- 1:10000
n.grid <- 50
x <- runif(n,-1,1)
s1 <- 2
s2 <- 2
alpha <- s1+sin(pi*x)
beta <- s2-cos(pi*x)
y <- rbeta(n,alpha,beta)

f.yx <- bkcde(x=x,y=y,n.grid=n.grid,proper=TRUE,progress=TRUE)

par(mfrow=c(2,3),cex=.75)
x.seq <- sort(unique(f.yx$x.eval))
y.seq <- sort(unique(f.yx$y.eval))
alpha <- s1+sin(pi*x.seq)
beta <- s2-cos(pi*x.seq)
dgp.seq <- alpha/(alpha+beta)
alpha <- s1+sin(pi*f.yx$x.eval)
beta <- s2-cos(pi*f.yx$x.eval)
f.dgp.mat <- matrix(dbeta(f.yx$y.eval,alpha,beta),n.grid,n.grid)
F.dgp.mat <- matrix(pbeta(f.yx$y.eval,alpha,beta),n.grid,n.grid)

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
      main="Estimated Conditional Density",
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
      main="Estimated Conditional Distribution",
      ticktype="detailed",
      expand=0.75,
      shade=.25)

plot(x[ii],y[ii],cex=.1,col="lightgrey",xlab="x",ylab="y",main="Estimated Conditional Mean")
lines(x.seq,f.yx$E.yx)

summary(f.yx)
