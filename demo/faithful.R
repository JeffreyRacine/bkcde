data(faithful)
library(bkcde)

## A data frame with 272 observations on 2 variables, waiting and eruptions (in
## minutes), the waiting time between eruptions and the duration of the eruption
## for the Old Faithful geyser in Yellowstone National Park, Wyoming, USA. This
## is a starkly bimodal bivariate distribution.

x <- faithful$waiting
y <- faithful$eruptions
f.yx <- bkcde(x=x,y=y,proper=TRUE)

## Create a 2x2 plot of the npcdens() and bkcde() figures

par(mfrow=c(2,2),cex=.6)

plot(y,x,ylab="waiting",xlab="eruptions",col = ifelse(y < 3,'black','red'))
plot(f.yx,n.grid=50,theta=30,phi=55,xlab="waiting",ylab="eruptions",expand=.75)

## Plot two 2D slices for different values of x (waiting = 50 and 85, the
## approximate centers of each mode)

plot(f.yx,persp=FALSE,x.eval=50,xlab="eruptions")
par(new=TRUE)
plot(f.yx,persp=FALSE,x.eval=85,lty=2,col=2,sub="",xlab="",axes=FALSE)
legend("topright",legend=c("f(y|x=50)","f(y|x=85)"),lty=c(1,2),col=c(1,2),bty="n")

summary(f.yx)
