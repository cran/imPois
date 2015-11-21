# This demo describes Bayesian updating for a limiting case.
# 

library(imPois)
# source("./../R/imPois.R")
# source("./../R/visualize.R")

e1 <- seq(from=0.1, to=10, by=0.2) # shape 
e0 <- seq(from=0.1, to=10, by=0.2) # rate 
nps <- expand.grid(e1=e1,e0=e0)

# 2d: log-gamma
nps$et <- apply(nps, 1, function(x) evfn(y=numeric(0), pars=c(0,x[1],x[2]), ztrunc=FALSE)$value)

lambda <- 1
n <- 1e1
set.seed(16979238) # 1,1,1,0,0,2,1,0,2,
y <- rpois(n=n, lambda=lambda)

contour(x=e1, y=e0, z=matrix(nps$et, nrow=length(e1), ncol=length(e0)), levels=seq(from=-3, to=3, by=0.5), method="edge", labcex=1.5, lty="dashed", xlab=expression(xi[1]), ylab=expression(xi[0]), cex.lab=1.5, xlim=c(0,10), ylim=c(0,10))
polygon(x=c(0,0,1,1),y=c(0,1,1,0), border="darkblue")
polygon(x=c(0,1),y=c(1,0), lwd=2)

convh <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 

for(i in 1:n){
	op <- update(convh, y=y[1:i])
	polygon(op$vtx1)
}
