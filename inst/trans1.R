# This demo describes Bayesian updating 
#
# Notes on 2015.11.21
# devel is required; then, it will be moved from ./inst to ./demo

rm(list=ls())

library(imPois)
library(lattice)

# source("./../R/imPois.R")
# source("./../R/visualize.R")

e2 <- seq(from=0.1, to=3, length=30)  # 0.01 is possible, 0.05 is good as well 
e1 <- seq(from=-5, to=10, length=31)
e0 <- seq(from=0, to=10, length=11)
nps <- expand.grid(e2=e2,e1=e1,e0=e0)
# nrow(nps)

z <- apply(nps,1,function(x){ 
	rv <- evfn(y=numeric(0), pars=c(x[1],x[2],x[3]), ztrunc=FALSE)$value
	# cat("x=[",x,"] val = ", rv, "\n")
	return(rv)
})

nps$et <- z
lambda <- 1
n <- 1e2

set.seed(16979238) # 1,1,1,0,0,2,1,0,2,
y <- rpois(n=n, lambda=lambda)[1:10]

nps2 <- nps[which(as.numeric(as.character(nps$e2))==0.5),]  
contour(x=e1, y=e0, z=matrix(nps2$et, nrow=length(e1), ncol=length(e0)), labcex=1.5, cex.lab=1.5, lty="dashed", xlab=expression(xi[1]), ylab=expression(xi[0]), levels=seq(-3,15,0.5), xlim=c(-3,10), ylim=c(0,10), method="edge")
polygon(x=c(-1,1),y=c(0,0), border="darkblue", lwd=2)

# xyplot(e1~e2|as.factor(e0), data=h0, group=et, type="l")
# wireframe(e1~e0*e2, data=nps, group=et, scales=list(arrows=FALSE), xlab=expression(eta[0]), ylab=expression(eta[2]), zlab=expression(eta[1]), zlim=c(-10,10), screen=list(x=-35,y=70,z=-25), auto.key=list(space="right"))

cumy <- cumsum(y)
ny <- seq_len(n)
xy_10 <- xy10 <- list()

for (i in 1:n){ # translation
	polygon(x=c(-1,1)+cumy[i], y=c(0,0)+ny[i], lwd=2, col="red3")
	text(x=0.2+cumy[i], y=0.2+ny[i], labels=bquote(y[.(i)]), cex=1.5)
	xy10[[i]] <- c(1,0) + c(cumy[i], ny[i])
	xy_10[[i]] <- c(-1,0) + c(cumy[i], ny[i])
}

# convh <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
# op <- iprior(ui=rbind(c(1,0,0), c(1,0,0), c(0,1,0), c(0,-1,0), c(0,0,1), c(0,0,-1)), ci=c(0.5, 0.5, -1, -1, 0, -0.5))
# scatterplot3d(op$vtx)
# op$vtx <- op$vtx*0.5 + 1
# update.impinf(op, y=numeric(0))

