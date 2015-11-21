# This demo describes level curves 
# foo2.R (in previous package)
# confirmed on 2015.11.21

# rm(list=ls())
library(imPois)
library(lattice)

# source("./../R/imPois.R")
# source("./../R/visualize.R")


# fails <- list()
e2 <- seq(from=0.05, to=4, by=0.25) # from=0.01, to=5, by=0.25 (if e2 <= 0.05, it fails. it is safe starting with 0.7)
e0 <- seq(from=0, to=2.5, by=0.5)   # from=0,    to=3, by=0.25
et <- -3:3
h0 <- expand.grid(e2=e2, e0=e0, et=et)
fn <- function(x){
	x <- as.vector(x)
	f0 <- function(xm, ...){ # xm = eta1
		val <- evfn(y=numeric(0), pars=c(x[1], xm, x[2]), ztrunc=FALSE)$value - x[3]
		return(val)
	}
	
	val <- tryCatch(uniroot(f0, lower=-10, upper=10, extendInt="yes", tol=1e-5)$root, error=function(e){ 
		# fails[[length(fails)+1]] <- x; 
		return(NaN) 
	})
# 	cat("[x=", x, "] val=", val, "\n")
	return(val)	
}
h0$e1 <- apply(h0, 1, fn)

wireframe(e1~e0*e2, data=h0, group=et, scales=list(arrows=FALSE), xlab=expression(eta[0]), ylab=expression(eta[2]), zlab=expression(eta[1]), zlim=c(-10,10), screen=list(x=-35,y=70,z=-25), auto.key=list(space="right"))

wireframe(e1~e0*e2, data=h0, group=et, scales=list(arrows=FALSE), xlab=expression(eta[0]), ylab=expression(eta[2]), zlab=expression(eta[1]), zlim=c(-10,10), screen=list(x=-40,y=-60,z=-25), auto.key=list(space="right"))

# opt1. screen=list(y=60, x=15, z=20), 
#	opt2. screen=list(y=105, x=15),
# opt3. screen=list(y=80, x=20, z=-20),
# xyplot(e1~e2|as.factor(e0), data=h0, group=et, type="l")

