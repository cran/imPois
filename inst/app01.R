##
## freezed on 2017-02-28
## 
## Resample 10 times from the existing dataset (Example 1, Howlader 2003)
## 
## Howlader (2003) Bayesian estimation of the distribution function of the 
## Poisson model 
## 

rm(list=ls())
source("./../R/imTools.R")
source("./../R/imPoisC.R")
source("./../R/imPoisM.R")
source("./../R/imGraphics.R")
# source("./../R/imZTPois.R")

rm(list=ls())
library(imPois)


x <- c(0,1,2,3,4,5,6,7,8,9)
nx <- c(75,103,121,54,30,13,2,1,0,1)  # sum(nx) = 400, sum(nx*x)=720
y <- rep(x, times=nx)

#
pdf(file="app1_resampling_10times.pdf")
#

xi1 <- seq(from=0, to=900, by=0.5)
xi0 <- seq(from=0, to=500, by=0.5)
hs <- expand.grid(x1=xi1, x0=xi0)
hs$f <- apply(hs, 1, function(x) x[1]/x[2])
hsf <- matrix(hs$f, nrow=length(xi1))

contour(x=xi1, y=xi0, z=hsf, method="edge", levels=seq(1,2.4,0.2), labcex=1.25, lwd=2)
title(xlab=expression(xi[1]~"(sum)"), cex.lab=1.25, line=2)
title(ylab=expression(xi[0]~"(number of observations)"), cex.lab=1.25, line=2)

cvh <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
set.seed(12345678)

for(j in 1:10){
	
	y <- sample(y)
	
	for (i in 1:length(y)) {
		
		yi <- y[1:i]
		vtx1 <- cvh$vtx
		vtx1[,1] <- vtx1[,1] + sum(yi)
		vtx1[,2] <- vtx1[,2] + length(yi)
		polygon(vtx1, border="grey", lwd=1.5)
		
		# keep the following codes
		# op <- update(object=cvh, y=y[1:i], wrt="mean")
		# polygon(op$vtx1, border="grey", lwd=1.5)
		# vtx1.m <- colMeans(op$vtx1)
		
	}
	
}

polygon(x=c(0,0,1,1)+sum(nx*x),y=c(0,1,1,0)+sum(nx), border="darkblue", lwd=2, cex=3)
points(x=720,y=400,col="darkblue",pch=17)

dev.off()



###
### The codes from this points have been used for exercise 
### 2017-02-28

##
## Approach 1. contour lines are found by using "uniroot()" function
##
rm(list=ls())
xi1 <- seq(from=620, to=820, by=1) 
mlvls <- seq(1,2,0.2)
hs <- expand.grid(xi1=xi1, em=mlvls)
fn <- function(x){
	f0 <- function(xm, ...) em.pdf(y=numeric(0), pars=c(0,x[1],xm))$value-x[2]
	val <- tryCatch(uniroot(f0, lower=200, upper=600, extendInt="yes", tol=1e-5)$root, error=function(e) return(NaN))
	return(val)
} 
hs$xi0 <- apply(hs, 1, fn) 
plot(x=hs$xi1, y=hs$xi0, type="n")
for(v in mlvls){
	hsi <- subset(hs, subset=(em==v))
	# lines(x=hsi$xi1, hsi$xi0, col="grey", lty="dashed", lwd=1.5)
	points(x=hsi$xi1, y=hsi$xi0, cex=0.3, pch=19)
}

##
## Approach 2. contour lines are found by using "contour()" function 
##
rm(list=ls())
xi1 <- seq(from=620, to=820, by=1) # T=720
xi0 <- seq(from=300, to=500, by=1) # n=400
hs <-expand.grid(x1=xi1, x0=xi0)
fn <- function(x){
	rv <- tryCatch(em.pdf(y=numeric(0), pars=c(0,x[1],x[2]))$value, error=function(e) NA)
	# rv <- tryCatch(et.rt(y=numeric(0), pars=c(0,x[1],x[2]))$value, error=function(e) NA)
	return(rv)
}
hs$f <- apply(hs,1,fn)
hsf <- matrix(hs$f, nrow=length(xi1))
contour(x=xi0, y=xi1, z=hsf)
points(x=400,y=720)

