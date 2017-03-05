## zero-truncated Poisson 
## application 

rm(list=ls())
source("./../R/imTools.R")

mpm.ztrunc <- function(m, xi2, xi1, xi0, log=FALSE){
	stopifnot(length(xi2)==1, length(xi1)==1, length(xi0)==1)
	# val <- exp(-xi2*log(m)^2 + (xi1-1)*log(m) - xi0*m) # -> 1.561661 
	# val <- exp(-xi2*log(m)^2 + (xi1-1)*log(m) - xi0*log(exp(m)-1)) # -> 0.9691044
	val <- exp(-xi2*log(m)^2 + (xi1-1)*log(m) - xi0*log(expm1(m))) # -> 0.9691044
	# val <- exp(-xi2*log(m)^2 + (xi1-1)*log(m) - xi0*(m+ppois(q=0, lambda=m, lower.tail=FALSE, log.p=TRUE))) # -> 0.9691044
	if(log) val <- log(val)
	return(val)
}

#' @rdname cpm
#' @export 
# 
mpm1.ztrunc <- function(m, xi2, xi1, xi0, log=FALSE){
	stopifnot(length(xi2)==1, length(xi1)==1, length(xi0)==1)
	val <- m*mpm.ztrunc(m, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)
	if(log) val <- log(val)
	return(val)
}


sumy <- 32*1+2*16+3*6+4*1
n <- 32+16+6+1 
rv0 <- integrate(mpm.ztrunc, lower=0, upper=Inf, xi2=0, xi1=sumy, xi0=n)$value
rv1 <- integrate(mpm1.ztrunc, lower=0, upper=Inf, xi2=0, xi=sumy, xi0=n)$value
rv1/rv0




#' @rdname ecpm
#' @export 
dmpm.ztrunc <- function(x, pars, log.p=FALSE){	

	stopifnot(length(pars) == 3, is.logical(log.p))
	
	p2 <- pars[1]
	p1 <- pars[2]
	p0 <- pars[3] 

	fx <- mpm.ztrunc(m=x, xi2=p2, xi1=p1, xi0=p0, log=FALSE)
	nconst <- stats::integrate(mpm.ztrunc, lower=0, upper=Inf, xi2=p2, xi1=p1, xi0=p0, log=FALSE)$value
	
	dx <- if (all(is.finite(fx), is.finite(nconst))) fx/nconst else 0
	if (log.p) dx <- log(dx) 

	return(dx)
}


#' @rdname ecpm
#' @export
em.pdf.ztrunc <- function(y=NULL, pars){
	# TODO:
	# y=NULL -> y=numeric(0)
	
	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]
	
	fm <- function(m, ...) m*dmpm.ztrunc(x=m, pars=pars1, log.p=FALSE)
	em <- stats::integrate(fm, lower=0, upper=Inf)$value
	
	robj <- list(y=y, pars=pars, pars1=pars1, value=em)
	return(robj)
}

sumy <- 32*1+2*16+3*6+4*1
n <- 32+16+6+1 
em.pdf.ztrunc(y=numeric(0), pars=c(0,sumy,n))   # 0.969102



## approach 2. by using ratio of two integrals 
##

#' @rdname ecpm
#' @export
em.rt.ztrunc <- function(y=NULL, pars){

	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]

	v0 <- stats::integrate(mpm.ztrunc, lower=1e-2, upper=Inf, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)$value
	v1 <- stats::integrate(mpm1.ztrunc, lower=1e-2, upper=Inf, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)$value
	et <- if (v0 != 0) v1/v0 else Inf

	robj <- list(y=y, pars=pars, pars1=pars1, value=et)
	return(robj)
}


sumy <- 32*1+2*16+3*6+4*1
n <- 32+16+6+1 
em.rt.ztrunc(y=numeric(0), pars=c(0,sumy,n))




#### 
rm(list=ls())
source("./../R/imZTPoisM.R")
source("./../R/imTools.R")

#####
rm(list=ls())
library(imPois)

pdf(file="app_cholera_ztruncpois_EY_lines.pdf")

##
## using contour plot
## 
sumy <- 32*1+2*16+3*6+4*1
n <- 32+16+6+1 

xi1 <- seq(from=0, to=100, by=0.5) # alpha: shape 	
xi0 <- seq(from=0, to=60, by=0.5) # beta : rate 
hs <- expand.grid(xi1=xi1, xi0=xi0)

hs$em <- apply(hs, 1, function(x){
	rv <- tryCatch(em.rt.ztrunc(y=numeric(0), pars=c(0,x[1],x[2]))$value, error=function(e) return(NaN))
	return(rv)
})

contour(x=xi1, y=xi0, z=matrix(hs$em, nrow=length(xi1), ncol=length(xi0)), levels=seq(0, 1.5, 0.5), labcex=1.5, xlim=c(0,100), ylim=c(0,60), method="flattest")
title(xlab=expression(xi[1]~"(sum)"), cex.lab=1.25, line=2)
title(ylab=expression(xi[0]~"(number of observations)"), cex.lab=1.25, line=2)
polygon(x=c(0,0,1,1)*2,y=c(0,1,1,0)*2, border="darkred", col="darkred") # imprecise prior 
polygon(x=c(0,0,1,1)*2+sumy,y=c(0,1,1,0)*2+n, border="darkred", col="darkred") # imprecise prior (updated) 


####
#### more detailed (using uniroot function)
#### 

xi1a <- seq(from=0, to=100, by=1) # alpha: shape
mlvls <- seq(0.7, 1.3, by=0.2)
hs1 <- expand.grid(xi1a=xi1a, em=mlvls)

fn <- function(x){
	f0 <- function(xm, ...) em.rt.ztrunc(y=numeric(0), pars=c(0,x[1],xm))$value-x[2]
	val <- tryCatch(uniroot(f0, lower=1e-2, upper=100)$root, error=function(e) return(NaN))
	return(val)
} 
hs1$xi0a <- apply(hs1, 1, fn) 

for(v in mlvls){
	hsi <- subset(hs1, subset=(em==v))
	lines(x=hsi$xi1a, hsi$xi0a, col="grey", lty="dashed", lwd=1.5)
}
text(x=c(80,90,100,100), y=c(60,60,60,55), labels=c("0.7", "0.9", "1.1", "1.3"))

dev.off()

em.rt.ztrunc(y=numeric(0), pars=c(0,0+sumy,1+n))$value # minimum from (0,1)
em.rt.ztrunc(y=numeric(0), pars=c(0,1+sumy,0+n))$value # maximum from (1,0)




## maximum likelihood estimate from Irwin (1958)
x <- c(1,2,3,4)
nx <- c(32,16,6,1)
xbar <- sum(nx*x)/sum(nx)

j <- 0
rv <- 0
epsilon <- 1
while(epsilon > 1e-8){
	j <- j+1
	rv1 <- rv + j^(j-1)/factorial(j)*(xbar*exp(-xbar))^j
	epsilon <- abs(rv1-rv)
	rv <- rv1
}
xbar - rv

# McKendrick (1925), Applications of mathematics to medical problems
# Dahiya and Gross (1973), Estimating the zero class from a truncated Poisson sample -> lambda = 0.970 
# Blumenthal et al. (1978), Estimating the complete sample size from an incomplete Poisson sample. -> not clear 
# Scollnik (1997) Inference concerning the size of the zero class from an incompelte poisson sample -- E(lambda)=0.993, 1.028, 0.982, 1.023 from various prior specification 
# Howlader (2003), lambda=0.9721 (Bayes)

