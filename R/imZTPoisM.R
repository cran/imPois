#' @rdname cpm
#' @export 
# 
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




#' @rdname cpm
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


#' @rdname cpm
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

#' @rdname cpm
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


