#' @rdname cpm
#' @param m mean parameter
#' @export
# 
mpm <- function(m, xi2, xi1, xi0, log=FALSE){
	stopifnot(length(xi2)==1, length(xi1)==1, length(xi0)==1)
	val <- exp(-xi2*log(m)^2 + (xi1-1)*log(m) - xi0*m)
	if(log) val <- log(val)
	return(val)
}

#' @rdname cpm
#' @export 
# 
mpm1 <- function(m, xi2, xi1, xi0, log=FALSE){
	stopifnot(length(xi2)==1, length(xi1)==1, length(xi0)==1)
	val <- m*mpm(m, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)
	if(log) val <- log(val)
	return(val)
}


#' @rdname cpm
#' @export 
dmpm <- function(x, pars, log.p=FALSE){	

	stopifnot(length(pars) == 3, is.logical(log.p))
	
	p2 <- pars[1]
	p1 <- pars[2]
	p0 <- pars[3] 

	fx <- mpm(m=x, xi2=p2, xi1=p1, xi0=p0, log=FALSE)
	nconst <- stats::integrate(mpm, lower=0, upper=Inf, xi2=p2, xi1=p1, xi0=p0, log=FALSE)$value
	
	dx <- if (all(is.finite(fx), is.finite(nconst))) fx/nconst else 0
	if (log.p) dx <- log(dx) 

	return(dx)
}


#' @rdname cpm
#' @export
em.pdf <- function(y=NULL, pars){
	# TODO:
	# y=NULL -> y=numeric(0)
	
	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]
	
	fm <- function(m, ...) m*dmpm(x=m, pars=pars1, log.p=FALSE)
	em <- stats::integrate(fm, lower=0, upper=Inf)$value
	
	robj <- list(y=y, pars=pars, pars1=pars1, value=em)
	return(robj)
}

## approach 2. by using ratio of two integrals 
##

#' @rdname cpm
#' @export
em.rt <- function(y=NULL, pars){

	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]

	v0 <- stats::integrate(mpm, lower=0, upper=Inf, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)$value
	v1 <- stats::integrate(mpm1, lower=0, upper=Inf, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)$value
	et <- if (v0 != 0) v1/v0 else Inf

	robj <- list(y=y, pars=pars, pars1=pars1, value=et)
	return(robj)
}

