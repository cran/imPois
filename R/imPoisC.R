##
## ./R/imPoisC.R
##
## Freezed on 2017-03-05
##

#' @rdname cpm
#' @title Conjugate Prior Measure
#'
#' @description A conjugate prior measure on the canonical parameter is defined in the form of three-parameter exponenitial family of probability measure.  \code{cpm1} computes a normalizing constant.   \code{dcpm} and \code{pcpm} give density and distribution functions for the canonical variable \code{t}.  Also see the \sQuote{Details}.
#'
#' @param t variable
#' @param xi2 precision  
#' @param xi1 linear component
#' @param xi0 effective sample size
#' @param log logical
#'
#' @details  A formal definition of this conjugate prior measure is given by
#' \deqn{e^(-\xi_2\theta^2 + \xi_1\theta - \xi_0\exp(\theta))},
#' where \eqn{\theta} is a canonical variable ranged from \code{-Inf} to \code{Inf}, \eqn{\xi_2} is a precision parameter, \eqn{\xi_1} is a linear component, and \eqn{\xi_0} is an effective sample size.
#'  
#' @references
#' Lee, C.H. (2014) Imprecise Prior for Imprecise Inference on Poisson Sampling Models, PhD Thesis, University of Saskatchewan
#'
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export
cpm <- function(t, xi2, xi1, xi0, log=FALSE){
	stopifnot(length(xi2)==1, length(xi1)==1, length(xi0)==1)
	val <- exp(-xi2*t^2 + xi1*t - exp(log(xi0)+t))
	if(log) val <- log(val)
	return(val)
}

#' @rdname cpm
#' @export 
cpm1 <- function(t, xi2, xi1, xi0, log=FALSE){
	stopifnot(length(xi2)==1, length(xi1)==1, length(xi0)==1)
	val <- t*cpm(t, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)
	if(log) val <- log(val)
	return(val)
}

#' @rdname cpm 
#' @param x quantile	 
#' @param pars parameters of length 3
#' @param log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}
#' @export 
dcpm <- function(x, pars, log.p=FALSE){	

	stopifnot(length(pars) == 3, is.logical(log.p))
	
	p2 <- pars[1]
	p1 <- pars[2]
	p0 <- pars[3] 

	fx <- cpm(t=x, xi2=p2, xi1=p1, xi0=p0, log=FALSE)
	nconst <- stats::integrate(cpm, lower=-Inf, upper=Inf, xi2=p2, xi1=p1, xi0=p0, log=FALSE)$value
	
	dx <- if (all(is.finite(fx), is.finite(nconst))) fx/nconst else 0
	if (log.p) dx <- log(dx) 

	return(dx)
}

#' @rdname cpm
#' @param q quantiles
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X\le x]} otherwise, \eqn{P[X>x]}
#' @export
pcpm <- function(q, pars, lower.tail=TRUE, log.p=FALSE){
	
	stopifnot(length(q)==1, length(pars)==3, is.logical(lower.tail), is.logical(log.p))
	
	if (lower.tail) px <- stats::integrate(dcpm, lower=-Inf, upper=q, pars=pars)$value
	# TODO: else rv <- stats:integrate(dcpm, lower=q, upper=Inf, pars=pars)$value  
	if (log.p) px <- log(px)
	
	return(px)
}


#' @rdname cpm
#' @param y a vector of observations, by default, NULL
#' @export
et.pdf <- function(y=NULL, pars){
	# TODO:
	# y=NULL -> y=numeric(0)
	
	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]
	
	ft <- function(t, ...) t*dcpm(x=t, pars=pars1, log.p=FALSE)
	et <- stats::integrate(ft, lower=-Inf, upper=Inf)$value
	
	robj <- list(y=y, pars=pars, pars1=pars1, value=et)
	return(robj)
}


#' @rdname cpm
#' @export
ey <- function(y=NULL, pars){
	v1 <- stats::integrate(cpm, lower=-Inf, upper=Inf, xi2=pars[1], xi1=pars[2]+1, xi0=pars[3], log=FALSE)$value
	v0 <- stats::integrate(cpm, lower=-Inf, upper=Inf, xi2=pars[1], xi1=pars[2], xi0=pars[3], log=FALSE)$value
	rv <- v1/v0
	robj <- list(y=y, pars=pars, value=rv)
	return(robj) 
}


#' @rdname cpm
#' @export
et.rt <- function(y=NULL, pars){

	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]

	v0 <- stats::integrate(cpm, lower=-Inf, upper=Inf, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)$value
	v1 <- stats::integrate(cpm1, lower=-Inf, upper=Inf, xi2=xi2, xi1=xi1, xi0=xi0, log=FALSE)$value
	et <- if (v0 != 0) v1/v0 else Inf

	robj <- list(y=y, pars=pars, pars1=pars1, value=et)
	return(robj)
}



## approach 3. laplace approximation 
## NOTES: numerical results are very close to those found from approach 2.
## tolerance is 1e-8 when xi0=0, otherwise tolerance is varied from 1e-2 when xi0 is not 0.  

#' @rdname cpm
#' @param const constant
#' @export
et.la <- function(y=NULL, pars, const=1e3){
	
	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]
	
	# denominator 
	fn0 <- function(x, ...) cpm(t=x, xi2=xi2, xi1=xi1, xi0=xi0, log=TRUE)
	gr0 <- function(x, ...) -2*xi2*x + xi1 - xi0*exp(x)

	# numerator 
	fn1 <- function(x, ...) log(x+const) + cpm(t=x, xi2=xi2, xi1=xi1, xi0=xi0, log=TRUE)
	gr1 <- function(x, ...) 1/(x+const) - 2*xi2*x + xi1 - xi0*exp(x)

	fit0 <- stats::optim(par=1, fn=fn0, gr=gr0, method="L-BFGS-B", control=list(fnscale=-1), hessian=TRUE)
	fit1 <- stats::optim(par=1, fn=fn1, gr=gr1, method="L-BFGS-B", control=list(fnscale=-1), hessian=TRUE)
	
	sigma0 <- -1/fit0$hessian
	sigma1 <- -1/fit1$hessian

  et <- sqrt(sigma1/sigma0)*exp(fit1$value - fit0$value) - const
  robj <- list(y=y, pars=pars, pars1=pars1, value=et, fit0=fit0, fit1=fit1, const=const)
	
	return(robj)
}
NULL

## approach 4. metropolis-hastings algorithm 
## approach 5. importance sampling 


## E(theta) has been studied until now.  
## Since we have dcpm(x, pars, log.p=FALSE) now, we can find E(theta^2) and E(e^theta) as well. 
## theta^2 is associated with xi2
## e^theta is associated with xi0

#' @rdname cpm
#' @export
et2.pdf <- function(y=NULL, pars){

	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]
	
	ft <- function(t, ...) t^2*dcpm(x=t, pars=pars1, log.p=FALSE)
	et <- stats::integrate(ft, lower=-Inf, upper=Inf)$value
	
	robj <- list(y=y, pars=pars, pars1=pars1, value=et)
	return(robj)
}

#' @rdname cpm 
#' @export 
et0.pdf <- function(y=NULL, pars){

	stopifnot(length(pars) == 3)
	if (!is.null(y)) stopifnot(is.vector(y))
	
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")

	xi2 <- pars1[1]
	xi1 <- pars1[2]
	xi0 <- pars1[3]
	
	ft <- function(t, ...) exp(t)*dcpm(x=t, pars=pars1, log.p=FALSE)
	et <- stats::integrate(ft, lower=-Inf, upper=Inf)$value
	
	robj <- list(y=y, pars=pars, pars1=pars1, value=et)
	return(robj)
}

