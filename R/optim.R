#' @rdname hparaOptim
#' @title Hyperparamter Optimization
#' @description Score and graident functions are defined
#' @param x parameter needs to be optimized
#' @export 
fn.evfn <- function(x){
	if (length(x) == 2) val <- et.pdf(y=numeric(0), pars=c(0,x[1],x[2]))$value
	if (length(x) == 3) val <- et.pdf(y=numeric(0), pars=c(x[1],x[2],x[3]))$value
	return(val)
}

#' @rdname hparaOptim
#' @export 
gr.evfn <- function(x){
	len.x <- length(x)

	# E[theta]	
	em1 <- function(p2,p1,p0){
		if(len.x == 2) p2 <- 0 
		f0 <- function(x) x*dcpm(x=x, pars=c(p2,p1,p0))
		val <- stats::integrate(f=f0, lower=-Inf, upper=Inf)$value
		return(val)
	}
	
	# E[theta^2]
	em2 <- function(p2,p1,p0){
		if (len.x == 2) p2 <- 0
		f0 <- function(x) x^2*dcpm(x=x, pars=c(p2,p1,p0))
		val <- stats::integrate(f=f0, lower=-Inf, upper=Inf)$value
		return(val)
	}
	
	# E[theta^3]
	em3 <- function(p2,p1,p0){
		if (len.x == 2) p2 <- 0
		f0 <- function(x) x^3*dcpm(x=x, pars=c(p2,p1,p0))
		val <- stats::integrate(f=f0, lower=-Inf, upper=Inf)$value
		return(val)
	}
	
	# dE[theta]/d eta_2
	deriv2.em <- function(p2,p1,p0){
		val <- em1(p2=p2,p1=p1,p0)*em2(p2=p2,p1=p1,p0=p0)-em3(p2=p2,p1=p1,p0=p0)
		return(val)
	}	
	
	# dE[theta]/d eta_1
	deriv1.em <- function(p2,p1,p0){
		v1 <- em1(p2=p2,p1=p1,p0=p0)
		val <- em2(p2=p2,p1=p1,p0=p0)-v1*v1
		return(val)
	}
	
	# dE[theta]/d eta_0
	deriv0.em <- function(p2,p1,p0){
		# pred.y <- cgf(xi2=p2,xi1=(p1+1),xi0=p0,log=FALSE)$value/cgf(xi2=p2,xi1=p1,xi0=p0,log=FALSE)$value
		pred.y <- ey(y=NULL, pars=c(p2,p1+1,p0))$value/ey(y=NULL, pars=c(p2,p1,p0))$value
		v1 <- em1(p2=p2,p1=p1,p0=p0)
		v1_1 <- em1(p2=p2,p1=(p1+1),p0=p0)
		val <- pred.y*(v1-v1_1)
		return(val)
	}
	
	if (len.x == 2) {
		x1 <- x[1]
		x0 <- x[2]
		val <- c(deriv1.em(0,p1=x1,p0=x0), deriv0.em(0,p1=x1,p0=x0))
	}

	if (len.x == 3) {
		x2 <- x[1]
		x1 <- x[2]
		x0 <- x[3]
		val <- c(deriv2.em(p2=x2,p1=x1,p0=x0), deriv1.em(p2=x2,p1=x1,p0=x0), deriv0.em(p2=x2,p1=x1,p0=x0)) 
	}
	
	return(val)
}

