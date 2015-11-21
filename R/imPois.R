#' @title Imprecise Inferential Framework for GLIM
#' 
#' @description The package 'ipeglim' (Imprecise Probability Estimates for GLIM)
#' is a collection of tools for conducting analysis of epistemic
#' uncertainty in the generalized linear model setup based on Walley's
#' (1991) Imprecise Probability Theory.
#'
#' @name ipeglim
#' @docType package
NULL

#' @title Bickis and Lee's Conjugate Formula for Imprecise Prior Measure
#' @description Bickis and Lee's conjugate formulation is described in the form of three-parameter exponential families.  See \sQuote{Details}.
#' @param t numeric value
#' @param xi2 parameter associated with precision 
#' @param xi1 parameter associated with linear combination
#' @param xi0 effective sample size
#' @param ztrunc a logical value incidcating whether standard Poisson is truncated at 0.
#'
#' @details  The formal definition of Bickis and Lee's conjugate formulation is
#' \deqn{exp(-\xi_2*\theta^2 + \xi_1*\theta + \exp(\theta))}
#' \eqn{t} is ranged from \code{-Inf} to \code{Inf}. 
#'
#' @references
#' Lee, C.H. (2014) Imprecise Prior for Imprecise Inference on Poisson Sampling Model, PhD Thesis, The Collaborative Biostatistics Program, School of Public Health, University of Saskatchewan
#' 
# @return
#' 
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export
kcpm <- function(t, xi2, xi1, xi0, ztrunc){
	if (ztrunc) {
		ln.1mp0 <- stats::ppois(q=0, lambda=exp(t), lower.tail=FALSE, log.p=TRUE)
		ft <- exp(-xi2*t^2 + xi1*t - exp(log(xi0)+t) - xi0*ln.1mp0)
	} else ft <- exp(-xi2*t^2 + xi1*t - exp(log(xi0)+t))
	return(ft)
}


#' @title Comupting Normalizing Constant of Bickis and Lee's Probability Distribution
#' @description Given parameters, a normalizing constant for Bickis and Lee's probability distribution is computed. 
#' @param xi2 parameter associated with precision 
#' @param xi1 parameter associated with linear combination
#' @param xi0 effective sample size
#' @param ztrunc a logical value incidcating whether standard Poisson is truncated at 0.
#' @param log a logical value; if TRUE, probabilities are given as \eqn{log(p)}.
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export 
cgf <- function(xi2, xi1, xi0, ztrunc, log=TRUE){
	
	# backward approach - from Inf to 0 (* by theory)
	op <- obj <- list()
	# obj$value <- Inf 	
	cutoff <- c(Inf, 1e9, 1e8, 1e7, 5e6, 1e6, 5e5, 1e5, 5e4, 1e4, 5e3, 1e3, 5e2, 3e2, 1e2, 5e1, 3e1, 2e1)
	
	# this integration should have the same value regardless of canonical parametrization or mean parameterization
	for(i in 1:length(cutoff)){
	
		limit <- cutoff[i]
		op[[i]] <- tryCatch(stats::integrate(kcpm, lower=-limit, upper=limit, xi2=xi2, xi1=xi1, xi0=xi0, ztrunc=ztrunc), 
			error=function(e){
				errmsg <- e$message
				if (grepl("non-finite function", errmsg)) e$value <- e$abs.error <- Inf # in fact, value is NULL
				else print(errmsg)
				return(e)
			})
			
		if(is.finite(op[[1]]$value)){
			# this value must be positive
			value <- op[[1]]$value
			value <- log(value)
			return(list(attempts=op, limit=limit, value=value, i=i))
		}
			
		if (i > 2) {
			v_2 <- op[[i-2]]$value
			v_1 <- op[[i-1]]$value
			v0 <- op[[i]]$value
			err01 <- abs(v0-v_1)
			err02 <- abs(v0-v_2)
			err12 <- abs(v_1-v_2)
			if(is.na(err01)) next
			if ( (err01<1e-8) & (err02<1e-8)) break
		}
	#	print(i)
	} # end of for
	
	vals <- do.call(c,lapply(op, "[[", "value"))
	value <- log(vals[i])
	names(vals) <- 1:i
	
	limit <- cutoff[1:i]
	names(limit) <- 1:i
	
	robj <- list(vals=vals, attempts=op, limit=limit, value=value, i=(i-2))
	return(robj)
	# if integration fails, NULL will be returned.
}

#' @rdname BL.dist
#' @title BL Probability Distribution
#' @description Density and distribution function for the Bickis and Lee's distribution with three parameters \eqn{xi_2}, \eqn{xi_1}, and \eqn{xi_0}. 
#' @param x quantiles
#' @param pars a numeric vector of parameters
#' @param ztrunc a logical value incidcating whether standard Poisson is truncated at 0.
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export 
dcpm <- function(x, pars, ztrunc){
	p2 <- pars[1]
	p1 <- pars[2]
	p0 <- pars[3]
	fx <- kcpm(t=x, xi2=p2, xi1=p1, xi0=p0, ztrunc=ztrunc)
	nconst <- cgf(xi2=p2, xi1=p1, xi0=p0, ztrunc=ztrunc)$value
	p <- fx/exp(nconst)
	return(p)
}

#' @rdname BL.dist
#' @param q quantiles
#' @export
pcpm <- function(q, pars, ztrunc){
	p2 <- pars[1]
	p1 <- pars[2]
	p0 <- pars[3]
	op <- cgf(xi2=p2, xi1=p1, xi0=p0, ztrunc=ztrunc)
	limit <- op$limit[op$i]
	rv <- stats::integrate(dcpm, lower=-op$limit[op$i], upper=q, pars=pars, ztrunc=ztrunc)$value
	return(rv)
}

#' @title Expected Value of Canonical Variable
#' @description Expected value of canonical variable is computed using \code{integrate} function.
#' @param y a vector of observations
#' @param pars a numeric vector of parameters 
#' @param ztrunc a logical value incidcating whether standard Poisson is truncated at 0.
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca} 
#' @export
evfn <- function(y, pars, ztrunc){
	xi2 <- pars[1]
	xi1 <- pars[2]
	xi0 <- pars[3]
	pars1 <- c(xi2, xi1+sum(y), xi0+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")
	
	op <- cgf(xi2=pars1[1], xi1=pars1[2], xi0=pars1[3], ztrunc=ztrunc)
	limit <- op$limit[op$i]
	ft <- function(t, ...) t*dcpm(t, pars=pars1, ztrunc=ztrunc)
	value <- stats::integrate(ft, lower=-limit, upper=limit, ztrunc=ztrunc)$value
	
	robj <- list(pars=pars, pars1=pars1, attempts=op, value=value)
	invisible(robj)
}
