#' @rdname bary2cart
#' @title Conversion from bary to cartesian
#' @description Conversion from barycentric coordinates to cartesian coordinates
#' @param x point coordinate information
#' @export 
#
bary2cart <- function(x){	
	
	# default cartesian coordinates
	r1 <- c(0,0)
	r2 <- c(1,0)
	r3 <- c(1/2,sqrt(3)/2)
	
	R <- rbind(cbind(r1,r2,r3),1)
	r <- R%*%x
	
	return(r)
}

#' @rdname bary2cart
#' @export
bary2cart3d <- function(x){	
	
	# default cartesian coordinates
	r1 <- c(0,sqrt(2/3),0)
	r2 <- c(1/sqrt(3),0,0)
	r3 <- c(-1/(2*sqrt(3)),0,-1/2)
	r4 <- c(-1/(2*sqrt(3)),0,1/2)
	
	R <- rbind(cbind(r1,r2,r3,r4),1)
	r <- R%*%x
	
	return(r)
}


# generating random variates from dirichlet distribution
# rdirichlet <- function(n, a){
# 	op <- lapply(1:n, function(x){
# 		y <- rgamma(n=length(a), shape=a, rate=1)
# 		return(y/sum(y))
# 	})
# 	op <- do.call(rbind, op)
# 	return(op)
# }
