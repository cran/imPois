#' @title Characterize Imprecise Prior
#' @description A set of linear inequalities is used for characterizing a polygonal convex set.  The \code{chull} function in the \code{grDevices} package and the \code{convhulln} function in the \code{geometry} package are used to search for extreme points constructing a convex hull. 
#' @param ui constraint matrix (k x p), see below.
#' @param ci constrain vector of length k, see below.
#' @examples 
#' # lc0 <- list(lhs=rbind(diag(2), -diag(2)), rhs=c(0,0,-1,-1))
#' # op <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
#' # op <- iprior(ui=rbind(c(1,0),c(0,1),c(-1,-1)), ci=c(0,0,-5)) 
#' op <- iprior(ui=rbind(c(1,0),c(0,1),c(0,-1),c(1,1),c(-2,-1)), 
#'    ci=c(1,2,-8,5,-14)) # (3,8),(1,8), (1,4),(3,2)(6,2)
#' 
#' @author Chel Hee Lee <gnustats@@gmail.com>
#' @export 
iprior <- function(ui, ci){ # in 2 dimensions
	
	allsys <- utils::combn(nrow(ui), ncol(ui))
	sol <- list()
	for (i in 1:ncol(allsys)) {
		idx <- c(allsys[,i])
		A <- ui[idx,]
		b <- ci[idx]
		sol[[i]] <- tryCatch(solve(A,b), error=function(e) return(NaN))
	}
	roots <- do.call(cbind, sol)
	idx <- which(colSums(ui %*% roots >= ci) == nrow(ui))
	vtx <- t(roots[,idx])
	
	tmp <- vtx
	
	if (ncol(vtx) == 2) {
		vtx <- vtx[grDevices::chull(vtx),]
	}
	
	if (ncol(vtx) >= 3) {
		# convex polyhedron - http://mathworld.wolfram.com/ConvexPolyhedron.html
		imat <- t(geometry::convhulln(vtx, options="Tv"))
		vtx <- vtx[imat,]
		vtx <- vtx[!duplicated(vtx),]
	}
	
	rownames(vtx) <- paste("p", 1:nrow(vtx), sep="")
	colnames(vtx) <- paste("d", 1:ncol(vtx), sep="")
	
	robj <- list(ui=ui, ci=ci, sol=sol, roots=roots, idx=idx, vtx=vtx, tmp=tmp)
	class(robj) <- "impinf"
	return(robj)
}


#' @title Applying Bayes Rule
#' @param object an object for which the Bayes rule is needed. 
#' @param y a vector of observations
#' @param ... further arguments passed to methods
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export 
update.impinf <- function(object, y, ...){ 
	
	stopifnot(inherits(object, "impinf"))
	
	vtx <- object$vtx # matrix
	
	op <- apply(vtx, 1, function(x){
		if (length(x) == 2) tryCatch(evfn(y=y, pars=c(0,x[1],x[2]), ztrunc=FALSE)$value, error=function(e) return(NaN))
		else tryCatch(evfn(y=y, pars=c(x[1],x[2],x[3]), ztrunc=FALSE)$value, error=function(e) return(NaN))
	})
	
	if (ncol(vtx) == 2) vtx1 <- cbind(vtx[,1]+sum(y), vtx[,2]+length(y))
	else vtx1 <- cbind(vtx[,1], vtx[,2]+sum(y), vtx[,3]+length(y))
	
	robj <- list(vtx=vtx, y=y, vtx1=vtx1, ev=op)
	class(robj) <- "impinf"
	return(robj)
}

#' @title Summary of \code{impinf} object
#' @param object an object for which a summary is needed.
#' @param ... additional arguments affecting the summary produced.
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export
summary.impinf <- function(object, ...){
	stopifnot(inherits(object, "impinf"))
	
	vtx <- object$vtx
	vtx1 <- object$vtx1
	ev <- object$ev
	
	inf <- ev[which.min(ev)]
	sup <- ev[which.max(ev)]
	delta <- sup - inf
	names(delta) <- NULL
	
	inf.p <- vtx[which.min(ev), ]
	sup.p <- vtx[which.max(ev), ]
	inf.p1 <- vtx1[which.min(ev), ]
	sup.p1 <- vtx1[which.max(ev), ]
	
	robj <- list(inf=inf, sup=sup, delta=delta, inf.p=inf.p, sup.p=sup.p, inf.p1=inf.p1, sup.p1=sup.p1)
	class(robj) <- "impinf"
	return(robj)
}
