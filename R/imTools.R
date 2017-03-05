#' @rdname iprior
#' @title Characterizing Imprecise Prior
#' @description A set of prior distributions is characterized as an imprecise prior for inference.  See \sQuote{Details}. 
#' @details A convex set of hyperparameters is charcterized using a set of linear inequalities.   \code{grDevices::chull} and \code{geometry::convhulln} functions are used to search for extreme points of this convex set. 
#'
#' @param ui constraint matrix (k x p), see below.
#' @param ci constrain vector of length k, see below.
#' @param pmat matrix (k x p) containig coordinate information in d-dimensions.
#'
#' @examples
#' ## 2-dims (xi2=0, xi1, xi0)
#' lc0 <- list(lhs=rbind(diag(2), -diag(2)), rhs=c(0,0,-1,-1))
#' op <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
#' op <- iprior(ui=rbind(c(1,0),c(0,1),c(-1,-1)), ci=c(0,0,-5)) 
#' op <- iprior(ui=rbind(c(1,0),c(0,1),c(0,-1),c(1,1),c(-2,-1)), 
#'              ci=c(1,2,-8,5,-14)) # (3,8),(1,8), (1,4),(3,2)(6,2)
#'
#' ## 3-dimes (xi2, xi1, xi0)
#' op <- iprior(ui=rbind(c(1,0,0), c(-1,0,0), c(0,1,0), c(0,-1,0), c(0,0,1)), 
#'              ci=c(-0.5, -1, -2, -2, 0))
#' op <- iprior(ui=rbind(c(1,0), c(-1,0), c(0,1), c(0,-1)), 
#'              ci=c(0.5, -1, -2, -2))
#' lc0 <- cbind(rbind(c(1,0,0), c(-1,0,0), c(0,1,0), c(0,-1,0), c(0,0,1), 
#'              c(0,0,-1)), c(0.5, -1, -2, -2,0,-1))
#' iprior(pmat=lc0)
#' lc0 <- rbind(c(-2,1,0), c(2,1,0), c(-2,0.5,0), c(2,0.5,0))
#' lc0 <- rbind(c(1,2,0), c(1,-2,0), c(0.5,-2,0), c(0.5,2,0)) 
#' iprior(pmat=lc0)
#' 
#' @export 
iprior <- function(ui, ci, pmat){ # in 2 dimensions
	
	if (missing(ui)) ui <- NULL
	if (missing(ci)) ci <- NULL
	if (missing(pmat)) pmat <- NULL
	
	robj <- list()
	
	# TODO:
	# Direct input of coordinate information
	# The input does not need to have a convex set
	# pmat must be in matrix form
	# 
	if (!is.null(pmat) & (is.null(ui) | is.null(ci))) {
		if (!is.null(ui) | !is.null(ci)) message("ui or ci is ignored; instead, uses NA.")

		if(FALSE){		
		# TODO:
		# if (ncol(pmat) == 1) stop("Matrix with at least two columns is needed.")
		# if (ncol(pmat) == 2) {
		#	imat <- grDevices::chull(pmat)
		#	vtx <- pmat[imat,]
		#	 # columns should be arranged in (xi1, xi0)
		# }
		# if (ncol(pmat) >= 3) {
		# 	robj$vtx0 <- pmat
		# 	imat <- t(geometry::convhulln(pmat, options="Tv"))
		# 	robj$imat <- imat
		# 	vtx <- pmat[imat,]
		# 	vtx <- vtx[!duplicated(vtx),]
		# 	# columns should be arranged in (xi2, xi1, xi0)
		# }
		} # end of if(FALSE)
		vtx <- pmat
	}
	
	# working with a set of linear inequalities 
	if(!is.null(ui) & !is.null(ci)){
		if( !is.null(pmat) ) message("pmat is ignored.")
		pmat <- NULL
	
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
		
		if (ncol(vtx) == 2) {
			vtx <- vtx[grDevices::chull(vtx),]
		}

		if (ncol(vtx) == 3) {
			vtx <- vtx[!duplicated(vtx),]
		}
		
		if (ncol(vtx) >= 4) {
			imat <- t(geometry::convhulln(vtx, options="Tv"))
			robj$imat <- imat
			robj$vtx0 <- vtx
			vtx <- vtx[imat,]
			vtx <- vtx[!duplicated(vtx),]
		}
		robj$ui <- ui
		robj$ci <- ci
		robj$sol <- sol
		robj$roots <- roots
		robj$idx <- idx 
	} # end of if ui and ci
	
	rownames(vtx) <- paste("p", 1:nrow(vtx), sep="")
	colnames(vtx) <- paste("d", 1:ncol(vtx), sep="")
	
	robj$vtx <- vtx
	class(robj) <- c("iprior","impinf")
	return(robj)
}



#' @rdname update
#' @title Applying Bayes Rule
#' @description The Bayes rule is applied to an imprecise prior and produce an imprecise posterior.
#' @param object an object for which an update is needed
#' @param y vector of observations
#' @param wrt parameterization method with respect to canonical or mean
#' @param ... further arguments passed to methods
#'
#' @examples
#'
#' # 2-dimensions
#' lc0 <- list(lhs=rbind(diag(2), -diag(2)), rhs=c(0,0,-1,-1))
#' op <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
#' op <- iprior(ui=rbind(c(1,0),c(0,1),c(-1,-1)), ci=c(0,0,-5)) 
#' op <- iprior(ui=rbind(c(1,0),c(0,1),c(0,-1),c(1,1),c(-2,-1)), 
#'              ci=c(1,2,-8,5,-14)) # (3,8),(1,8), (1,4),(3,2)(6,2)
#' op1 <- update(op, y=NULL)
#'
#' # 3-dimensions
#' lc0 <- rbind(c(1,2,0), c(1,-2,0), c(0.5,-2,0), c(0.5,2,0)) 
#' op <- iprior(pmat=lc0)
#' op1 <- update(op, y=NULL)
#'
#' @export 
update.impinf <- function(object, y=NULL, wrt=c("canonical", "mean"), ...){ 
	
	## TODO: 'wrt' argument for entire program 
	## TODO: add some argument to clearly show which prior is used. (2 dim = lgamma, 3 dim=normal)
	## 
	stopifnot(inherits(object, "impinf"))
	
	vtx <- as.matrix(object$vtx) # matrix
	# TODO: colnames(vtx1) <- colnames(vtx)

	op <- apply(vtx, 1, function(x){
	# TODO: optimize code here in terms of 'pars' and 'pars1'
		
		if (length(x) == 2) {  # TODO: iprior <- "lgamma"

			if (wrt == "canonical") {
			
				tryCatch(et.pdf(y=y, pars=c(0, x[1],x[2])), 
					error=function(e) {
						robj <- list(y=y, pars=c(0, x[1],x[2]), pars1=c(0, x[1]+sum(y),x[2]+length(y)), value=NA)
						return(robj)
						}) 

			} else if (wrt == "mean") { 
				
				tryCatch(em.pdf(y=y, pars=c(0, x[1],x[2])), 
					error=function(e) {
			 			robj <- list(y=y, pars=c(0, x[1],x[2]), pars1=c(0, x[1]+sum(y),x[2]+length(y)), value=NA)
						return(robj)
			 		}) 

			}

		} else if (length(x) == 3) {  # TODO: iprior <- "normal"

			if (wrt == "canonical") { 

				tryCatch(et.pdf(y=y, pars=c(x[1],x[2],x[3])), 
					error=function(e) {
						robj <- list(y=y, pars=c(x[1],x[2],x[3]), pars1=c(x[1],x[2]+sum(y),x[3]+length(y)), value=NA)
						return(robj)
					} ) 

			} else if (wrt == "mean") { 
				tryCatch(em.pdf(y=y, pars=c(x[1],x[2],x[3])), 
					error=function(e) {
						robj <- list(y=y, pars=c(x[1],x[2],x[3]), pars1=c(x[1],x[2]+sum(y),x[3]+length(y)), value=NA)
					return(robj)
					} ) 
	
			}

		}
	
	}) 

	vtx1 <- do.call(rbind, lapply(op, "[[", "pars1"))
	if (ncol(vtx) == 2) vtx1 <- vtx1[,-1]
	# print(vtx1)

	et <- do.call(c, lapply(op, "[[", "value"))

	robj <- list(y=y, vtx=vtx, vtx1=vtx1, op=op, et=et)
	class(robj) <- c("update", "impinf")
	return(robj)
}



#' @rdname summary
#' @title Summaries of \code{impinf} object
#' @description Summarizing outputs produced from \code{update.impinf}. 
#' @param object an object for which a summary is needed.
#' @param ... additional arguments affecting the summary produced.
#' @export
summary.impinf <- function(object, ...){
	stopifnot(inherits(object, "impinf"))
	
	y <- object$y
	vtx <- object$vtx
	vtx1 <- object$vtx1
	et <- object$et
	em <- object$em
	if(!is.null(et)) ev <- et
	else ev <- em
	
	# infimum and supremum 
	inf <- ev[which.min(ev)]
	sup <- ev[which.max(ev)]
	delta <- sup - inf
	names(delta) <- NULL
	
	# coordinates
	inf.p <- vtx[which.min(ev), ]
	sup.p <- vtx[which.max(ev), ]
	inf.p1 <- vtx1[which.min(ev), ]
	sup.p1 <- vtx1[which.max(ev), ]
	
	robj <- list(y=y, inf=inf, sup=sup, delta=delta, inf.p=inf.p, sup.p=sup.p, inf.p1=inf.p1, sup.p1=sup.p1, ev=ev)
	class(robj) <- c("summary.impinf")
	invisible(robj)
}


#' @rdname print
#' @title Print Summarized Objects
#' @description Printing \code{impinf} objects
#' @param x an object used to select a method 
#' @param ... further arguments passed to or from other methods 
#' @export
print.summary.impinf <- function(x, ...){
	# stopifnot(inherits(x, "impinf"), inherits(x, "summary"))
	stopifnot(inherits(x, "summary.impinf"))
	
	m <- c(x$inf, x$sup, x$delta)	
	names(m) <- c("min.E(X|y)", "max.E(X|y)", "imprecision")
	print(m)
	
	cat("\n  min.E(X|y) and max.E(X|y) are found \nat", names(x$inf), "and", names(x$sup), "respectively.\n\n  Coordinates:\n" )
	m <- rbind(x$inf.p, x$sup.p)
	rownames(m) <- c(names(x$inf), names(x$sup))
	print(m)
	
	if(any(is.na(x$ev))) message("\nNotes:\nNaN is produced from one of evaluation points.\nmin.E(X|y) or max.E(X|y) may be incorrect.\n")
	
}

# TODO:
# model <- function(formula, data, dist=c("pois", "binom", "geo", "exp", "multinom"), ztrunc=FALSE){
# }




## Diaconis and Yxxxxx parameterization
## XXX: Depreciated
# 
#' @rdname update
#' @export 
update2.impinf <- function(object, y=NULL, ...){
	## this is a test module
	## 
	# based on diaconis and yyyyy parameterization 
	vtx <- as.matrix(object$vtx)
	colnames(vtx) <- c("xi1", "xi0")
	vtx1 <- vtx
	vtx1[,2] <- vtx1[,2] + rep(length(y), nrow(vtx1))
	vtx1[,1] <- vtx1[,1] + rep(sum(y), nrow(vtx1))
#	vtx1[,1] <- vtx1[,1]/vtx1[,2]
	em <- as.vector(vtx1[,1]/vtx1[,2])
	robj <- list(y=y, vtx=vtx, vtx1=vtx1, em=em)
	class(robj) <- c("impinf", "update2")
	return(robj)
}


