#' @rdname plot
#' @title Plotting Imprecise Objects
#' @description Generic function for plotting of \code{imprecise} objects
#' @param x an object for which a plot is needed.
#' @param after.obs logical indicating imprecise prior or posterior. 
#' @param ... additional arguments affecting the plot produced. 
#' @export 
plot.impinf <- function(x, after.obs=FALSE, ...){

	## TODO (2016.04.30):
	##  - support of lattice presentation
	##  - support of rgl presentation
	##  - support of both imprecise prior and posterior 

	
	if ( !after.obs ) vtx <- x$vtx
	else vtx <- x$vtx1

	# 2 dimension
	if ( ncol(vtx) == 2 ){
		graphics::plot(x=vtx, ...)
		graphics::polygon(x=vtx, col="azure2", border="darkblue", lty="dashed", lwd=3, ...)
		graphics::points(x=vtx, col="red", pch=19, cex=0.5)
		graphics::text(x=vtx, labels=rownames(vtx), pos=4)
	}
	
	# 3 dimension requires rgl package 
	if ( ncol(vtx) == 3 ) {
		vtx0 <- x$vtx0
		imat <- x$imat
		suppressWarnings(rgl::triangles3d(x=vtx0[imat,1],y=vtx0[imat,2],z=vtx0[imat,3], color="green"))
		rgl::decorate3d()
		rgl::texts3d(x=vtx[,1],y=vtx[,2],z=vtx[,3], text=rownames(vtx), col="red")
	}
	
}


#' @rdname plot
#' @title Probabilty Box for Imprecise Objects
#' @description Generating a probability box 
#' @param xmin lower limit of quantile
#' @param xmax upper limit of quantile
#' @param x.by argument used in seq on cdf
#' @param ncdfs number of CDFs
#' @param xmar adjustment for label
#' @export
pbox <- function(x, xmin=-4, xmax=4, ncdfs=5, x.by=0.05, xmar, ...){
	
	## TODO:
	## use the object produced from update() or summary()?
	## 

	vtx1 <- x$vtx1
	y <- x$y
	## TODO:
	## support mean parameterization
	inf.idx <- which.min(x$et)
	sup.idx <- which.max(x$et)
	inf.p1 <- vtx1[inf.idx,]
	sup.p1 <- vtx1[sup.idx,]
	pcpm <- Vectorize(pcpm, "q")

	t <- seq(from=xmin, to=xmax, by=x.by)
	col.cdf <- grDevices::rainbow(ncdfs)

	graphics::plot(0,0,type="n", ylim=c(0,1), xlim=c(xmin,xmax), xlab="", ylab="", ...)
	graphics::title(xlab=expression(theta), ylab="CDF", mgp=c(2,2,3), cex.lab=1.25, cex.main=1.25)
	graphics::abline(h=c(0,1), col="grey", lty="dashed")
	
	# get coords having extreme cdfs
	xtm <- rbind(inf.p1, sup.p1)
	colnames(xtm) <- c("xi2", "xi1", "xi0")
	xtm <- as.data.frame(xtm)

	# get coords of internal cdfs encolosed by two extreme cdfs
	pars.cdf <- cbind(xi2=seq(min(xtm$xi2), max(xtm$xi2), length.out=ncdfs), xi1=seq(min(xtm$xi1), max(xtm$xi1), length.out=ncdfs), xi0=seq(min(xtm$xi0), max(xtm$xi0), length.out=ncdfs))
	
	op <- list()
	for (i in 1:nrow(pars.cdf)) {
		ft <- pcpm(q=t, pars=c(pars.cdf[i,]))
		if (i %in% c(1,nrow(pars.cdf))) graphics::lines(x=t, y=ft, col=col.cdf[i], lwd=2)
		else graphics::lines(x=t, y=ft, col=col.cdf[i], lwd=1, lty=2)
		op[[i]] <- cbind(t=t, ft=ft)
	}

	# computing area of band enclosed by two cdfs
	# 
	ft.left <- unlist(op[[1]][,2])
	ft.right <- unlist(op[[nrow(pars.cdf)]][,2])
	ab <- sum(ft.left-ft.right)*x.by 

	# degree of imprecision 
	delta <- max(x$et) - min(x$et)

	# TODO:
	# xmar <- c(left, right)
	# XXX:
	# graphics::text(x=xmax-2, y=0.27, labels=bquote(n==.(length(y))), adj=1, cex=1.25)
	graphics::text(x=xmax-xmar[2], y=0.2, labels=bquote(paste(Delta,"[",E,"] = ", .(round(delta,4)), sep="")), adj=1, cex=1.25)
	graphics::text(x=xmax-xmar[2], y=0.15, labels=bquote(paste("area = ", .(round(ab,4)), sep="")), adj=1, cex=1.25)

	graphics::text(x=xmin+xmar[1], y=0.75, labels=bquote(paste("sup at (", .(paste(sup.p1, collapse=",")), ")", sep="")), adj=1, cex=1.25)
	graphics::text(x=xmin+xmar[1], y=0.7, labels=bquote(paste("inf at (", .(paste(inf.p1, collapse=",")), ")", sep="")), adj=1, cex=1.25)
	
	robj <- list(cdf=op, cdf.inf=ft.left, cdf.sup=ft.right, area=ab, delta=delta)
	return(robj)
}

# ui <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(-1,0,0),c(0,-1,0),c(0,0,-1))
# ci <- c(0,0,0,-1,-1,-1)
# op <- iprior(ui=ui, ci=ci)
# op

