## ./inst/fig03.R
## 
## Figure 3 shows a behaviour of soft-linearity on the expected values of 
## canonical parametr using contour and surface plots. 
## 
## imprecise normal prior specification
## rbind(c(0.2,-1,0), c(0.2,1,0), c(1,-1,0), c(1,1,0)) 
##
## seed: 16979238 
##
## This program will be moved to ./demo/softlinear.R
##
## freezed on 2016.04.30


rm(list=ls())
source("./../R/imTools.R")
source("./../R/imPoisC.R")
source("./../R/imPoisM.R")
source("./../R/imGraphics.R")
# source("./../R/imZTPois.R")
# library(imPois)


## contour plot 
panel.3d.contour <- function(x, y, z, rot.mat, distance, nlevels = 10, zlim.scaled, ...) { 
	add.line <- trellis.par.get("add.line") 
	panel.3dwire(x, y, z, rot.mat, distance, zlim.scaled = zlim.scaled, ...) 
	clines <- contourLines(x, y, matrix(z, nrow = length(x), byrow = TRUE), nlevels = nlevels) 
	for (ll in clines) { 
		m <- ltransform3dto3d(rbind(ll$x, ll$y, zlim.scaled[2]), rot.mat, distance) 
		panel.lines(m[1,], m[2,], col = add.line$col, lty = add.line$lty, lwd = add.line$lwd) 
	} 
}

set.seed(16979238) # 1,1,1,0,0,2,1,0,2,
y <- rpois(n=1e1, lambda=1)

# imprecise normal prior specification 
vtx <- rbind(c(0.2,-1,0), c(0.2,1,0), c(1,-1,0), c(1,1,0))

# TODO:
# can we develop this optim function?
# ui <- rbind(c(1,0,0),c(-1,0,0),c(0,1,0),c(0,-1,0),c(0,0,1),c(0,0,-1))
# ci <- c(0.2, -1,-1,-1,0,0)
# op <- iprior(ui=ui, ci=ci)
# minimum at (1,1,2) $value [1] -0.3261499
# constrOptim(theta=c(0.5,0.5,1), f=fn.em.normal, grad=gr.em.normal, ui=ui, ci=ci, method="BFGS", hessian=TRUE) 

hs0 <- expand.grid(xi2=seq(min(vtx[,1]),max(vtx[,1]),len=9), xi1=seq(min(vtx[,2]),max(vtx[,2]),len=9), xi0=0)
hs0$et <- apply(hs0, 1, function(x) et.pdf(y=numeric(0), pars=x)$value)

hs1 <- expand.grid(xi2=seq(min(vtx[,1]),max(vtx[,1]),len=9), xi1=seq(min(vtx[,2]),max(vtx[,2]),len=9)+sum(y[1:1]), xi0=length(y[1:1]))
hs1$et <- apply(hs1, 1, function(x) et.pdf(y=numeric(0), pars=x)$value)

hs2 <- expand.grid(xi2=seq(min(vtx[,1]),max(vtx[,1]),len=9), xi1=seq(min(vtx[,2]),max(vtx[,2]),len=9)+sum(y[1:2]), xi0=length(y[1:2]))
hs2$et <- apply(hs2, 1, function(x) et.pdf(y=numeric(0), pars=x)$value)

hs5 <- expand.grid(xi2=seq(min(vtx[,1]),max(vtx[,1]),len=9), xi1=seq(min(vtx[,2]),max(vtx[,2]),len=9)+sum(y[1:5]), xi0=length(y[1:5]))
hs5$et <- apply(hs5, 1, function(x) et.pdf(y=numeric(0), pars=x)$value)

minz <- -3.0
maxz <- 2.5
mainy <- -10


## surface and contour plots when n=0
plot0 <- wireframe(et ~ xi1*xi2, data=hs0, 
	xlab=list(expression(xi[1]), cex=1.25), ylab=list(expression(xi[2]), cex=1.5), zlab=list(expression(paste("E[",theta,"]",sep="")), cex=1.25), zlim=c(minz, maxz),
 	panel.3d.wireframe=panel.3d.contour, 
	main=expression("n = 0"), 
	par.settings=list(
		par.main.text=list(y=grid::unit(mainy,"mm")), 
		axis.line = list(col = "transparent")
	), 	
	scales=list(arrows=FALSE, col="black", tck=0.6)
)


## surface and contour plots when n=1
plot1 <- wireframe(et ~ xi1*xi2, data=hs1, 
	xlab=list(expression(xi[1]), cex=1.25), ylab=list(expression(xi[2]), cex=1.5), zlab=list(expression(paste("E[",theta,"]",sep="")), cex=1.25), zlim=c(minz, maxz),
 	panel.3d.wireframe=panel.3d.contour, 
	main=expression("n = 1"), 
	par.settings=list(
		par.main.text=list(y=grid::unit(mainy,"mm")), 
		axis.line = list(col = "transparent")
	), 	
	scales=list(arrows=FALSE, col="black", tck=0.6)
)


## surface and contour plots when n=2
plot2 <- wireframe(et ~ xi1*xi2, data=hs2, 
	xlab=list(expression(xi[1]), cex=1.25), ylab=list(expression(xi[2]), cex=1.5), zlab=list(expression(paste("E[",theta,"]",sep="")), cex=1.25), zlim=c(minz, maxz),
 	panel.3d.wireframe=panel.3d.contour, 
	main=expression("n = 2"), 
	par.settings=list(
		par.main.text=list(y=grid::unit(mainy,"mm")), 
		axis.line = list(col = "transparent")
	), 	
	scales=list(arrows=FALSE, col="black", tck=0.6)
)


## surface and contour plots when n=5
plot5 <- wireframe(et ~ xi1*xi2, data=hs5, 
	xlab=list(expression(xi[1]), cex=1.25), ylab=list(expression(xi[2]), cex=1.5), zlab=list(expression(paste("E[",theta,"]",sep="")), cex=1.25), zlim=c(minz, maxz),
 	panel.3d.wireframe=panel.3d.contour, 
	main=expression("n = 5"), 
	par.settings=list(
		par.main.text=list(y=grid::unit(mainy,"mm")), 
		axis.line = list(col = "transparent")
	), 	
	scales=list(arrows=FALSE, col="black", tck=0.6)
)


## presenting multiple plots 
pdf(file="figure3.pdf", width=15, height=6)
gridExtra::grid.arrange(plot0, plot1, plot2, plot5, ncol=4)
dev.off()


##
## Codes from this point are in experiment
## 


