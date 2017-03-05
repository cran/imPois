## TODO: 2016.04.28
##
## Fig02.R shows a level set running from -4 to -4 by 1 on the expected 
## canonical parameter, E(theta).  What if we show the identical level based
## on the values running from -4 to -4 by 1 on the mean parameter? 

rm(list=ls())
source("./../R/imTools.R")
source("./../R/imPoisC.R")
source("./../R/imPoisM.R")
source("./../R/imGraphics.R")
# source("./../R/imZTPois.R")
# library(imPois)


# Grid on the hyperparameter space 
xi2 <- seq(from=0.05, to=3, by=0.1)
xi0 <- seq(from=0, to=5, by=0.25)
em <- seq(0,4,0.25)
hs <- expand.grid(xi2=xi2, xi0=xi0, em=em)

fn <- function(x){
	f0 <- function(xm, ...) em.pdf(y=numeric(0), pars=c(x[1], xm, x[2]))$value - x[3]
	val <- tryCatch(uniroot(f0, lower=0, upper=10, extendInt="yes", tol=1e-5)$root, error=function(e) return(NaN))
	return(val)
}
hs$xi1 <- apply(hs, 1, fn)

## 
## Codes from this point are in experiment
## 

## TODO: use 'rgl' package for dynamic presentation
## 
library(rgl)
m <- matrix(hs$xi1, ncol=length(xi0), nrow=length(xi2))
rgl.open()
rgl.clear()
rgl.bg(color="white")
rgl.points(x=hs$xi0, y=hs$xi2, z=hs$xi1, color="grey")
decorate3d()
aspect3d(1,1.5,2)

cols <- rainbow(length(em), s=0.7)
for(i in em){
	h1 <- subset(hs, subset=(em==i))
	m1 <- matrix(h1$xi1, ncol=length(xi0), nrow=length(xi2))
	surface3d(x=xi0, y=xi2, z=t(m1), color=cols[i+4])
}







# colors are used for reference levels E(t) = -4, -3, -2, -1, 0, 1, 2, 3, 4
# mycolors <- c("blue", "chocolate", "green", "red", "maroon", "orange", "purple", "sienna", "yellow")
mycolors <- grDevices::rainbow(length(em))

plot1 <- lattice::wireframe(xi1~xi0*xi2, data=hs, group=em, 
	xlab=list(label=expression(xi[0]), cex=1.25), 
	ylab=list(label=expression(xi[2]), cex=1.25),
	zlab=list(label=expression(xi[1]), cex=1.25), 
	zlim=c(0,10), xlim=c(0,4), ylim=c(0,4), 
	# screen=list(x=-35,y=-50,z=25), 
	screen=list(x=-55,y=-45,z=45), 
	key=list(
		text=list(as.character(em), col=mycolors), 
		points=list(pch=16, col=mycolors), 
		space="bottom", 
		columns=5
	), 
	col.groups=mycolors, 
	par.settings = list(
		axis.line = list(col = "transparent") 
	), 
	scales=list(arrows=FALSE, col="black", tck=0.6), margin=TRUE)


# Snapshot on each sample 
hs1 <- subset(hs, subset=(xi0 %in% c(0,1,2,5)))
# set.seed(18248439) # 1 1 2 0 1 2 2 1 0 1
set.seed(16979238) # 1,1,1,0,0,2,1,0,2,
y <- rpois(n=1e1, lambda=1)

# lc0 <- rbind(c(1,2,0), c(1,-2,0), c(0.5,-2,0), c(0.5,2,0)) # c(xi2, xi1, xi0) # E(theta) = 1, -1, -2, 2 
lc0 <- rbind(c(1,1,0), c(1,-1,0), c(0.1,-1,0), c(0.1,1,0)) # c(xi2, xi1, xi0)
op <- iprior(pmat=lc0)
op0 <- update(op, y=NULL, wrt="mean")
op1 <- update(op, y=y[1:1], wrt="mean")
op2 <- update(op, y=y[1:2], wrt="mean")
op5 <- update(op, y=y[1:5], wrt="mean")

rbind(op0$et, op1$et, op2$et, op5$et)


xleft <- c(min(op0$vtx1[,2]), min(op1$vtx1[,2]), min(op2$vtx1[,2]), min(op5$vtx1[,2]))
xright <- c(max(op0$vtx1[,2]), max(op1$vtx1[,2]), max(op2$vtx1[,2]), max(op5$vtx1[,2]))
ybottom <- c(min(op0$vtx1[,1]), min(op1$vtx1[,1]), min(op2$vtx1[,1]), min(op5$vtx1[,1]))
ytop <- c(max(op0$vtx1[,1]), max(op1$vtx1[,1]), max(op2$vtx1[,1]), max(op5$vtx1[,1]))

plot2 <- lattice::xyplot(xi2~xi1|as.factor(xi0), data=hs1, 
		group=em, 
		xlab=list(label=expression(xi[1]), cex=1.25), 
		ylab=list(label=expression(xi[2]), cex=1.25), 
		type="l", lwd=2, 
		layout=c(2,2), 
		auto.key=FALSE, 
		panel=function(x,y, ...){
			panel.rect(xleft=xleft[panel.number()], xright=xright[panel.number()], ybottom=ybottom[panel.number()], ytop=ytop[panel.number()], col="lightgray")
			panel.xyplot(x,y,...)
		}, 
	par.settings=list(
		superpose.line=list(col=mycolors),
		layout.heights=list(top.padding=3, bottom.padding=5), 
		layout.widths=list(left.padding=3, right.padding=3)
	), 
	between=list(x=2, y=2), 
	strip=strip.custom(var.name="n = ", bg="white", factor.levels=c("n = 0", "n = 1", "n = 2", "n = 5")) 
)

gridExtra::grid.arrange(plot1, plot2, ncol=2)

