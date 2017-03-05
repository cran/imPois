## developed on 2017-02-26
##
## ./inst/fig01.R
##
## Figure 1 shows a translation behaviour of imprecise log-gamma prior (i.e.
## xi2=0) on both canonical and mean parameter spaces.  This example is an 
## instance of no confliction between data and prior. 
## 
##

##########################################
##
## Contour lines of E(Y), mean parameter
##
##########################################


rm(list=ls())
source("./../R/imTools.R")
source("./../R/imPoisC.R")
source("./../R/imPoisM.R")
source("./../R/imGraphics.R")
# source("./../R/imZTPois.R")
# library(imPois)


set.seed(16979238) # 1,1,1,0,0,2,1,0,2,
y <- rpois(n=1e1, lambda=1)

xi1 <- seq(from=0, to=10, by=0.05) # alpha: shape
mlvls <- c(0.05, 0.1, 0.3, 0.5, 0.7,1,1.5,2.5,5,10,21)

hs <- expand.grid(xi1=xi1, em=mlvls)

fn <- function(x){
	f0 <- function(xm, ...) em.pdf(y=numeric(0), pars=c(0,x[1],xm))$value-x[2]
	val <- tryCatch(uniroot(f0, lower=1e-2, upper=10, extendInt="yes", tol=1e-5)$root, error=function(e) return(NaN))
	return(val)
} 
hs$xi0 <- apply(hs, 1, fn) # the estimated xi0 using given 



pdf(file="./contours_EY.pdf")


plot(x=hs$xi1, y=hs$xi0, cex=0.5, pch=19, type="n", xlim=c(0,10.5), ylim=c(0,10.5), xlab="", ylab="")
title(xlab=expression(xi[1]~"(sum)"), cex.lab=1.25, line=2)
title(ylab=expression(xi[0]~"(number of observations)"), cex.lab=1.25, line=2)
axis(side=2, at=1:11, labels=1:11)
abline(a=0,b=1, lty="dashed", lwd=1.5, col=rgb(0,1,0,0.5))
polygon(x=c(0,0,1,1),y=c(0,1,1,0), border="darkblue", lwd=2)
cvh <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 

# all 
for (i in 1:10) {
	op <- update(object=cvh, y=y[1:i], wrt="mean")
	polygon(op$vtx1, border="grey", lwd=1.5)
	vtx1.m <- colMeans(op$vtx1)
	text(x=0.2+vtx1.m[1], y=0.2+vtx1.m[2], labels=bquote(.(y[i])), cex=1.5)
}

# see expectations
for(v in mlvls){
	hsi <- subset(hs, subset=(em==v))
	lines(x=hsi$xi1, hsi$xi0, col="grey", lty="dashed", lwd=1.5)
	# points(x=hsi$xi1, y=hsi$xi0, cex=0.3, pch=19)
}
text(x=c(0.1,3,5,7,10),y=10, labels=c(0.05, 0.3, 0.5, 0.7, 1), pos=3)
text(x=10,y=c(7,4,2,1,0.5), labels=c(1.5,2.5,5,10,21), pos=4)

tmp <- list()
for (i in 1:length(y)) {
	op <- update(object=cvh, y=y[1:i], wrt="mean")
	sop <- summary(op)
	text(x=(sop$inf.p1)[1], y=(sop$inf.p1)[2], labels=round(sop$inf,3), col="darkblue", cex=1, pos=2)
	text(x=(sop$sup.p1)[1], y=(sop$sup.p1)[2], labels=round(sop$sup,3), col="red", cex=1, pos=4)
	tmp[[i]] <- sop
}
legend("bottomright", legend=c("min", "max"), ncol=2, text.col=c("darkblue", "red"), bty="n")
polygon(x=c(0,1),y=c(1,0), col="darkblue", lwd=2)
opc <- tmp # for numerical summary in the manuscript 

for (i in 1:length(y)) {
	op <- update(object=cvh, y=y[1:i], wrt="mean")
	sop <- summary(op)
	polygon(x=c(sop$inf.p1[1], sop$sup.p1[1]), y=c(sop$inf.p1[2], sop$sup.p1[2]), lwd=2)
}
mtext(side=1, line=3, text="Number of cases", cex=1.5)
mtext(side=2, line=4, text="Number of patients that you have seen", cex=1.5)


dev.off()

###############################################
## 
## contour lines of E(t), canonical parameter 
## 
###############################################

tlvls <- c(-3,-2,-1,-0.5,0,0.5,1,2,3)
hs.t <- expand.grid(xi1=xi1, em=tlvls)
fn1 <- function(x){
	f0 <- function(xt, ...) et.pdf(y=numeric(0), pars=c(0,x[1],xt))$value-x[2]
	val <- tryCatch(uniroot(f0, lower=1e-2, upper=10, extendInt="yes", tol=1e-5)$root, error=function(e) return(NaN))
	return(val)
}
hs.t$xi0 <- apply(hs.t, 1, fn1)


pdf(file="contours_ET.pdf")


plot(x=hs.t$xi1, y=hs.t$xi0, cex=0.5, pch=19, type="n", xlim=c(0,10.5), ylim=c(0,10.5), xlab="", ylab="")
title(xlab=expression(xi[1]~"(sum)"), cex.lab=1.25, line=2)
title(ylab=expression(xi[0]~"(number of observations)"), cex.lab=1.25, line=2)
axis(side=2, at=1:11, labels=1:11)
abline(a=0,b=1, lty="dashed", lwd=1.5, col=rgb(0,1,0,0.5))
polygon(x=c(0,0,1,1),y=c(0,1,1,0), border="darkblue", lwd=2)
cvh <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 

for (i in 1:10) {
	op <- update(object=cvh, y=y[1:i], wrt="canonical")
	polygon(op$vtx1, border="grey", lwd=1.5)
	vtx1.m <- colMeans(op$vtx1)
	text(x=0.2+vtx1.m[1], y=0.2+vtx1.m[2], labels=bquote(.(y[i])), cex=1.5)
}

# see expectations
for(v in tlvls){
	hsi <- subset(hs.t, subset=(em==v))
	lines(x=hsi$xi1, hsi$xi0, col="grey", lty="dashed", lwd=1.5)
	# points(x=hsi$xi1, y=hsi$xi0, cex=0.3, pch=19)
}
text(x=c(1,2,4.5,7,10.5),y=10, labels=c(-3,-2,-1,-0.5,0), pos=3)
text(x=10,y=c(6,3.5,1.5,0.5), labels=c(0.5,1,2,3), pos=4)

tmp <- list()
for (i in 1:length(y)) {
	op <- update(object=cvh, y=y[1:i], wrt="canonical")
	sop <- summary(op)
	text(x=(sop$inf.p1)[1], y=(sop$inf.p1)[2], labels=round(sop$inf,3), col="darkblue", cex=1, pos=2)
	text(x=(sop$sup.p1)[1], y=(sop$sup.p1)[2], labels=round(sop$sup,3), col="red", cex=1, pos=4)
	tmp[[i]] <- sop
}
legend("bottomright", legend=c("min", "max"), ncol=2, text.col=c("darkblue", "red"), bty="n")
polygon(x=c(0,1),y=c(1,0), col="darkblue", lwd=2)

for (i in 1:length(y)) {
	op <- update(object=cvh, y=y[1:i], wrt="canonical")
	sop <- summary(op)
	polygon(x=c(sop$inf.p1[1], sop$sup.p1[1]), y=c(sop$inf.p1[2], sop$sup.p1[2]), lwd=2)
}
mtext(side=1, line=3, text="Number of cases", cex=1.5)
mtext(side=2, line=4, text="Number of patients that you have seen", cex=1.5)

opm <- tmp ## for numerical summary in the manuscript 


dev.off()



## Codes from this point are used for numerical summaries 
## use Hmisc::latex(object, file="", digitis=3, ) . 

r1c <- do.call(c, lapply(opc, "[[", "sup"))
r2c <- do.call(c, lapply(opc, "[[", "inf"))
r3c <- do.call(c, lapply(opc, "[[", "delta"))
tbc <- rbind(r1c, r2c, r3c)

r1m <- do.call(c, lapply(opm, "[[", "sup"))
r2m <- do.call(c, lapply(opm, "[[", "inf"))
r3m <- do.call(c, lapply(opm, "[[", "delta"))
tbm <- rbind(r1m, r2m, r3m)

idx <- c(1,2,5,10)
tbc <- tbc[,idx]
tbm <- tbm[,idx]

tb <- cbind(tbc, tbm)

##
## codes from this points are in experiment
##


##################################################################
##
## Comparing E(t) and E(Y) over the same hyperparameters (e1,e0)
## 
##################################################################

rm(list=ls())
source("./../R/imTools.R")
source("./../R/imPoisC.R")
source("./../R/imPoisM.R")
source("./../R/imGraphics.R")
# source("./../R/imZTPois.R")
# library(imPois)


set.seed(16979238) # 1,1,1,0,0,2,1,0,2,
y <- rpois(n=1e1, lambda=1)

xi1 <- seq(from=0, to=10, by=0.05) # alpha: shape
mlvls <- c(0.05, 0.1, 0.3, 0.5, 0.7, 1, 1.5, 2.5, 5, 10, 21)

hs <- expand.grid(xi1=xi1, em=mlvls)

fn <- function(x){
	f0 <- function(xm, ...) em.pdf(y=numeric(0), pars=c(0,x[1],xm))$value-x[2]
	val <- tryCatch(uniroot(f0, lower=1e-2, upper=10, extendInt="yes", tol=1e-5)$root, error=function(e) return(NaN))
	return(val)
} 
hs$xi0 <- apply(hs, 1, fn) # the estimated xi0 using given mean levels and 'xi1'. 

hs$et <- apply(hs, 1, function(x){
	val <- tryCatch(et.pdf(y=numeric(0), pars=c(0,x[1],x[3]))$value, error=function(e) return(NaN))
	return(val)
})


require(rgl)
rgl.open()
rgl.bg(color="white")
plot3d(x=hs$xi1, y=hs$xi0, z=hs$em)

rgl.open()
rgl.bg(color="white")
plot3d(x=hs$xi1, y=hs$xi0, z=hs$et)
