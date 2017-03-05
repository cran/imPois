## ./inst/fig04.R
## 
## Figure 4 shows two types of confliction between data and prior during
## a learning process.  A concept of agreement between two intentional units 
## who are learning the same observations.  
## 
## The concept of agreement should be discussed in the Discuss section. 
##
## This program will be moved to ./demo/conflict.R
##
## freezed on 2016.07.15
## 

rm(list=ls())
source("./../R/imTools.R")
source("./../R/imPoisC.R")
source("./../R/imPoisM.R")
source("./../R/imGraphics.R")
# source("./../R/imZTPois.R")
# library(imPois)



pdf(file="./figure4.pdf", width=10, height=6)
layout(matrix(c(1,2,2,3,4,4), nrow=2, byrow=TRUE), widths=c(1,2,1,2), height=c(1,1))
par(mar=c(4,4,2,0), oma=c(2,2,2,2))


# type 1
set.seed(16979238)
n <- 2e1
y <- rpois(n=n, lambda=1)

cvh0 <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
cvh1 <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(2,2,-3,-3)) 

plot(0, 0, type="n", xlab="", ylab="", xlim=c(0,5), ylim=c(0,5), cex.axis=1.25)
title(xlab=expression(xi[1]), ylab=expression(xi[0]), mgp=c(2,2,2), cex.lab=1.5, cex.main=1.75, face=2)
abline(a=0,b=1,col="grey",lwd=2)

polygon(cvh0$vtx[chull(cvh0$vtx),], lwd=2, lty="dashed", col=rgb(0,1,0,0.5))
mid0 <- colMeans(cvh0$vtx)
text(x=mid0[1], y=mid0[2], labels=expression(C[0]), cex=2)

polygon(cvh1$vtx[chull(cvh1$vtx),], lwd=2)
mid1 <- colMeans(cvh1$vtx)
text(x=mid1[1], y=mid1[2], labels=expression(C[1]), cex=2)


tmp0 <- list()
for(i in 1:length(y)) tmp0[[i]] <- summary(update(object=cvh0, y=y[1:i], wrt="mean"))
xm0 <- do.call(c, lapply(lapply(tmp0, "[[", "y"), mean))
xinf0 <- do.call(c, lapply(tmp0, "[[", "inf"))
xsup0 <- do.call(c, lapply(tmp0, "[[", "sup"))

tmp1 <- list()
for(i in 1:length(y)) tmp1[[i]] <- summary(update(object=cvh1, y=y[1:i], wrt="mean"))

xm1 <- do.call(c, lapply(lapply(tmp1, "[[", "y"), mean))
xinf1 <- do.call(c, lapply(tmp1, "[[", "inf"))
xsup1 <- do.call(c, lapply(tmp1, "[[", "sup"))

plot(x=0, y=0, type="n", ylim=c(0,2), xlim=c(0,n), xlab="", ylab="", cex.axis=1.25)
title(xlab="", ylab="E(Y)", mgp=c(2,2,2), cex.lab=1.5, cex=1.5)
abline(h=1, col="grey")
segments(x0=1:n, y0=xinf0, x1=1:n, y1=xsup0, lwd=12, col=rgb(0,1,0,0.5))
segments(x0=1:n, y0=xinf1, x1=1:n, y1=xsup1, lwd=2)
points(x=1:n, y=xm1, pch="-", cex=3, col="red")
text(x=1:n, y=rep(0,n), labels=y, cex=1.5)
legend("topright", legend=c(expression("upper and lower bounds from"~C[0]), expression("upper and lower bounds from"~C[1]), "sample mean"), lwd=2, col=c(rgb(0,1,0,0.5), "black", "red"), cex=1.25, bty="n")



# type 2.
#
# 20 random variates
n <- 2e1 
set.seed(16979238) # 1,1,1,0,0,2,1,0,2,
y <- rpois(n=n, lambda=1)

cvh0 <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
cvh2 <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(2,0,-3,-1)) 
cvh3 <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,2,-1,-3)) 

plot(0, 0, type="n", xlab="", ylab="", xlim=c(0,5), ylim=c(0,5), cex.axis=1.25)
title(xlab=expression(xi[1]), ylab=expression(xi[0]), mgp=c(2,2,2), cex.lab=1.5, cex.main=1.5, face=2)
abline(a=0,b=1,col="grey",lwd=2)

polygon(cvh0$vtx[chull(cvh0$vtx),], lwd=2, lty="dashed")
mid0 <- colMeans(cvh0$vtx)
text(x=mid0[1], y=mid0[2], labels=expression(C[0]), cex=2)

polygon(cvh2$vtx[chull(cvh2$vtx),], lwd=2, col=rgb(1,0,0,0.4))
mid2 <- colMeans(cvh2$vtx)
text(x=mid2[1], y=mid2[2], labels=expression(C[2]), cex=2)

polygon(cvh3$vtx[chull(cvh3$vtx),], lwd=2, col=rgb(0,0,1,0.4))
mid3 <- colMeans(cvh3$vtx)
text(x=mid3[1], y=mid3[2], labels=expression(C[3]), cex=2)

tmp2 <- list()
for(i in 1:length(y)) tmp2[[i]] <- summary(update(object=cvh2, y=y[1:i], wrt="mean"))
xm2 <- do.call(c, lapply(lapply(tmp2, "[[", "y"), mean))
xinf2 <- do.call(c, lapply(tmp2, "[[", "inf"))
xsup2 <- do.call(c, lapply(tmp2, "[[", "sup"))

tmp3 <- list()
for(i in 1:length(y)) tmp3[[i]] <- summary(update(object=cvh3, y=y[1:i], wrt="mean"))
xm3 <- do.call(c, lapply(lapply(tmp3, "[[", "y"), mean))
xinf3 <- do.call(c, lapply(tmp3, "[[", "inf"))
xsup3 <- do.call(c, lapply(tmp3, "[[", "sup"))

plot(x=0, y=0, type="n", ylim=c(0,2), xlim=c(0,n), xlab="", ylab="", cex.axis=1.25)
title(xlab="", ylab="E(Y)", mgp=c(2,2,2), cex.lab=1.5, cex=1.5)
abline(h=1, col="grey")
segments(x0=1:n, y0=xinf2, x1=1:n, y1=xsup2, lwd=6, col=rgb(0,0,1,0.4))
segments(x0=1:n, y0=xinf3, x1=1:n, y1=xsup3, lwd=6, col=rgb(1,0,0,0.4))
points(x=1:n, y=xm2, pch="-", cex=3, col="red")
text(x=1:n, y=rep(0,n), labels=y, cex=1.5)
legend("topright", legend=c(expression("upper and lower bounds"~C[2]), expression("upper and lower bounds from"~C[3]), "sample mean"), lwd=2, col=c(rgb(1,0,1,0.4), rgb(0,0,1,0.4), rgb(1,0,0,1)), cex=1.25, bty="n")

dev.off()

