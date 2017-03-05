## ./inst/fig05.R
## 
## Figure 5 illustrates a focusing behaviour of imprecise probabilities using
## a probability box.  Note that a degree of imprecision has the same amount 
## of quantity that the band emclosed by two extreme probabilities has. 
## 
##
## 
## This program will be moved to ./demo/focus.R
## 
## freezed on 2016.04.30
## 


rm(list=ls())
source("./../R/imTools.R")
source("./../R/imPoisC.R")
source("./../R/imPoisM.R")
source("./../R/imGraphics.R")
# source("./../R/imZTPois.R")
# library(imPois)


set.seed(16979238)
lambda <- 1
n <- 1e2
y <- rpois(n=n, lambda=lambda)

# lc0 <- rbind(c(0.2,1,0), c(0.2,-1,0), c(0.5,-1,0), c(0.5,1,0)) # c(xi2, xi1, xi0)
lc0 <- rbind(c(0.2,-1,0), c(0.2,1,0), c(1,-1,0), c(1,1,0)) # the same as figure 2,3,4
op <- iprior(pmat=lc0)

# hs <- expand.grid(xi2=seq(min(lc0[,1]),max(lc0[,1]),0.1), xi1=seq(min(lc0[,2]), max(lc0[,2]), 0.1), xi0=c(0,1,2,5,10))
# hs$et <- apply(hs, 1, function(x) et.pdf(y=numeric(0), pars=x)$value)
# library(lattice)
# wireframe(et~xi1*xi2, data=hs, group=xi0, scales=list(arrows=FALSE))

op0 <- update(op, y=NULL, wrt="canonical")
op1 <- update(op, y=y[1:1], wrt="canonical")
op2 <- update(op, y=y[1:2], wrt="canonical")
op5 <- update(op, y=y[1:5], wrt="canonical")
# op10 <- update(op, y=y[1:10], wrt="canonical")

pdf(file="./figure5.pdf", height=4, width=12)
par(mfrow=c(1,4), mar=c(4,4,4,1))
pb0 <- pbox(op0, xmin=-10, xmax=10, xmar=c(8,1), main="n = 0", cex.main=1.5)
pb1 <- pbox(op1, xmin=-8, xmax=5, xmar=c(6,0.5), main="n = 1", cex.main=1.5)
pb2 <- pbox(op2, xmin=-8, xmax=5, xmar=c(6,0.5), main="n = 2", cex.main=1.5)
pb5 <- pbox(op5, xmin=-8, xmax=5, xmar=c(6,0.5), main="n = 5", cex.main=1.5)
dev.off()


##
## Codes from this points are in experiment
## 

## normal prior
xtms <- expand.grid(x=c(1,2), y=c(1,2), z=c(1,2)) # overlapped
xtms <- rbind(c(1,1,2), c(1,1.2,1.8), c(1,1.4,1.6), c(1,1.6,1.4), c(1,1.8,1.2), c(1,2,1)) # ok when xi2 is fixed
# xtms <- rbind(c(1,2,1), c(1.2,1.8,1), c(1.4,1.6,1), c(1.6,1.4,1), c(1.8,1.2,1), c(2,1,1)) # ok when xi0 is fixed
# xtms <- rbind(c(1,1,2), c(1.2,1,1.8), c(1.4,1,1.6), c(1.6,1,1.4), c(1.8,1,1.2), c(2,1,1)) # ok when xi1 is fixed


op <- iprior(pmat=xtms)
op1 <- update(op, wrt="canonical")
t <- seq(-3,3,0.1)

plot(0,0,type="n", ylim=c(0,1), xlim=c(-3,3))
tmp <- list()
for(i in 1:nrow(pmat)){
	tmp[[i]] <- ft <- pcpm(q=t, pars=unlist(xtms[i,]))
	lines(t,ft, col=i+1)
}
max(op1$et) - min(op1$et)
sum(tmp[[1]] - tmp[[6]])*0.1

# computation of area in the band - width*height
# width = 0.1
# height = ft1 - ft6



rm(list=ls())
source("./../R/imPois.R")
source("./../R/imTools.R")
source("./../R/imPoism.R")
source("./../R/imGraphics.R")

set.seed(16979238)
lambda <- 1
n <- 1e2
y <- rpois(n=n, lambda=lambda)

## log-gamma prior (when xi2=0)
pmat <- rbind(c(1.0,2.0),c(1.2,1.8),c(1.4,1.6),c(1.6,1.4),c(1.8,1.2),c(2.0,1.0)) # is it normal? gamma?
op <- iprior(pmat=pmat)
op1 <- update(op, y=y[1:5], wrt="canonical")

pcpm <- Vectorize(pcpm, "q")
t <- seq(-3,1,0.1)

plot(0,0,type="n", ylim=c(0,1), xlim=c(-3,1))
tmp <- list()
for(i in 1:nrow(pmat)){
	ft <- pcpm(q=t, pars=c(0,op1$vtx1[i,]))
	tmp[[i]] <- ft
	lines(t,ft, col=i+1)
}
max(op1$et) - min(op1$et)
sum(tmp[[1]]-tmp[[6]])*0.1
# mean(tmp[[1]]-tmp[[6]])


# tmp1 <- do.call(cbind, tmp)
# dif <- apply(tmp1, 1, function(x) max(x)-min(x))
# mean(dif)

