# Numerical Tests
# ----------------
# 
# Modified on 2015-NOV-17
# Chel Hee Lee 

rm(list=ls())
# source("./../R/imPois.R")
# source("./../R/visualize.R")

## TEST EXAMPLE 1. 
p2 <- 4; p1 <- 1; p0 <- 0.5; x <- 1
# embl.pois.normal(x=numeric(0), v=p2, m=p1,n=p0, iprior="normal", wrt="canonical")$et  # [1] 0.05500347

op <- cgf(xi2=p2, xi1=p1, xi0=p0, ztrunc=TRUE)
op$value

# kcpm(t=x, xi2=p2, xi1=p1, xi0=p0, ztrunc=FALSE)/integrate(kcpm, lower=-5e2, upper=5e2, xi2=p2, xi1=p1, xi0=p0, ztrunc=FALSE)$value # [1] 0.02423068
kcpm(t=x, xi2=p2, xi1=p1, xi0=p0, ztrunc=TRUE)/integrate(kcpm, lower=-5e2, upper=5e2, xi2=p2, xi1=p1, xi0=p0, ztrunc=TRUE)$value # [1] 0.01995356

dcpm(x=x, pars=c(p2,p1,p0), ztrunc=TRUE) # [1] 0.02423068 (ztrunc = [1] 0.01995356)
pcpm(q=x, pars=c(p2,p1,p0), ztrunc=TRUE) # = 0.9973794 (from mathematica) # [1] 0.9973939 (from R), [1] 0.9978708 (for ztrunc)
evfn(y=x, pars=c(p2=0,p1,p0), ztrunc=TRUE) # [1] 0.01731923
# digamma(p1+sum(x)) - log(p0+length(x)) # [1] 0.01731923

p2 <- 0.25; p1 <- 1; p0 <- 1; x <- 2; y <- c(0,1);
dcpm(x=x, pars=c(p2,p1,p0), ztrunc=FALSE) # [1] 0.002220832
pcpm(q=x, pars=c(p2,p1,p0), ztrunc=FALSE) # [1] 0.9997336
evfn(y=numeric(0), pars=c(p2,p1,p0), ztrunc=FALSE)$value # [1] -0.2108353
evfn(y=numeric(0), pars=c(p2,p1+sum(y),p0+length(y)), ztrunc=FALSE)$value # [1] -0.4783651
evfn(y=y, pars=c(p2=p2, p1=p1, p0=p0), ztrunc=FALSE)$value # [1] -0.4783651

