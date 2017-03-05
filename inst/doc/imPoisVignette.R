## ------------------------------------------------------------------------
library(imPois)

## ----include=FALSE-------------------------------------------------------
sumy <- 32*1+2*16+3*6+4*1
n <- 32+16+6+1 

em.min <- em.rt.ztrunc(y=numeric(0), pars=c(0,0+sumy,1+n))$value # minimum from (0,1) #$
em.max <- em.rt.ztrunc(y=numeric(0), pars=c(0,1+sumy,0+n))$value # maximum from (1,0) #$

## ----include=FALSE-------------------------------------------------------
## MLE using the explicit formula given by Irwin 1959.

x <- c(1,2,3,4)
nx <- c(32,16,6,1)
xbar <- sum(nx*x)/sum(nx)

j <- 0
rv <- 0
epsilon <- 1
while(epsilon > 1e-8){
	j <- j+1
	rv1 <- rv + j^(j-1)/factorial(j)*(xbar*exp(-xbar))^j
	epsilon <- abs(rv1-rv)
	rv <- rv1
}
lambda <- xbar - rv

