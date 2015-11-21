# test of convh and visualization
# model() -> iprior() -> update() -> summary() -> plot() -> print()

# source("./../R/imPois.R")
# source("./../R/visualize.R")

# in 2 dimensions
op <- iprior(ui=rbind(c(1,0),c(0,1),c(0,-1),c(1,1),c(-2,-1)), ci=c(1,2,-8,5,-14)) # (3,8),(1,8), (1,4),(3,2)(6,2)
update.impinf(op, y=numeric(0)) 



# in 3 dimensions
# source("./../R/imPois.R")
# source("./../R/visualize.R")
op <- iprior(ui=rbind(diag(3), -diag(3)), ci=-rep(1,6))
# scatterplot3d(op$vtx)
op$vtx <- op$vtx*0.5 + 1
update.impinf(op, y=numeric(0))

