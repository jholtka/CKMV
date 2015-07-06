library(mvtnorm)
k <- 3
R <- diag(3)
R[2,1] <- 3/5
R[3,1] <- 1/3
R[3,2] <- 11/15
pmvnorm(mean=rep(0,k), R, lower=rep(-Inf,k), upper=c(1,4,2))