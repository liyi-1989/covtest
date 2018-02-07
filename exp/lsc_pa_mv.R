# mean-variance
library(Rcpp)
source("utils.R")
Rcpp::sourceCpp('src/utils.cpp')
library(fields)
library(mvtnorm)
library(matrixcalc)

p=100; a=0.9; b=0.25; sigma=1; i0=5; j0=95
# true NULL model
A0=mdiag(p,c(1,a))
Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p))
# true alternative model
A1=A0
A1[i0,j0]=A1[j0,i0]=b #A1[(j0-1):(j0+1),j0]=0; A1[(j0-1):(j0+1),i0]=c(a,1,a)
Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p))

####### 1. Estimation problem ######
n=1000
N=100
i0j0=array(NA,c(3,3,N))
for(k in 1:N){
  # H0
  data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma0)
  S=cov(data)
  i0j0[,,k]=S[(i0-1):(i0+1),(j0-1):(j0+1)]
}

cor(i0j0[1,2,],i0j0[2,2,])

par(mfrow=c(3,3))
C0=matrix(NA,3,3)
for(i in 1:3){
  for(j in 1:3){
    C0[i,j]=cor(i0j0[i,j,],i0j0[2,2,])
    plot(i0j0[i,j,],i0j0[2,2,])
  }
}

XX=matrix(NA,N,5)
XX[,1]=i0j0[2,2]

XX=cbind(i0j0[2,2,],i0j0[1,2,],i0j0[2,3,],i0j0[3,2,],i0j0[2,1,])
C1=cov(XX)









