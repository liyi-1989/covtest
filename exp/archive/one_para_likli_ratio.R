library(Rcpp)
source("utils.R")
Rcpp::sourceCpp('src/utils.cpp')
library(fields)
library(mvtnorm)
library(matrixcalc)

p=100; a=0.5; b=1; sigma=1; i0=5; j0=95
# true NULL model
A0=mdiag(p,c(1,a))
Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p))
# true alternative model
A1=A0
A1[i0,j0]=A1[j0,i0]=b
#A1[(j0-1):(j0+1),j0]=0
#A1[(j0-1):(j0+1),i0]=c(a,1,a)
Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p))

printM(A0)
printM(Sigma0)
printM(A1)
printM(Sigma1)

n=1000
data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma1)

# just to take a look at likelihood function (as a function of parameter a)
aa=seq(-5,5,by=0.1)
n1=length(aa)
laa=laa0=laa1=rep(0,n1)
for(i in 1:n1){
  laa0[i]=loglikelihood(data,aa[i],i0=i0,j0=j0,type=0)
  laa1[i]=loglikelihood(data,aa[i],i0=i0,j0=j0,type=1)
}
plot(aa,laa0/n,type="l",xlab="parameter a",ylab="log likelihood",main="")
points(aa,laa1/n,type="l",col="blue")
abline(v=a,col="red")

plot(aa,laa0/pmax(laa0,laa1),type="l",xlab="parameter a",ylab="log likelihood",main="")
abline(v=a,col="red")

# likelihood ratio test TS sampling distribution
res_lrt_0=res_lrt_1=rep(0,100)
for(k in 1:100){
  print(k)
  data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma0)
  res_lrt_0[k]=lrt_a(data,i0,j0)
  data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma1)
  res_lrt_1[k]=lrt_a(data,i0,j0)
  
}

plot(density(res_lrt_0))
points(density(res_lrt_1),col="red")

#######################
###### Functions ######
#######################
loglikelihood=function(X,a,sigma=1,i0,j0,type=0){
  n=dim(X)[1]
  p=dim(X)[2]
  S=cov(X)
  A0=mdiag.r(p,c(1,a))
  Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p))
  A1=A0
  A1[(j0-1):(j0+1),j0]=0
  A1[(j0-1):(j0+1),i0]=c(a,1,a)
  Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p))
  temp0=-2*p*log(2*pi)-log(det(Sigma0))-sum(diag(solve(Sigma0)%*%S))
  temp1=-2*p*log(2*pi)-log(det(Sigma1))-sum(diag(solve(Sigma1)%*%S))
  if(type==0){
    return((n/2)*temp0)
  }else{
    return((n/2)*temp1)
  }
}

lrt_a=function(data,i0,j0){
  fit0=optimize(loglikelihood,c(-5, 5), tol = 0.001,maximum = T, X=data,i0=i0,j0=j0,type=0)
  fit1=optimize(loglikelihood,c(-5, 5), tol = 0.001,maximum = T, X=data,i0=i0,j0=j0,type=1)
  
  l0=fit0$objective
  l1=fit1$objective
  
  return(2*(max(l1,l0)-l0))
}

sigmabar=function(X){
  n=dim(X)[1]
  p=dim(X)[2]
  S=cov(X)
  S1=S
  ahat=mean(diag(S[-p,-1]))/2
  for(i in 2:(p-1)){
    for(j in i:(p-1)){
      if(abs(i-j)>=1){
        S1[i,j]=S1[j,i]=(2*S[i,j]+ahat*(S[i,j+1]+S[i,j-1]+S[i+1,j]+S[i-1,j]))/(2+2*ahat^2)
      }
    }
  }
  
  return(S1)
}
