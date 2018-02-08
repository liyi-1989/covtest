library(Rcpp)
source("utils.R")
source('wishart_fun.R')
Rcpp::sourceCpp('src/utils.cpp')
library(fields)
library(mvtnorm)
library(matrixcalc)
library(quadprog)
library(Matrix)


p=20; a=0.9; r=0.8; sf=1; se=2; i0=3; j0=13
M1=Mp(p)
A0=mdiag(p,c(1,a)); Sigma0=(sf^2)*A0%*%t(A0)+se^2*diag(rep(1,p))
A1=A0; 
A1[j0,i0]=r; A1[j0,j0]=sqrt(1-r^2); 
A1[j0+1,i0]=a*r; A1[j0+1,j0]=a*sqrt(1-r^2)
A1[j0-1,i0]=a*r; A1[j0-1,j0]=a*sqrt(1-r^2)
Sigma1=(sf^2)*A1%*%t(A1)+se^2*diag(rep(1,p))
C1=t(base::chol(Sigma1))
# by formula
Dmat=matrix(NA,5,5)
dfi0j0=data.frame(i0=c(i0,i0-1,i0,i0+1,i0),j0=c(j0,j0,j0+1,j0,j0-1))
for(i in 1:5){
  for(j in 1:5){
    i1=dfi0j0[i,"i0"]; j1=dfi0j0[i,"j0"]
    i2=dfi0j0[j,"i0"]; j2=dfi0j0[j,"j0"]
    Dmat[i,j]=ijthcovcov(C=C1,M=M1,df=1,i1=i1,j1=j1,i2=i2,j2=j2)
  }
}
dvec=rep(0,5)
Amat=matrix(c(1,a,a,a,a),5,1)
bvec=1
fit.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
mu=fit.qp$solution
mu
# by hand
Dmat1=dmat(a,r,sf,se,i0,j0)$dmat
fit1.qp=solve.QP(Dmat1,dvec,Amat,bvec,meq=1)
mu1=fit1.qp$solution
mu1

(Dmat-Dmat1)>0.0001

l2=dmat(a,r,sf,se,i0,j0)$lambda2#
l2=((4*a^5-a)*sf^4+4*a^3*sf^2*se^2+a*se^4)/((16*a^6-10*a^4+a^2+1)*sf^4+(16*a^4-3*a^2+2)*sf^2*se^2+(4*a^2+1)*se^4)
c(1-4*a*l2,l2,l2,l2,l2)


n1=31
ai=li=li2=li3=rep(0,n1)
for(i in 1:n1){
  ai[i]=(i-1)/30
  li[i]=dmat(ai[i],r,sf,se-1,i0,j0)$lambda2
  li2[i]=dmat(ai[i],r,sf,se,i0,j0)$lambda2
  li3[i]=dmat(ai[i],r,sf,se+1,i0,j0)$lambda2
}
plot(ai,1-4*ai*li,type="l",xlim=c(0,1),ylim=c(0,1),xlab="a",ylab="lambda star",col="blue",lwd=1)
lines(ai,li,col="orange",lwd=1)
lines(ai,1-4*ai*li2,col="blue",lwd=2)
lines(ai,li2,col="orange",lwd=2)
lines(ai,1-4*ai*li3,col="blue",lwd=3)
lines(ai,li3,col="orange",lwd=3)
legend("topright",c("lambda for center","lambda for others"),col=c("blue","orange"),lty=c(1,1),cex=0.8)

dmat=function(a,r,sf,se,i0,j0){
  dmat=matrix(NA,5,5)
  #k0
  dmat[1,1]=(4*a^4+4*a^2+r^2+1)*sf^4+(4*a^2+2)*se^2*sf^2+se^4
  dmat[2,2]=dmat[3,3]=dmat[4,4]=dmat[5,5]=(4*a^4+4*a^2+a^2*r^2+1)*sf^4+(4*a^2+2)*se^2*sf^2+se^4
  # row 1
  dmat[1,2]=dmat[1,3]=dmat[1,4]=dmat[1,5]=(4*a^3+2*a+a*r^2)*sf^4+2*a*se^2*sf^2
  #k1,3
  dmat[2,3]=dmat[3,4]=dmat[4,5]=dmat[2,5]=(4*a^2+a^2*r^2)*sf^4
  #k2
  dmat[2,4]=dmat[3,5]=(2*a^4+a^2+a^2*r^2)*sf^4+a^2*se^2*sf^2
  for(i in 2:5){
    for(j in 1:(i-1)){
      dmat[i,j]=dmat[j,i]
    }
  }
  #--------------------------
  d1=(4*a^4+4*a^2+r^2+1)*sf^4+(4*a^2+2)*se^2*sf^2+se^4
  d2=(4*a^4+4*a^2+a^2*r^2+1)*sf^4+(4*a^2+2)*se^2*sf^2+se^4
  r1=(4*a^3+2*a+a*r^2)*sf^4+2*a*se^2*sf^2
  s1=(4*a^2+a^2*r^2)*sf^4
  s2=(2*a^4+a^2+a^2*r^2)*sf^4+a^2*se^2*sf^2
  return(list(dmat=dmat,lambda2=(a*d1-r1)/(d2+2*s1+s2-8*a*r1+4*a^2*d1)))
}

