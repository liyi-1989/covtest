# long short covariance model 
# bar: average surroundings
library(Rcpp)
source('utils_lscm.R')
Rcpp::sourceCpp('src/utils.cpp')
library(fields)
library(mvtnorm)
library(matrixcalc)



####### -1. True Covariance Model ######
####### -1.1 Basic Model ######
p=100; a=0.9; b=0.9; sigma=3; i0=5; j0=95
M1=Mp(p)
# true NULL model
A0=mdiag(p,c(1,a))
Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p))
# true alternative model
A1=A0
A1[i0,j0]=A1[j0,i0]=b #A1[(j0-1):(j0+1),j0]=0; A1[(j0-1):(j0+1),i0]=c(a,1,a)
Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p))

C=t(base::chol(Sigma1))
####### -1.2 Ding's (small) Model ######
# small p size, just to check if Ding's derivation is the same as quadratic programming
p=30; a=0.9; b=0.8; sigma=3; i0=12; j0=22
M1=Mp(p)
A0=mdiag(p,c(1,a)); Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p))
A1=A0; A1[j0,i0]=sqrt(b); A1[j0,j0]=sqrt(1-b); A1[j0+1,i0]=a*sqrt(b); A1[j0+1,j0]=a*sqrt(1-b)
Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p))
C1=t(base::chol(Sigma1))
# best lambda based on Ding's derivation
lambda_star_ding=(sigma^2+1)/((a^2-1)^2+(a^2+sigma^2)*(a^2+1))
c(lambda_star_ding,(1-lambda_star_ding)/a)
# best lambda base on true covcov quad-programming
d11=ijthcovcov(C=C1,M=M1,df=n,i1=i0,j1=j0,i2=i0,j2=j0)
d12=d21=ijthcovcov(C=C1,M=M1,df=n,i1=i0,j1=j0,i2=i0,j2=j0+1) 
d22=ijthcovcov(C=C1,M=M1,df=n,i1=i0,j1=j0+1,i2=i0,j2=j0+1) 

Dmat=matrix(c(d11,d21,d12,d22),2,2)
dvec=rep(0,2)
Amat=matrix(c(1,a),2,1)
bvec=1
fit.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
mu=fit.qp$solution
mu

# --- Conclusion: Ding and quad-prog have the same result (C is on true sigma1)
# best lambda base on sample covcov quad-programming
X=mvrnorm(n, mu = rep(0,p), Sigma = Sigma1)
C2=t(base::chol(cov(X)))
d11=ijthcovcov(C=C2,M=M1,df=n,i1=i0,j1=j0,i2=i0,j2=j0)
d12=d21=ijthcovcov(C=C2,M=M1,df=n,i1=i0,j1=j0,i2=i0,j2=j0+1) 
d22=ijthcovcov(C=C2,M=M1,df=n,i1=i0,j1=j0+1,i2=i0,j2=j0+1) 
Dmat=matrix(c(d11,d21,d12,d22),2,2)
fit2.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
mu2=fit2.qp$solution
mu2
# --- Conclusion: sampled version covcov will not differ too much for best lambda
printM(A0)
printM(Sigma0)
printM(A1)
printM(Sigma1)


####### -1.3 Ding's (large) Model ######
# small p size, just to check if Ding's derivation is the same as quadratic programming
p=100; a=0.9; b=0.8; sigma=1; i0=12; j0=92
M1=Mp(p)
A0=mdiag(p,c(1,a)); Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p))
A1=A0; 
A1[j0,i0]=sqrt(b); A1[j0,j0]=sqrt(1-b); 
A1[j0+1,i0]=a*sqrt(b); A1[j0+1,j0]=a*sqrt(1-b)
A1[j0-1,i0]=a*sqrt(b); A1[j0-1,j0]=a*sqrt(1-b)
Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p))
C1=t(base::chol(Sigma1))
# best lambda based on Ding's derivation
# lambda_star_ding=(sigma^2+1)/((a^2-1)^2+(a^2+sigma^2)*(a^2+1))
# c(lambda_star_ding,(1-lambda_star_ding)/a)
# best lambda base on true covcov quad-programming
Dmat=matrix(NA,5,5)
dfi0j0=data.frame(i0=c(i0,i0-1,i0,i0+1,i0),j0=c(j0,j0,j0+1,j0,j0-1))
for(i in 1:5){
  for(j in 1:5){
    i1=dfi0j0[i,"i0"]; j1=dfi0j0[i,"j0"]
    i2=dfi0j0[j,"i0"]; j2=dfi0j0[j,"j0"]
    Dmat[i,j]=ijthcovcov(C=C1,M=M1,df=n,i1=i1,j1=j1,i2=i2,j2=j2)
  }
}
dvec=rep(0,5)
Amat=matrix(c(1,a,a,a,a),5,1)
bvec=1
fit.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
mu=fit.qp$solution
mu


printM(A0)
printM(Sigma0)
printM(A1)
printM(Sigma1)
####### -1.4 2Dim Ding's (large) Model ######
# small p size, just to check if Ding's derivation is the same as quadratic programming
p=10; a=0.9; b=0.8; sigma=2; i0=23; j0=73
#M1=Mp(p^2)
A0=mdiag(p^2,c(1,a))
for(i in 1:(p^2)){
  if(i+p<=p^2) A0[i,i+p]=a
  if(i-p>=0) A0[i,i-p]=a
}
Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p^2))

A1=A0
A1[j0,i0]=sqrt(b); A1[j0,j0]=sqrt(1-b); 
A1[j0+1,i0]=a*sqrt(b); A1[j0+1,j0]=a*sqrt(1-b)
A1[j0-1,i0]=a*sqrt(b); A1[j0-1,j0]=a*sqrt(1-b)
Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p^2))
C1=t(base::chol(Sigma1))
Dmat=matrix(NA,5,5)
dfi0j0=data.frame(i0=c(i0,i0-1,i0,i0+1,i0),j0=c(j0,j0,j0+1,j0,j0-1))
for(i in 1:5){
  for(j in 1:5){
    i1=dfi0j0[i,"i0"]; j1=dfi0j0[i,"j0"]
    i2=dfi0j0[j,"i0"]; j2=dfi0j0[j,"j0"]
    Dmat[i,j]=ijthcovcov(C=C1,M=M1,df=n,i1=i1,j1=j1,i2=i2,j2=j2)
  }
}
dvec=rep(0,5)
Amat=matrix(c(1,a,a,a,a),5,1)
bvec=1
fit.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
mu=fit.qp$solution
mu

####### 0. Estimate Covariance of surrounding cov by simulation ######

n=1000

# # covcov by simulation
# N=1000
# i0j0=array(NA,c(3,3,N))
# for(k in 1:N){
#   # H0
#   data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma1)
#   S=cov(data)
#   i0j0[,,k]=S[(i0-1):(i0+1),(j0-1):(j0+1)]
# }
# XX=cbind(i0j0[2,2,],i0j0[1,2,],i0j0[2,3,],i0j0[3,2,],i0j0[2,1,])

# covcov by formula (gaussian)
Dmat=matrix(NA,5,5)
dfi0j0=data.frame(i0=c(i0,i0-1,i0,i0+1,i0),j0=c(j0,j0,j0+1,j0,j0-1))
for(i in 1:5){
  for(j in 1:5){
    i1=dfi0j0[i,"i0"]; j1=dfi0j0[i,"j0"]
    i2=dfi0j0[j,"i0"]; j2=dfi0j0[j,"j0"]
    Dmat[i,j]=ijthcovcov(C=C,M=M1,df=n,i1=i1,j1=j1,i2=i2,j2=j2)
  }
}

#Dmat=cov(XX)
dvec=rep(0,5)
Amat=matrix(c(2,a,a,a,a),5,1)
bvec=2
fit.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
mu=fit.qp$solution
mu
# check optimal variance (if variance reduced significantly)
Variance=t(mu)%*%Dmat%*%mu
mu0=c(1,0,0,0,0)
Variance0=t(mu0)%*%Dmat%*%mu0

Variance
Variance0
# check constraints
t(mu)%*%Amat
t(mu0)%*%Amat

# mean((XX%*%mu-2*b)^2)
# mean((XX%*%mu0-2*b)^2)
# 
# sd(XX%*%mu)
# sd(XX%*%mu0)

# solve quad-programming by hand
NUM=(2*solve(Dmat)%*%Amat)
DEN=(t(Amat)%*%solve(Dmat)%*%Amat)
NUM/DEN[1,1]

####### 1. Estimation problem ######
n=1000
N=100
sbar_0=shat_0=ahat_0=ahat2_0=ahat22_0=bbar_0=bhat_0=rep(NA,N)
sbar_1=shat_1=ahat_1=ahat2_1=ahat22_1=bbar_1=bhat_1=rep(NA,N)
for(k in 1:N){
  # H0
  data = mvrnorm(n, mu = rep(0,p^2), Sigma = Sigma0)
  S=cov(data)
  temp=sigmabar(data,S,a=NULL,mu=mu) # sigmabar
  Sbar=temp$S
  ahat_0[k]=temp$a
  ahat2_0[k]=temp$a2
  ahat22_0[k]=temp$a22
  sbar_0[k]=base::norm(Sbar-Sigma0)
  shat_0[k]=base::norm(S-Sigma0)
  bhat_0[k]=S[i0,j0]/2
  bbar_0[k]=Sbar[i0,j0]/2
  # H1
  data = mvrnorm(n, mu = rep(0,p^2), Sigma = Sigma1)
  S=cov(data)
  temp=sigmabar(data,S,a=NULL,mu=mu) # sigmabar
  Sbar=temp$S
  ahat_1[k]=temp$a
  ahat2_1[k]=temp$a2
  ahat22_1[k]=temp$a22
  sbar_1[k]=base::norm(Sbar-Sigma1,"F")
  shat_1[k]=base::norm(S-Sigma1,"F")
  bhat_1[k]=S[i0,j0]/2
  bbar_1[k]=Sbar[i0,j0]/2
}
mean(ahat_0)
sd(ahat_0)
mean(ahat2_0)
sd(ahat2_0)
mean(ahat22_0)
sd(ahat22_0)
mean(shat_0)
sd(shat_0)
mean(sbar_0)
sd(sbar_0)
mean(bhat_0)
sd(bhat_0)
mean(bbar_0)
sd(bbar_0)

mean(ahat_1)
sd(ahat_1)
mean(ahat2_1)
sd(ahat2_1)
mean(ahat22_1)
sd(ahat22_1)
mean(shat_1)
sd(shat_1)
mean(sbar_1)
sd(sbar_1)
mean(bhat_1)
sd(bhat_1)
mean(bbar_1)
sd(bbar_1)

par(mfrow=c(3,2))
plot(density(ahat_0),xlim=c(min(c(ahat_0,ahat2_0)),max(c(ahat_0,ahat2_0))),xlab="a",main="H0: est of a")
points(density(ahat2_0),col="red") # mean first, then sqrt
points(density(ahat22_0),col="blue") # sqrt fisrt (each elements), then mean [not good, bias+NaN]
legend("topright",c("est with 2a","est with a^2 (mean->sqrt)","est with a^2 (sqrt->mean)"),col=c("black","red","blue"),lty=c(1,1,1),cex=0.5)

plot(density(ahat_1),xlim=c(min(c(ahat_1,ahat2_1)),max(c(ahat_1,ahat2_1))),xlab="a",main="H1: est of a")
points(density(ahat2_1),col="red") # mean first, then sqrt
points(density(ahat22_1),col="blue") # sqrt fisrt (each elements), then mean [not good, bias+NaN]
legend("topright",c("est with 2a","est with a^2 (mean->sqrt)","est with a^2 (sqrt->mean)"),col=c("black","red","blue"),lty=c(1,1,1),cex=0.5)

plot(density(shat_0),xlim=c(min(c(shat_0,sbar_0)),max(c(shat_0,sbar_0))),xlab="MSE of cov estimation",main="H0: est of cov")
points(density(sbar_0),col="red") 
legend("topright",c("MSE of S","MSE of Sbar"),col=c("black","red"),lty=c(1,1),cex=0.75)

plot(density(shat_1),xlim=c(min(c(shat_1,sbar_1)),max(c(shat_1,sbar_1))),xlab="MSE of cov estimation",main="H1: est of cov")
points(density(sbar_1),col="red") 
legend("topright",c("MSE of S","MSE of Sbar"),col=c("black","red"),lty=c(1,1),cex=0.75)

plot(density(bhat_0),xlim=c(min(c(bhat_0,bbar_0)),max(c(bhat_0,bbar_0))),xlab="b",main="H0: est of b")
points(density(bbar_0),col="red") 
legend("topright",c("est b with S","est b with Sbar"),col=c("black","red"),lty=c(1,1),cex=0.5)

plot(density(bhat_1),xlim=c(min(c(bhat_1,bbar_1)),max(c(bhat_1,bbar_1))),xlab="b",main="H1: est of b")
points(density(bbar_1),col="red") 
legend("topright",c("est b with S","est b with Sbar"),col=c("black","red"),lty=c(1,1),cex=0.5)


####### 2. Test problem ######
######## 2.1. H0+H1 #########
# 10: H0 for cov
# 20: H0 for average cov
# true covariance under H0: Sigma0
n=1000 # sample size
N=100 # number of replicates
lmax10=lmax20=lmax30=rep(NA,N) # 1: coherence of S, 2: coherence of Sbar, 3: likelihood ratio, *0: null model 
lmax11=lmax21=lmax31=rep(NA,N) # 1: coherence of S, 2: coherence of Sbar, 3: likelihood ratio, *1: alternative model 

t0=proc.time()
for(k in 1:N){
  print(k)
  # H0
  data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma0)
  Shat=cov(data)
  lmax10[k]=max(abs(Shat)*(Sigma0==0))
  Sbar0=sigmabar(data,Shat,mu=mu)$S
  lmax20[k]=maxd.r(abs(Sbar0[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(Sigma0==0))
  lmax30[k]=lrt_a(data,Shat,i0,j0)
  # H1
  data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma1)
  Shat=cov(data)
  lmax11[k]=max(abs(Shat)*(Sigma0==0))
  Sbar=sigmabar(data,Shat,mu=mu)$S
  lmax21[k]=maxd.r(abs(Sbar[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(Sigma0==0))
  lmax31[k]=lrt_a(data,Shat,i0,j0)
}
t1=proc.time()-t0
t1

alpha=0.05
############ cut-off ############
# 10
TS10=n*lmax10^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
CUT.simu10=quantile(TS10,1-alpha)
CUT.asym10=-log(8*pi)-2*log(log(1/(1-alpha)))#y_alpha
# 20
TS20=n*lmax20^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
CUT.simu20=quantile(TS20,1-alpha)
CUT.asym20=-log(8*pi)-2*log(log(1/(1-alpha)))#y_alpha
# 30
TS30=lmax30
CUT.simu30=quantile(TS30,1-alpha)
CUT.asym30= qchisq(1-alpha,df=1)# chi-sq 
############ type I error ############
# size 10
mean(TS10>CUT.simu10) # by quantile definition, it is always 0.05
ES10=mean(TS10>CUT.asym10) # Empirical size
ES10
# size 20
mean(TS20>CUT.simu20) # by quantile definition, it is always 0.05
ES20=mean(TS20>CUT.asym20) # Empirical size
ES20
# size 30
mean(TS30>CUT.simu30) # by quantile definition, it is always 0.05
ES30=mean(TS30>CUT.asym30) # Empirical size
ES30

############ 1-type II error: power ############
# 11
TS11=n*lmax11^2-4*log(p)+log(log(p)) # coherence of S
EP11=mean(TS11>CUT.simu10) # Empirical power
EP11
# 21
TS21=n*lmax21^2-4*log(p)+log(log(p)) # coherence of Sbar
EP21=mean(TS21>CUT.simu20) # Empirical power
EP21
# 31
TS31=lmax31 # test statistics of likelihood ratio test
EP31=mean(TS31>CUT.simu30) # Empirical power
EP31

# size: is empirical size based on asy result (if wrong -> bad asy result)
par(mfrow=c(1,2))
# 11 plot rejection region
plot(density(TS10),xlim=c(min(density(TS10)$x,density(TS11)$x),max(density(TS10)$x,density(TS11)$x)),
     main=paste("Lmax: power",EP11,"size",ES10))
lines(density(TS11),col="red")
abline(v=CUT.simu10,lty=1,col="red")
abline(v=CUT.asym10,lty=2,col="blue")

# 21 plot rejection region
plot(density(TS20),xlim=c(min(density(TS20)$x,density(TS21)$x),max(density(TS20)$x,density(TS21)$x)),
     main=paste("Position: power",EP21,"size",ES20))
lines(density(TS21),col="red")
abline(v=CUT.simu20,lty=1,col="red")
abline(v=CUT.asym20,lty=2,col="blue")

mm=abs(Sbar0[2:(p-1),2:(p-1)])
which(mm == maxd.r(mm,4), arr.ind = TRUE)
mm=abs(Sbar[2:(p-1),2:(p-1)])
which(mm == maxd.r(mm,4), arr.ind = TRUE)
# # 31 plot rejection region
# plot(density(TS30),xlim=c(min(density(TS30)$x,density(TS31)$x),max(density(TS30)$x,density(TS31)$x)),
#      main=paste("Position: power",EP31,"size",ES30))
# lines(density(TS31),col="red")
# abline(v=CUT.simu30,lty=1,col="red")
# abline(v=CUT.asym30,lty=2,col="blue")



# check shrinkage
par(mfrow=c(1,2))
# H0
plot(density(TS10),xlim=c(min(density(TS10)$x,density(TS11)$x),max(density(TS10)$x,density(TS11)$x)),
     main=paste("Lmax: power",EP11,"size",ES10))
lines(density(TS20),col="red")

# H1
plot(density(TS11),xlim=c(min(density(TS20)$x,density(TS21)$x),max(density(TS20)$x,density(TS21)$x)),
     main=paste("Position: power",EP21,"size",ES20))
lines(density(TS21),col="red")

sd(TS10)
sd(TS20)

sd(TS11)
sd(TS21)


# loglikelihood=function(ab,X,S=NULL,sigma=1,i0,j0,type=0){
#   a=ab[1]
#   b=ab[2]
#   n=dim(X)[1]
#   p=dim(X)[2]
#   if(is.null(S)){
#     S=cov(X)
#   }
#   A0=mdiag.r(p,c(1,a))
#   Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p))
#   temp0=-2*p*log(2*pi)-log(det(Sigma0))-sum(diag(solve(Sigma0)%*%S))
#   
#   if(type==0){
#     return((n/2)*temp0)
#   }else{
#     A1=A0
#     A1[i0,j0]=A1[j0,i0]=b
#     Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p))
#     temp1=-2*p*log(2*pi)-log(det(Sigma1))-sum(diag(solve(Sigma1)%*%S))
#     return((n/2)*temp1)
#   }
# }
# 
# lrt_a=function(data,S,i0,j0){
#   # fit0=optimize(loglikelihood,c(-5, 5), tol = 0.001,maximum = T, X=data,S=S,i0=i0,j0=j0,b=b,type=0)
#   # fit1=optimize(loglikelihood,c(-5, 5), tol = 0.001,maximum = T, X=data,S=S,i0=i0,j0=j0,b=b,type=1)
#   # 
#   # l0=fit0$objective
#   # l1=fit1$objective
#   
#   fit0=optim(c(0,0), loglikelihood, control=list(fnscale=-1), X=data,S=S,i0=i0,j0=j0,type=0)
#   fit1=optim(c(0,0), loglikelihood, control=list(fnscale=-1), X=data,S=S,i0=i0,j0=j0,type=1)
#   
#   l0=fit0$value
#   l1=fit1$value
#   
#   return(2*(max(l1,l0)-l0))
# }
# 
