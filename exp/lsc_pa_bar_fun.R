library(Rcpp)
source("utils.R")
source('wishart_fun.R')
Rcpp::sourceCpp('src/utils.cpp')
library(fields)
library(mvtnorm)
library(matrixcalc)
library(quadprog)
library(Matrix)

###### simulation function ######

simu=function(B,bi,M1){
  b=0.2
  n=B[bi]#1000 # sample size
  ######### 0. setup #########
  p=100; a=0.9; sigma=1; i0=5; j0=95
  # true NULL model
  A0=mdiag(p,c(1,a))
  Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p))
  A1=A0
  # --- true alternative model (basic model)
  #A1[i0,j0]=A1[j0,i0]=b #A1[(j0-1):(j0+1),j0]=0; A1[(j0-1):(j0+1),i0]=c(a,1,a)
  # --- true alternative model (Ding model)
  A1[j0,i0]=sqrt(b); A1[j0,j0]=sqrt(1-b);
  A1[j0+1,i0]=a*sqrt(b); A1[j0+1,j0]=a*sqrt(1-b)
  A1[j0-1,i0]=a*sqrt(b); A1[j0-1,j0]=a*sqrt(1-b)
  # --- true alternative model (2D Ding model)
  # p=10; a=0.9; b=0.8; sigma=2; i0=23; j0=73
  # #M1=Mp(p^2)
  # A0=mdiag(p^2,c(1,a))
  # for(i in 1:(p^2)){
  #   if(i+p<=p^2) A0[i,i+p]=a
  #   if(i-p>=0) A0[i,i-p]=a
  # }
  # Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p^2))
  # 
  # A1=A0
  # A1[j0,i0]=sqrt(b); A1[j0,j0]=sqrt(1-b); 
  # A1[j0+1,i0]=a*sqrt(b); A1[j0+1,j0]=a*sqrt(1-b)
  # A1[j0-1,i0]=a*sqrt(b); A1[j0-1,j0]=a*sqrt(1-b)
  # p=100
  
  Sigma1=A1%*%t(A1)+sigma^2*diag(rep(1,p))
  C=t(base::chol(Sigma1))
  
  # best lambda (mu): covcov by formula (gaussian)
  Dmat=matrix(NA,5,5)
  dfi0j0=data.frame(i0=c(i0,i0-1,i0,i0+1,i0),j0=c(j0,j0,j0+1,j0,j0-1))
  for(i in 1:5){
    for(j in 1:5){
      i1=dfi0j0[i,"i0"]; j1=dfi0j0[i,"j0"]
      i2=dfi0j0[j,"i0"]; j2=dfi0j0[j,"j0"]
      Dmat[i,j]=ijthcovcov(C=C,M=M1,df=n,i1=i1,j1=j1,i2=i2,j2=j2)
    }
  }
  
  dvec=rep(0,5)
  Amat=matrix(c(2,a,a,a,a),5,1)
  bvec=2
  fit.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
  mu=fit.qp$solution
  mu
  ######### 1. simulation #########
  
  N=100 # number of replicates
  lmax10=lmax20=lmax30=rep(NA,N) # 1: coherence of S, 2: coherence of Sbar, 3: likelihood ratio, *0: null model 
  lmax11=lmax21=lmax31=rep(NA,N) # 1: coherence of S, 2: coherence of Sbar, 3: likelihood ratio, *1: alternative model 
  
  
  for(k in 1:N){
    print(k)
    # H1
    data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma1)
    Shat=cov(data)
    #------------------------------
    C=t(base::chol(Shat))
    for(i in 1:5){
      for(j in 1:5){
        i1=dfi0j0[i,"i0"]; j1=dfi0j0[i,"j0"]
        i2=dfi0j0[j,"i0"]; j2=dfi0j0[j,"j0"]
        Dmat[i,j]=ijthcovcov(C=C,M=M1,df=n,i1=i1,j1=j1,i2=i2,j2=j2)
      }
    }
    fit.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
    mu=fit.qp$solution
    #------------------------------
    lmax11[k]=max(abs(Shat)*(Sigma0==0))
    Sbar=sigmabar(data,Shat,mu=mu)$S
    lmax21[k]=maxd.r(abs(Sbar[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(Sigma0==0))
    lmax31[k]=0#lrt_a(data,Shat,i0,j0)
    # H0
    data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma0)
    Shat=cov(data)
    #------------------------------
    C=t(base::chol(Shat))
    for(i in 1:5){
      for(j in 1:5){
        i1=dfi0j0[i,"i0"]; j1=dfi0j0[i,"j0"]
        i2=dfi0j0[j,"i0"]; j2=dfi0j0[j,"j0"]
        Dmat[i,j]=ijthcovcov(C=C,M=M1,df=n,i1=i1,j1=j1,i2=i2,j2=j2)
      }
    }
    fit.qp=solve.QP(Dmat,dvec,Amat,bvec,meq=1)
    mu=fit.qp$solution
    #------------------------------
    lmax10[k]=max(abs(Shat)*(Sigma0==0))
    Sbar0=sigmabar(data,Shat,mu=mu)$S
    lmax20[k]=maxd.r(abs(Sbar0[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(Sigma0==0))
    lmax30[k]=0#lrt_a(data,Shat,i0,j0)
    
  }
  
  
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
  ######### 2. save results #########
  fname=paste0("mis/lsc_pa_ab_bar_b_",bi)
  save(list = ls(all=TRUE), file = paste0(fname,".RData"))
  
  png(paste0(fname,".png"))
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
  dev.off()
  
  return(c(EP11,EP21,EP31))
}

###### run simulation ######
p=100
M1=Mp(p)
#B=seq(0.2,0.7,length=10)
B=seq(200,2000,length=9)
EP1=EP2=EP3=rep(NA,length(B))
for(bi in 1:length(B)){
  print(bi)
  temp=simu(B,bi,M1)
  EP1[bi]=temp[1]
  EP2[bi]=temp[2]
  EP3[bi]=temp[3]
}

par(mfrow=c(1,1))
plot(B,EP1,type="b",ylim=c(0,1),xlab="sparse signal strength",ylab="power",main="power curve")
points(B,EP2,col="red",type="b")
#points(B,EP3,col="blue",type="b")
#legend("topleft",c("max of S","max of Sbar","LRT"),col=c("black","red","blue"),lty=c(1,1,1),cex=0.75)
legend("topleft",c("max of S","max of Sbar"),col=c("black","red"),lty=c(1,1),cex=0.75)


