# law of large number of the coherence
source('utils_lscm.R')
library(mvtnorm)

para1=list(p=100,n=1000,a=0.9,r=0.8,sf=1,se=2,i0=NULL,j0=NULL,M=NULL,d=1,dd=5,py=NULL,method="equalvar",N=100) # one-dimensional model
para1$i0=round(para1$p*0.23); para1$j0=round(para1$p*0.77)

covmodel=TCM(para1) # generating true covariance model
S0=covmodel$S0; S1=covmodel$S1

printM(S0)

#n=para1$n; p=para1$p; i0=para1$i0; j0=para1$j0; N=para1$N 

Ns=(1:10)*100
N=length(Ns)
lmax10=lmax20=lmaxlln10=lmaxlln20=lmaxlln30=lmaxlln30=rep(0,N)

for(k in 1:N){
  n=Ns[k]
  para1=list(p=n,n=n,a=0.9,r=0.8,sf=1,se=2,i0=NULL,j0=NULL,M=NULL,d=1,dd=5,py=NULL,method="equalvar",N=100) # one-dimensional model
  para1$i0=round(para1$p*0.23); para1$j0=round(para1$p*0.77); p=para1$p
  S0=TCM(para1)$S0; 
  #S0=S0/diag(S0)
  #============== H0 ==============
  data = mvrnorm(n, mu = rep(0,n), Sigma = S0)
  Shat=cor(data)
  
  hatfit=a_efratio_hat(S=Shat,p=p,d=para1$d) # estimate a, se/sf
  a_hat=hatfit[1]
  efratio_hat=hatfit[2]
  lstar=LS5m(a=a_hat,ef_ratio=efratio_hat,d=para1$d)$ls5 # estimate lambda star
  #lstar=c(0.26,0.22,0.17,0.19,0.25)
  lstar=c(sin(30),cos(30),0,0,0)
  Sbar=sigmabar(S=Shat,mu=lstar) 
  lstar3=c(sin(45),cos(45),0,0,0)
  Sbar3=sigmabar(S=Shat,mu=lstar3) 
  
  lmax10[k]=max_hat=max(abs(Shat)*(S0==0))
  lmax20[k]=max_bar=maxd.r(abs(Sbar[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(S0==0))
  lmaxlln10[k]=sqrt(n/log(n))*lmax10[k]
  lmaxlln20[k]=sqrt(n/log(n))*lmax20[k]
  
  lmax30[k]=max_bar=maxd.r(abs(Sbar3[2:(p-1),2:(p-1)]),4)
  lmaxlln30[k]=sqrt(n/log(n))*lmax30[k]
  
  # mm_hat=abs(Shat[2:(p-1),2:(p-1)])
  # mm_bar=abs(Sbar[2:(p-1),2:(p-1)])

  #R0[k,]=c(a_hat,efratio_hat,sij_hat,sij_bar,s_hat_F,s_bar_F,s_hat_2,s_bar_2,max_hat,max_bar,mm_hat_ji[1,],mm_bar_ji[1,])
}

plot(lmax10,ylim=c(0,1))
points(lmax20,col="blue")
points(lmax30,col="red")

plot(lmaxlln10,ylim=c(0,6))
points(lmaxlln20,col="blue")
points(lmaxlln30,col="red")


data=scale(data)

ii=10;jj=50
data1=cbind(data[,ii]*data[,jj],data[,ii-1]*data[,jj],data[,ii]*data[,jj+1],data[,ii+1]*data[,jj],data[,ii]*data[,jj-1])

var(data1%*%lstar)



