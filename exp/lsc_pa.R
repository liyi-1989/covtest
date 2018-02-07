# long short covariance model 
# (with position alternative):
# H1: only changes sparse signal locations
library(Rcpp)
source("utils.R")
Rcpp::sourceCpp('src/utils.cpp')
library(fields)
library(mvtnorm)
library(matrixcalc)
############## functions ##############
pa=function(data,shortband=2, longband=1){
  n=dim(data)[1]
  p=dim(data)[2]
  Shat=cov(data)
  # l0
  S0=banding(Shat,k=shortband)
  l0=sum(dmvnorm(data, mean=rep(0,p), sigma=S0, log=T))
  # l1
  L1=matrix(NA,p,p) # loglikelihood matrix under H1
  M1=NA # max of loglikelihood over parameter space Theta_1
  for(i in 1:p){
    for(j in 1:p){
      if(abs(i-j)>shortband+longband+2 & longband<i & i<p-longband & longband<j & j<p-longband){
        S1=S0
        idi=(i-longband):(i+longband)
        idj=(j-longband):(j+longband)  
        S1[idi,idj]=Shat[idi,idj]
        S1[idj,idi]=Shat[idj,idi]
        if(is.positive.semi.definite(S1)){
          L1[i,j]=sum(dmvnorm(data, mean=rep(0,p), sigma=S1, log=T))
        }
        
      }
    }
  }
  M1=max(L1,na.rm = T)
  l1=max(M1,l0,na.rm = T)
  # transformed loglikelihood ratio
  return(2*(l1-l0))
}

############## code ##############

# H0: Sigma=B
# true cov
p=100; r=0.5; s=1; tele=0; sigma=1
Ip=diag(rep(1,p))
Type=13; CONST=10; amp=0.75

######## 1. H0 #########
# 10: H0 for cov
# 20: H0 for average cov
S=mdiag(p,c(1,r,r^2)) # true covariance under H0
n=1000 # sample size
N=100 # number of replicates
lmax10=lmax30=rep(NA,N) # coherence (3, S==0)
shortband=2;longband=1
t0=proc.time()
for(k in 1:N){
  print(k)
  data = mvrnorm(n, mu = rep(0,p), Sigma = S)
  Shat=cor(data)
  lmax10[k]=max(abs(Shat)*(S==0))
  lmax30[k]=pa(data,shortband=shortband, longband=longband) # known bandwidth of true covariance. Only long position unknown
}
t1=proc.time()-t0
t1
save(lmax10,lmax30,file="mis/pa_lmax10_30.RData")
# 10
TS10=n*lmax10^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
alpha=0.05
CUT.simu10=quantile(TS10,1-alpha)
CUT.asym10=-log(8*pi)-2*log(log(1/(1-alpha)))#y_alpha

TS30=lmax30
alpha=0.05
CUT.simu30=quantile(TS30,1-alpha)
CUT.asym30= qchisq(1-alpha,df=(1+p-shortband)*(p-shortband)/2)# chi-sq 

# size 10
mean(TS10>CUT.simu10) # by quantile definition, it is always 0.05
ES10=mean(TS10>CUT.asym10) # Empirical size
ES10
# size 30
mean(TS30>CUT.simu30) # by quantile definition, it is always 0.05
ES30=mean(TS30>CUT.asym30) # Empirical size
ES30

# 10 plot: pdf
plot(density(TS10),main="distribution (density) of the TS under H0")
abline(v=CUT.simu10,lty=1,col="red")
abline(v=CUT.asym10,lty=2,col="blue")
# 10 plot: cdf
x=seq(min(TS10),max(TS10),length=100)
y=exp(-exp(-x/2)/sqrt(8*pi))
plot(ecdf(TS10),ylim=c(0,1),main="distribution of the TS: empirical vs asymptotic under H0")
lines(x,y,col="red")


######## 2. H1 #########

S1=S;s=0.5
# S1[4:6,94:96]=S1[94:96,4:6]=2*(sqrt(p/n))*matrix(c(0,s^2,0,
#                        s^2,s,s^2,
#                        0,s^2,0),3,3)
S1[4:6,94:96]=S1[94:96,4:6]=amp*(sqrt(p/n))*matrix(c(s^2,s^2,s^2,
                                                     s^2,s,s^2,
                                                     s^2,s^2,s^2),3,3)
lmax11=lmax31=rep(NA,N) # coherence (3, S==0)
t0=proc.time()
for(k in 1:N){
  print(k)
  data = mvrnorm(n, mu = rep(0,p), Sigma = S1)
  Shat=cor(data)
  lmax11[k]=max(abs(Shat)*(S==0))
  lmax31[k]=pa(data,shortband=shortband, longband=longband)
}
t1=proc.time()-t0
t1
save(lmax11,lmax31,file="mis/pa_lmax11_31_0_75.RData")
load("mis/pa_lmax11_31_0_75.RData")
# 11
TS11=n*lmax11^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
EP11=mean(TS11>CUT.simu10) # Empirical power
EP11
# 31
TS31=lmax31 # test statistics of likelihood ratio test
EP31=mean(TS31>CUT.simu30) # Empirical power
EP31

par(mfrow=c(1,2))
# 11 plot rejection region
plot(density(TS10),xlim=c(min(density(TS10)$x,density(TS11)$x),max(density(TS10)$x,density(TS11)$x)),
     main=paste("Lmax: power",EP11,"size",ES10))
lines(density(TS11),col="red")
abline(v=CUT.simu10,lty=1,col="red")
abline(v=CUT.asym10,lty=2,col="blue")

# 21 plot rejection region
plot(density(TS30),xlim=c(min(density(TS30)$x,density(TS31)$x),max(density(TS30)$x,density(TS31)$x)),
     main=paste("Position: power",EP31,"size",ES30))
lines(density(TS31),col="red")
abline(v=CUT.simu30,lty=1,col="red")
abline(v=CUT.asym30,lty=2,col="blue")

####################
##### amp=1 #####
####################
amp=1
S1=S;s=0.5
S1[4:6,94:96]=S1[94:96,4:6]=amp*(sqrt(p/n))*matrix(c(s^2,s^2,s^2,
                                                     s^2,s,s^2,
                                                     s^2,s^2,s^2),3,3)
lmax11=lmax31=rep(NA,N) # coherence (3, S==0)
for(k in 1:N){
  print(k)
  data = mvrnorm(n, mu = rep(0,p), Sigma = S1)
  Shat=cor(data)
  lmax11[k]=max(abs(Shat)*(S==0))
  lmax31[k]=pa(data,shortband=shortband, longband=longband)
}
save(lmax11,lmax31,file="mis/pa_lmax11_31_1_00.RData")
load("mis/pa_lmax11_31_1_00.RData")
# 11
TS11=n*lmax11^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
EP11=mean(TS11>CUT.simu10) # Empirical power
EP11
# 31
TS31=lmax31 # test statistics of likelihood ratio test
EP31=mean(TS31>CUT.simu30) # Empirical power
EP31

par(mfrow=c(1,2))
# 11 plot rejection region
plot(density(TS10),xlim=c(min(density(TS10)$x,density(TS11)$x),max(density(TS10)$x,density(TS11)$x)),
     main=paste("Lmax: power",EP11,"size",ES10))
lines(density(TS11),col="red")
abline(v=CUT.simu10,lty=1,col="red")
abline(v=CUT.asym10,lty=2,col="blue")

# 21 plot rejection region
plot(density(TS30),xlim=c(min(density(TS30)$x,density(TS31)$x),max(density(TS30)$x,density(TS31)$x)),
     main=paste("Position: power",EP31,"size",ES30))
lines(density(TS31),col="red")
abline(v=CUT.simu30,lty=1,col="red")
abline(v=CUT.asym30,lty=2,col="blue")
####################
##### amp=1.1 #####
####################
amp=1.1
S1=S;s=0.5
S1[4:6,94:96]=S1[94:96,4:6]=amp*(sqrt(p/n))*matrix(c(s^2,s^2,s^2,
                                                     s^2,s,s^2,
                                                     s^2,s^2,s^2),3,3)
lmax11=lmax31=rep(NA,N) # coherence (3, S==0)
for(k in 1:N){
  print(k)
  data = mvrnorm(n, mu = rep(0,p), Sigma = S1)
  Shat=cor(data)
  lmax11[k]=max(abs(Shat)*(S==0))
  lmax31[k]=pa(data,shortband=shortband, longband=longband)
}
save(lmax11,lmax31,file="mis/pa_lmax11_31_1_10.RData")
load("mis/pa_lmax11_31_1_10.RData")
# 11
TS11=n*lmax11^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
EP11=mean(TS11>CUT.simu10) # Empirical power
EP11
# 31
TS31=lmax31 # test statistics of likelihood ratio test
EP31=mean(TS31>CUT.simu30) # Empirical power
EP31

par(mfrow=c(1,2))
# 11 plot rejection region
plot(density(TS10),xlim=c(min(density(TS10)$x,density(TS11)$x),max(density(TS10)$x,density(TS11)$x)),
     main=paste("Lmax: power",EP11,"size",ES10))
lines(density(TS11),col="red")
abline(v=CUT.simu10,lty=1,col="red")
abline(v=CUT.asym10,lty=2,col="blue")

# 21 plot rejection region
plot(density(TS30),xlim=c(min(density(TS30)$x,density(TS31)$x),max(density(TS30)$x,density(TS31)$x)),
     main=paste("Position: power",EP31,"size",ES30))
lines(density(TS31),col="red")
abline(v=CUT.simu30,lty=1,col="red")
abline(v=CUT.asym30,lty=2,col="blue")
####################
##### amp=1.25 #####
####################
amp=1.25
S1=S;s=0.5
S1[4:6,94:96]=S1[94:96,4:6]=amp*(sqrt(p/n))*matrix(c(s^2,s^2,s^2,
                                                     s^2,s,s^2,
                                                     s^2,s^2,s^2),3,3)
lmax11=lmax31=rep(NA,N) # coherence (3, S==0)
for(k in 1:N){
  print(k)
  data = mvrnorm(n, mu = rep(0,p), Sigma = S1)
  Shat=cor(data)
  lmax11[k]=max(abs(Shat)*(S==0))
  lmax31[k]=pa(data,shortband=shortband, longband=longband)
}
save(lmax11,lmax31,file="mis/pa_lmax11_31_1_25.RData")
load("mis/pa_lmax11_31_1_25.RData")
# 11
TS11=n*lmax11^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
EP11=mean(TS11>CUT.simu10) # Empirical power
EP11
# 31
TS31=lmax31 # test statistics of likelihood ratio test
EP31=mean(TS31>CUT.simu30) # Empirical power
EP31

par(mfrow=c(1,2))
# 11 plot rejection region
plot(density(TS10),xlim=c(min(density(TS10)$x,density(TS11)$x),max(density(TS10)$x,density(TS11)$x)),
     main=paste("Lmax: power",EP11,"size",ES10))
lines(density(TS11),col="red")
abline(v=CUT.simu10,lty=1,col="red")
abline(v=CUT.asym10,lty=2,col="blue")

# 21 plot rejection region
plot(density(TS30),xlim=c(min(density(TS30)$x,density(TS31)$x),max(density(TS30)$x,density(TS31)$x)),
     main=paste("Position: power",EP31,"size",ES30))
lines(density(TS31),col="red")
abline(v=CUT.simu30,lty=1,col="red")
abline(v=CUT.asym30,lty=2,col="blue")

