# long short covariance model
library(Rcpp)
source("utils.R")
Rcpp::sourceCpp('src/utils.cpp')
library(fields)
# H0: Sigma=B
# true cov
p=100; r=0.5; s=1; tele=0; sigma=1
Ip=diag(rep(1,p))
Type=13; CONST=10; amp=0.75
# V=genV(p=p,r=r,s=s,tele=tele)
# S=V%*%t(V)+diag(rep(1,p))

######## 1. H0 #########
# 10: H0 for cov
# 20: H0 for average cov
S=mdiag(p,c(1,r,r^2)) # true covariance under H0
n=1000 # sample size
N=100 # number of replicates
lmax10=lmax20=rep(NA,N) # coherence (3, S==0)
for(k in 1:N){
  data = mvrnorm(n, mu = rep(0,p), Sigma = S)
  Shat=cor(data)
  lmax10[k]=max(abs(Shat)*(S==0))
  S0norm=MnormM(Shat,h=1,type=Type,RM=F,Lmax = lmax10[k],const=CONST)#mnorm(Shat,h=1,type=Type)#
  lmax20[k]=maxd(S0norm,4)#max(abs(S0norm)*(mdiag.r(p,c(1,1,1,1))==0))
}
# 10
TS10=n*lmax10^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
alpha=0.05
CUT.simu10=quantile(TS10,1-alpha)
CUT.asym10=-log(8*pi)-2*log(log(1/(1-alpha)))#y_alpha

mean(TS10>CUT.simu10) # by quantile definition, it is always 0.05
ES10=mean(TS10>CUT.asym10) # Empirical size
ES10

# 20
TS20=n*lmax20^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
alpha=0.05
CUT.simu20=quantile(TS20,1-alpha)
CUT.asym20=-log(8*pi)-2*log(log(1/(1-alpha)))#y_alpha

mean(TS20>CUT.simu20) # by quantile definition, it is always 0.05
ES20=mean(TS20>CUT.asym20) # Empirical size
ES20


# 10 plot: pdf
plot(density(TS10),main="distribution (density) of the TS under H0")
abline(v=CUT.simu10,lty=1,col="red")
abline(v=CUT.asym10,lty=2,col="blue")
# 10 plot: cdf
x=seq(min(TS10),max(TS10),length=100)
y=exp(-exp(-x/2)/sqrt(8*pi))
plot(ecdf(TS10),ylim=c(0,1),main="distribution of the TS: empirical vs asymptotic under H0")
lines(x,y,col="red")

# 20 plot: pdf
plot(density(TS20),main="distribution (density) of the TS under H0")
abline(v=CUT.simu20,lty=1,col="red")
abline(v=CUT.asym20,lty=2,col="blue")
# 20 plot: cdf
x=seq(min(TS20),max(TS20),length=100)
y=exp(-exp(-x/2)/sqrt(8*pi))
plot(ecdf(TS20),ylim=c(0,1),main="distribution of the TS: empirical vs asymptotic under H0")
lines(x,y,col="red")

# 10 20 compare
plot(density(TS10),xlim=c(min(TS10,TS20),max(TS10,TS20)),main="H0 compare")
lines(density(TS20),col="red")
######## 2. H1 #########
S1=S;s=0.5
# S1[4:6,94:96]=S1[94:96,4:6]=2*(sqrt(p/n))*matrix(c(0,s^2,0,
#                        s^2,s,s^2,
#                        0,s^2,0),3,3)
S1[4:6,94:96]=S1[94:96,4:6]=amp*(sqrt(p/n))*matrix(c(s^2,s^2,s^2,
                                                   s^2,s,s^2,
                                                   s^2,s^2,s^2),3,3)
lmax11=lmax21=rep(NA,N) # coherence (3, S==0)
for(k in 1:N){
  data = mvrnorm(n, mu = rep(0,p), Sigma = S1)
  Shat=cor(data)
  lmax11[k]=max(abs(Shat)*(S==0))
  S1norm=MnormM(Shat,h=1,type=Type,RM=F,Lmax = lmax11[k],const=CONST) #mnorm(Shat,h=1,type=Type)
  lmax21[k]=max(abs(S1norm)*(mdiag.r(p,c(1,1,1,1))==0))
}
# 11
TS11=n*lmax11^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
EP11=mean(TS11>CUT.simu10) # Empirical power
EP11
# 21
TS21=n*lmax21^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
EP21=mean(TS21>CUT.simu20) # Empirical power
EP21


par(mfrow=c(1,2))
# 11 plot rejection region
plot(density(TS10),xlim=c(min(density(TS10)$x,density(TS11)$x),max(density(TS10)$x,density(TS11)$x)),
     main=paste("Lmax: power",EP11,"size",ES10))
lines(density(TS11),col="red")
abline(v=CUT.simu10,lty=1,col="red")
abline(v=CUT.asym10,lty=2,col="blue")

# 21 plot rejection region
plot(density(TS20),xlim=c(min(density(TS20)$x,density(TS21)$x),max(density(TS20)$x,density(TS21)$x)),
     main=paste("Average: power",EP21,"size",ES20))
lines(density(TS21),col="red")
abline(v=CUT.simu20,lty=1,col="red")
abline(v=CUT.asym20,lty=2,col="blue")








# 11 plot rejection region
plot(density(lmax10))
lines(density(lmax11),col="red")


plot(density(lmax20))
lines(density(lmax21),col="red")

