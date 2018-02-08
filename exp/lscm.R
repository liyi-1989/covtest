# long short covariance model 
# bar: average surroundings
#library(Rcpp)
source('utils_lscm.R')
#Rcpp::sourceCpp('src/utils.cpp')
#library(fields)
library(mvtnorm)
#library(matrixcalc)
#library(tikzDevice)
####### 0. True Covariance Model ######
# parameters (N: number of replication)
para1=list(p=100,n=1000,a=0.9,r=0.8,sf=1,se=2,i0=NULL,j0=NULL,M=NULL,d=1,method="equalvar",N=100) # one-dimensional model
para1$i0=round(para1$p*0.23); para1$j0=round(para1$p*0.77)
#para1$M=Mp(para1$p)
#para2=para1; para2$d=2 # two-dimensional model 

covmodel=TCM(para1) # generating true covariance model
S0=covmodel$S0; S1=covmodel$S1
ls5fit=LS5(S1,para1)
ls1=ls5fit$ls5 # find optimal lambda star #Dmat=ls5fit$Dmat; Amat=ls5fit$Amat

####### 0.1 small simulation #######
# # double check optimal variance, i.e. objective (if variance reduced significantly)
# Variance=t(ls1)%*%Dmat%*%ls1
# ls0=c(1,0,0,0,0)
# Variance0=t(ls0)%*%Dmat%*%ls0
# Variance/Variance0
# # # check constraints
# # t(ls1)%*%Amat
# # t(ls0)%*%Amat

####### 0.2 simulation: lambda star #######
# The lambda star for (1) 5 by 5 case; (2) in the exact factor model with a; have manually derive solution
# Even though we can do that by LS5 function for all, due to the C1=chol() step, 
# we cannot find C1 explicitly in terms of the parameters (any nearest neighour cov of cov, too time consuming, and complex when not exact factor model).
# So we need:
# (1). Double check our manually derive lambda star for the 5 by 5 case
# (2). study the behaviour of the lambda star i.e. the variance-reduction combination of the surroundings. 
# Main work: derive Dmat manually as a function of a,r,sf,se (true parameters)

# ls5mfit=LS5m(para1) # lambda star by hand in 1-dim model
# ls1m=ls5mfit$ls5
# c(ls1,ls1m)
# ls5mfit=LS5m(para2) # lambda star by hand in 2-dim model
# ls2m=ls5mfit$ls5
# c(ls2,ls2m)

####### 1. Estimation and Test Simulation ######
# 1. estimate a
# 2. estimate r
# 3. estimate covariance matrix
n=para1$n; p=para1$p; i0=para1$i0; j0=para1$j0; N=para1$N 

R0=as.data.frame(matrix(NA,N,14))
colnames(R0)=c("a_hat","efratio_hat","sij_hat","sij_bar","s_hat_F","s_bar_F","s_hat_2","s_bar_2","max_hat","max_bar","max_hat_i","max_hat_j","max_bar_i","max_bar_j")
R1=R0

t0=proc.time()
set.seed(1)
for(k in 1:N){
  #============== H0 ==============
  data = mvrnorm(n, mu = rep(0,p), Sigma = S0)
  Shat=cov(data)
  
  hatfit=a_efratio_hat(S=Shat,d=para1$d) # estimate a, se/sf
  a_hat=hatfit[1]
  efratio_hat=hatfit[2]
  lstar=LS5m(a=a_hat,ef_ratio=efratio_hat,d=para1$d)$ls5 # estimate lambda star
  #lstar=ls1
  Sbar=sigmabar(S=Shat,mu=lstar) # sigmabar
  
  sij_hat=Shat[i0,j0]
  sij_bar=Sbar[i0,j0]
  s_hat_F=base::norm(Shat-S0,"F")
  s_bar_F=base::norm(Sbar-S0,"F")
  s_hat_2=base::norm(Shat-S0,"2")
  s_bar_2=base::norm(Sbar-S0,"2")
  max_hat=max(abs(Shat)*(S0==0))
  max_bar=maxd.r(abs(Sbar[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(S0==0))
  mm_hat=abs(Shat[2:(p-1),2:(p-1)])
  mm_hat_ji=which(mm_hat == maxd.r(mm_hat,4), arr.ind = TRUE)
  mm_bar=abs(Sbar[2:(p-1),2:(p-1)])
  mm_bar_ji=which(mm_bar == maxd.r(mm_bar,4), arr.ind = TRUE)
  
  R0[k,]=c(a_hat,efratio_hat,sij_hat,sij_bar,s_hat_F,s_bar_F,s_hat_2,s_bar_2,max_hat,max_bar,mm_hat_ji[1,],mm_bar_ji[1,])
  #============== H1 ==============
  data = mvrnorm(n, mu = rep(0,p), Sigma = S1)
  Shat=cov(data)
  
  hatfit=a_efratio_hat(S=Shat,d=para1$d) # estimate a, se/sf
  a_hat=hatfit[1]
  efratio_hat=hatfit[2]
  lstar=LS5m(a=a_hat,ef_ratio=efratio_hat,d=para1$d)$ls5 # estimate lambda star
  #lstar=ls1
  Sbar=sigmabar(S=Shat,mu=lstar) # sigmabar
  
  sij_hat=Shat[i0,j0]
  sij_bar=Sbar[i0,j0]
  s_hat_F=base::norm(Shat-S0,"F")
  s_bar_F=base::norm(Sbar-S0,"F")
  s_hat_2=base::norm(Shat-S0,"2")
  s_bar_2=base::norm(Sbar-S0,"2")
  max_hat=max(abs(Shat)*(S0==0))
  max_bar=maxd.r(abs(Sbar[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(S0==0))
  mm_hat=abs(Shat[2:(p-1),2:(p-1)])
  mm_hat_ji=which(mm_hat == maxd.r(mm_hat,4), arr.ind = TRUE)
  mm_bar=abs(Sbar[2:(p-1),2:(p-1)])
  mm_bar_ji=which(mm_bar == maxd.r(mm_bar,4), arr.ind = TRUE)
  
  R1[k,]=c(a_hat,efratio_hat,sij_hat,sij_bar,s_hat_F,s_bar_F,s_hat_2,s_bar_2,max_hat,max_bar,mm_hat_ji[1,],mm_bar_ji[1,])
}
t1=proc.time()-t0
t1

#a_true=para$a
sij_0=S0[i0,j0]
sij_1=S1[i0,j0]
####### 2. Collect and Analyze results ######
# apply(R0,2,mean)
# apply(R1,2,mean)
# apply(R0,2,sd)
# apply(R1,2,sd)
lmax10=R0[,"max_hat"]
lmax20=R0[,"max_bar"]
lmax11=R1[,"max_hat"]
lmax21=R1[,"max_bar"]

alpha=0.05
mean(lmax11>quantile(lmax10,1-alpha))
mean(lmax21>quantile(lmax20,1-alpha))
############ cut-off ############
TS10=n*lmax10^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
CUT.simu10=quantile(TS10,1-alpha) # CUT.asym10=-log(8*pi)-2*log(log(1/(1-alpha)))#y_alpha
TS20=n*lmax20^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
CUT.simu20=quantile(TS20,1-alpha) # type I error: mean(TS10>CUT.simu10) # by quantile definition, it is always 0.05
############ type II error: power ############
TS11=n*lmax11^2-4*log(p)+log(log(p)) # coherence of S
EP11=mean(TS11>CUT.simu10) # Empirical power
EP11
TS21=n*lmax21^2-4*log(p)+log(log(p)) # coherence of Sbar
EP21=mean(TS21>CUT.simu20) # Empirical power
EP21

# filepath="./mis/"
# filename="simu"
# save(R0,R1,sij_0,sij_1,EP11,EP21,para1,file=paste0(filepath,filename,".RData"))


###### 3. Plots ######
###### 3.1 Plots for Covariance Model ######
####### make plot for 1-dim model #######
n1=31
ai=li=li2=li3=rep(0,n1)
p1t=para1
for(i in 1:n1){
  p1t$a=ai[i]=(i-1)/30
  p1t$se=1; li[i]=LS5m(p1t)$ls5[2]
  p1t$se=2; li2[i]=LS5m(p1t)$ls5[2]
  p1t$se=3; li3[i]=LS5m(p1t)$ls5[2]
  
}
tikz("./fig/lambdastar1d.tex", width = 3.25, height = 3.25)
plot(ai,1-4*ai*li,type="l",xlim=c(0,1),ylim=c(0,1),xlab="$a$",ylab="$\\lambda$",col="blue",lwd=1)
lines(ai,li,col="orange",lwd=1,lty=2)
lines(ai,1-4*ai*li2,col="blue",lwd=2)
lines(ai,li2,col="orange",lwd=2,lty=2)
lines(ai,1-4*ai*li3,col="blue",lwd=3)
lines(ai,li3,col="orange",lwd=3,lty=2)
legend("topright",c("$\\lambda$ for center","$\\lambda$ for others"),col=c("blue","orange"),lty=c(1,2),cex=0.8)
dev.off()
####### make plot for 2-dim model #######
n1=31
ai=li=li2=li3=rep(0,n1)
p1t=para2
for(i in 1:n1){
  p1t$a=ai[i]=(i-1)/30
  p1t$se=1; li[i]=LS5m(p1t)$ls5[2]
  p1t$se=2; li2[i]=LS5m(p1t)$ls5[2]
  p1t$se=3; li3[i]=LS5m(p1t)$ls5[2]
  
}
tikz("./fig/lambdastar2d.tex", width = 3.25, height = 3.25)
plot(ai,1-4*ai*li,type="l",xlim=c(0,1),ylim=c(0,1),xlab="$a$",ylab="$\\lambda$",col="blue",lwd=1)
lines(ai,li,col="orange",lwd=1,lty=2)
lines(ai,1-4*ai*li2,col="blue",lwd=2)
lines(ai,li2,col="orange",lwd=2,lty=2)
lines(ai,1-4*ai*li3,col="blue",lwd=3)
lines(ai,li3,col="orange",lwd=3,lty=2)
legend("topright",c("$\\lambda$ for center","$\\lambda$ for others"),col=c("blue","orange"),lty=c(1,2),cex=0.8)
dev.off()
###### 3.2 Plots for Estimation Problem ######
######## make plot for estimation results ######## 
# Estimation of a
tikz("./fig/est_a_h0.tex", width = 3.25, height = 3.25)
hist(R0[,"a_bar"],ylim=c(0,20), main="$H_0$: estimation of a",xlab="a",ylab="density",breaks=25)
lines(density(R0[,"a_bar"]),lty=1,col="blue")
abline(v=para1$a,col="red",lty=4)
abline(v=mean(R0[,"a_bar"]),col="blue",lty=4)
dev.off()

tikz("./fig/est_a_h1.tex", width = 3.25, height = 3.25)
hist(R1[,"a_bar"],ylim=c(0,20), main="$H_1$: estimation of a",xlab="a",ylab="density",breaks=25)
lines(density(R1[,"a_bar"]),lty=1,col="blue")
abline(v=para1$a,col="red",lty=4)
abline(v=mean(R1[,"a_bar"]),col="blue",lty=4)
dev.off()

# Estimation of S_{i0j0}
tikz("./fig/h0h1dis.tex", width = 3.25, height = 3.25)
plot(density(R0[,"sij_hat"]),col="orange",lty=2,lwd=1,xlab="$S_{i_0j_0}$",main="$H_0$ and $H_1$ distributions",ylim=c(0,4),xlim=range(c(R0[,"sij_hat"],R0[,"sij_bar"],R1[,"sij_hat"],R1[,"sij_bar"])))
lines(density(R0[,"sij_bar"]),col="blue",lty=2,lwd=1)
lines(density(R1[,"sij_hat"]),col="orange",lty=1,lwd=1)
lines(density(R1[,"sij_bar"]),col="blue",lty=1,lwd=1)
abline(v=S0[i0,j0],col="red",lty=1)
abline(v=S1[i0,j0],col="red",lty=1)
legend("topright",c("$H_0$: $\\hat{S}_{i_0j_0}$","$H_0$: $\\bar{S}_{i_0j_0}$","$H_1$: $\\hat{S}_{i_0j_0}$","$H_1$: $\\bar{S}_{i_0j_0}$"),col=c("orange","blue","orange","blue"),lty=c(2,2,1,1),lwd=c(1,1,1,1),cex=0.5)
dev.off()

# Estimation of S
tikz("./fig/est_s.tex", width = 3.25, height = 3.25)
plot(density(R0[,"s_hat"]),col="orange",lty=2,lwd=1,xlab="Error in Frobenius norm",main="Estimation of Covariance",ylim=c(0,2),xlim=range(c(R0[,"s_hat"],R0[,"s_bar"],R1[,"s_hat"],R1[,"s_bar"])))
lines(density(R0[,"s_bar"]),col="blue",lty=2,lwd=1)
lines(density(R1[,"s_hat"]),col="orange",lty=1,lwd=1)
lines(density(R1[,"s_bar"]),col="blue",lty=1,lwd=1)
# abline(v=S0[i0,j0],col="red",lty=1)
# abline(v=S1[i0,j0],col="red",lty=1)
legend("topright",c("$H_0$: $\\hat{S}$","$H_0$: $\\bar{S}$","$H_1$: $\\hat{S}$","$H_1$: $\\bar{S}$"),col=c("orange","blue","orange","blue"),lty=c(2,2,1,1),lwd=c(1,1,1,1),cex=0.5)
dev.off()
###### 3.3 Plots for Test Problem ######
###### make plot of distribution of test statistics ######
par(mfrow=c(1,1))
tikz("./fig/test_stat_h0h1.tex", width = 3.25, height = 3.25)
plot(density(TS10),col="orange",lty=2,lwd=1,xlab="TS",main="$H_0$ and $H_1$ distributions",xlim=range(c(TS10,TS11,TS20,TS21)))
lines(density(TS20),col="blue",lty=2,lwd=1)
lines(density(TS11),col="orange",lty=1,lwd=1)
lines(density(TS21),col="blue",lty=1,lwd=1)
abline(v=CUT.simu10,lty=1,col="grey")
abline(v=CUT.simu20,lty=1,col="grey")
legend("topright",c("$H_0$: TS with $\\hat{S}$","$H_0$: TS with $\\bar{S}$","$H_1$: TS with $\\hat{S}$","$H_1$: TS with $\\bar{S}$"),col=c("orange","blue","orange","blue"),lty=c(2,2,1,1),lwd=c(1,1,1,1),cex=0.5)
dev.off()

apply(cbind(TS10,TS11,TS20,TS21),2,mean) # check variance reduction 
apply(cbind(TS10,TS11,TS20,TS21),2,sd)

###### make plot of distribution of k-coherence (raw test statistics) ######
tikz("./fig/coherence_h0h1.tex", width = 3.25, height = 3.25)
plot(density(lmax10),col="orange",lty=2,lwd=1,xlab="Lmax",main="$H_0$ and $H_1$ distributions",xlim=range(c(lmax10,lmax11,lmax20,lmax21)))
lines(density(lmax20),col="blue",lty=2,lwd=1)
lines(density(lmax11),col="orange",lty=1,lwd=1)
lines(density(lmax21),col="blue",lty=1,lwd=1)
abline(v=quantile(lmax10,1-alpha),lty=1,col="grey")
abline(v=quantile(lmax20,1-alpha),lty=1,col="grey")
legend("topright",c("$H_0$: TS with $\\hat{S}$","$H_0$: TS with $\\bar{S}$","$H_1$: TS with $\\hat{S}$","$H_1$: TS with $\\bar{S}$"),col=c("orange","blue","orange","blue"),lty=c(2,2,1,1),lwd=c(1,1,1,1),cex=0.5)
dev.off()

apply(cbind(lmax10,lmax11,lmax20,lmax21),2,mean) # check variance reduction 
apply(cbind(lmax10,lmax11,lmax20,lmax21),2,sd)
apply(cbind(lmax10,lmax11,lmax20,lmax21)^2,2,sd)

mm=abs(Sbar0[2:(p-1),2:(p-1)])
which(mm == maxd.r(mm,4), arr.ind = TRUE)
mm=abs(Sbar[2:(p-1),2:(p-1)])
which(mm == maxd.r(mm,4), arr.ind = TRUE)

mm=abs(Shat[2:(p-1),2:(p-1)])
which(mm == maxd.r(mm,4), arr.ind = TRUE)
mm=abs(Sbar[2:(p-1),2:(p-1)])
which(mm == maxd.r(mm,4), arr.ind = TRUE)