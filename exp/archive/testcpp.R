library(rbenchmark)

d=1000
X=runif(d*d)
A=matrix(t(X)%*%X,d,d)
maxd.r(A,2)
maxd(A,2)


benchmark(maxd.r(A,2),maxd(A,2),replications = 1)


t0=proc.time()
maxd.r(A,2)
proc.time()-t0

t0=proc.time()
maxd(A,2)
proc.time()-t0

# mdiag

mdiag(10,c(1,2))
mdiag.r(10,c(1,2))

d=1000
v=rnorm(100)
benchmark(mdiag(d,v),mdiag.r(d,v),replications = 100)

t0=proc.time()
A=mdiag(d,v);
proc.time()-t0

t0=proc.time()
B=mdiag.r(d,v);
proc.time()-t0

#
t0=proc.time()
Shat=cor(data)
proc.time()-t0

t0=proc.time()
max(abs(Shat)*(S==0))
proc.time()-t0

t0=proc.time()
maxd(abs(Shat),3)
proc.time()-t0

# 

A=matrix(1:10000,100,100)
MnormM(A,h=1,type=4)
mnorm(A,h=1,type=1)

benchmark(MnormM(A,h=1,type=4),mnorm(A,h=1,type=1),replications = 1)

t0=proc.time()
M1=MnormM(A,h=1,type=4)
proc.time()-t0

t0=proc.time()
M2=mnorm(A,h=1,type=1)
proc.time()-t0

# 

t0=proc.time()
for(k in 1:10){
  data = mvrnorm(n, mu = rep(0,p), Sigma = S)
  Shat=cor(data)
  S0norm=mnorm(Shat,h=1,type=1)#MnormM(Shat,h=1,type=4,RM=F)
  lmax10[k]=max(abs(Shat)*(S==0))
  lmax20[k]=maxd(S0norm,4)#max(abs(S0norm)*(mdiag.r(p,c(1,1,1,1))==0))
}
proc.time()-t0


t0=proc.time()
for(k in 1:10){
  data = mvrnorm(n, mu = rep(0,p), Sigma = S)
  Shat=cor(data)
  S0norm=MnormM(Shat,h=1,type=4,RM=F)
  lmax10[k]=max(abs(Shat)*(S==0))
  lmax20[k]=maxd(S0norm,4)#max(abs(S0norm)*(mdiag.r(p,c(1,1,1,1))==0))
}
proc.time()-t0




