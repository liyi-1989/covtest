# A test variance covariance matrix:
Sig <- matrix(c(1, .7, .6, 
               .7, 1, .4,
               .6, .4, 1), nrow = 3)
p<- dim(Sig)[1]
df<- 4
## EMPIRICAL (USING SIMULATION) ##
set.seed(901254)
# for replication
reps<- 100000 # number of obs in our sampling dist
W.empir <- matrix( nrow = reps, ncol = length( c(Sig) ) )
for(i in 1:reps)
  #W.empir[i, ] <- c(rwish(v = df, S = Sig))
  W.empir[i, ] <- c(rWishart(n=1,df = df, Sigma = Sig))
## THEORETICAL (USING EQUATION) ##
# The Cholesky decomposition of Sig:
C<- t( chol(Sig) )
# The strange M matrix:
M <- matrix(c(rep(c(rep(c(1, rep(0, times = p*p+(p-1))), times = p-1),
                              1, rep(0, times = p)),
                            times = p-1),
                        rep(c( 1, rep(0, times = p*p+(p-1)) ), times = p-1),
                        1),
                      nrow = p^2)
M


M1=Matrix(0,p^2,p^2)
for(idy in 1:(p^2)){
  idx=((idy+p-1)%/%p)+((idy-1)%%p)*p
  M1[idx,idy]=1
}
for(i in 1:(p^2)){
  M1[i,i]=M1[i,i]+1
}
M1


t0=proc.time()
t(1:p^2)%*%M1%*%(1:p^2)
proc.time()-t0

t0=proc.time()
t(1:p^2)%*%M%*%(1:p^2)
proc.time()-t0

#benchmark(t(1:p^2)%*%M1%*%(1:p^2),t(1:p^2)%*%M1%*%(1:p^2),replications = 1)

# And looks like the following:
M
W.theor <- {df * ( kronecker(Sig, Sig) + kronecker(C, C) %*% M %*% kronecker(t(C), t(C)) ) }
## COMPARING VAR/COV MATRICES ##
## 1. EMPIRICAL ##
round(var(W.empir), digits = 2)
## 2. THEORETICAL ##
round(W.theor, digits = 3)



#------------------------

X=mvrnorm(n, mu = rep(0,p), Sigma = Sigma1)
S=cov(X)
C=t(base::chol(S))
C=t(base::chol(Sigma1))

ijthcovcov(C=C,M=M1,df=n,i1=i0,j1=j0,i2=i0+1,j2=j0)
ijthcovcov(C=C,M=M1,df=n,i1=i0,j1=j0,i2=i0,j2=j0) 
