source('utils.R')
#library(Matrix)
library(quadprog)
############## 1. Wishart Distribution ##############
# used to compute covariance of sample covariance matrix
ithrowkron=function(C,i0,j0){
  # compute the ith row of kronecker(C,C)
  # where i=(j0-1)p+i0, p is the dimension of C
  p=dim(C)[1]
  #i=(j0-1)*p+i0
  ithrow=c(as.matrix(C[j0,])%*%t(as.matrix(C[i0,])))
  return(ithrow)
}

# used to compute covariance of sample covariance matrix
# compute the covariance between sigma(i1,j1) and sigma(i2,j2)
ijthcovcov=function(C,M,df,i1,j1,i2,j2){
  # compute the (i,j)th element of cov(vec(S))
  # this matrix is too large to store the whole thing, that's why we calculate one by one
  # i=(j1-1)p+i1; j=(j2-1)p+i2; p is dimension of S
  # C is chol decomp of S, S=CC' i.e. C=t(chol(S))
  p=dim(C)[1]
  
  vi=as.matrix(ithrowkron(C,i1,j1))
  vj=as.matrix(ithrowkron(C,i2,j2))
  ijcovcov=t(vi)%*%M%*%vj
  return(df*ijcovcov[1,1])
}

# The strange matrix used in the covariance of sample covariance matrix
Mp=function(p){
  M1=Matrix::Matrix(0,p^2,p^2)
  for(idy in 1:(p^2)){
    idx=((idy+p-1)%/%p)+((idy-1)%%p)*p
    M1[idx,idy]=1
  }
  return(M1+Diagonal(p^2))
}

############## 2. Used for simulation ##############
# Ture Covariance Model
# a: nearest neighbour effect
# r: teleconnection effect 
# sf: standard deviation of factors f_i
# se: standard deviation of errors e_i
# p: dimension of data (the covariance matrix)
# n: sample size
# d: dimension of the structure of the model (1: features can be arranged in a line, neighbours are the left and right; 
#                                             2: two-dim grid, neighbours are left, right, up, down)
# M: strange matrix used in covariance of sample covariance matrix
TCM=function(para){
  a=para$a; r=para$r; sf=para$sf; se=para$se; p=para$p; n=para$n; i0=para$i0; j0=para$j0; d=para$d; M=para$M
  if(para$method=="unequalvar"){
    A0=mdiag.r(p,c(1,a)); Sigma0=sf^2*A0%*%t(A0)+se^2*diag(rep(1,p))
    A1=A0;
    A1[i0,j0]=A1[j0,i0]=r; Sigma1=sf^2*A1%*%t(A1)+se^2*diag(rep(1,p))
    return(list(S0=Sigma0,S1=Sigma1))
  }
  
  if(d==1){
    A0=mdiag.r(p,c(1,a)); Sigma0=sf^2*A0%*%t(A0)+se^2*diag(rep(1,p))
    A1=A0; 
    A1[j0,i0]=r; A1[j0,j0]=sqrt(1-r^2); 
    A1[j0+1,i0]=a*r; A1[j0+1,j0]=a*sqrt(1-r^2)
    A1[j0-1,i0]=a*r; A1[j0-1,j0]=a*sqrt(1-r^2)
  }else if(d==2){
    p=sqrt(para$p)
    A0=mdiag(p^2,c(1,a))
    for(i in 1:(p^2)){
      if(i+p<=p^2) A0[i,i+p]=a
      if(i-p>=0) A0[i,i-p]=a
    }
    Sigma0=A0%*%t(A0)+sigma^2*diag(rep(1,p^2))
    
    A1=A0
    A1[j0,i0]=r; A1[j0,j0]=sqrt(1-r^2); 
    A1[j0+1,i0]=a*r; A1[j0+1,j0]=a*sqrt(1-r^2)
    A1[j0-1,i0]=a*r; A1[j0-1,j0]=a*sqrt(1-r^2)
    A1[j0+p,i0]=a*r; A1[j0+p,j0]=a*sqrt(1-r^2)
    A1[j0-p,i0]=a*r; A1[j0-p,j0]=a*sqrt(1-r^2)
  }
  Sigma1=sf^2*A1%*%t(A1)+se^2*diag(rep(1,p))
  #C1=t(base::chol(Sigma1))
  return(list(S0=Sigma0,S1=Sigma1))
}

# Find optimal lambda (5 by 5)
# S: sample covariance used to create covariance of sample covariance (of which Dmat is a 5 by 5 submatrix)
# Amat: coeeficients for (center, up, right, down, left) for unbiasedness (=1)
# para$i0, j0: position of the teleconnection (to calculate lambda star)
LS5=function(S,para,Amat=NULL){
  n=para$n; i0=para$i0; j0=para$j0; d=para$d; M=para$M
  C1=t(base::chol(S))
  Dmat=matrix(NA,5,5)
  dfi0j0=data.frame(i0=c(i0,i0-1,i0,i0+1,i0),j0=c(j0,j0,j0+1,j0,j0-1))
  for(i in 1:5){
    for(j in 1:5){
      i1=dfi0j0[i,"i0"]; j1=dfi0j0[i,"j0"]
      i2=dfi0j0[j,"i0"]; j2=dfi0j0[j,"j0"]
      Dmat[i,j]=ijthcovcov(C=C1,M=M,df=n,i1=i1,j1=j1,i2=i2,j2=j2)
    }
  }
  dvec=rep(0,5)
  if(is.null(Amat)){
    a=para$a
    Amat=matrix(c(1,a,a,a,a),5,1)
  }
  bvec=1
  fit.qp=quadprog::solve.QP(Dmat,dvec,Amat,bvec,meq=1)
  # # solve quad-programming by hand
  # NUM=(2*solve(Dmat)%*%Amat)
  # DEN=(t(Amat)%*%solve(Dmat)%*%Amat)
  # ls5=NUM/DEN[1,1]
  return(list(ls5=fit.qp$solution,obj=fit.qp$value,Dmat=Dmat,Amat=Amat))
}

# Modify the sample covariance with lambda star filtering
# X: data matrix
# S: sample covariance matrix
# mu: lambda star
# bandwidth: filtering S starting away from "bandwidth"th diagonal
sigmabar=function(X,S=NULL,a=NULL,mu,bandwidth=2){
  n=dim(X)[1]
  p=dim(X)[2]
  if(is.null(S)){
    S=cov(X)
  }
  S1=S
  if(is.null(a)){
    ahat=mean(diag(S[-p,-1]))/2
    # ahat2=sqrt(mean(diag(S[-c(p-1,p),-c(1,2)])))
    # ahat22=mean(sqrt(diag(S[-c(p-1,p),-c(1,2)])),na.rm=T)
  }else{
    ahat=a
  }
  
  for(i in bandwidth:(p-1)){
    for(j in i:(p-1)){
      if(abs(i-j)>=2){
        vs=c(S[i,j],S[i-1,j],S[i,j+1],S[i+1,j],S[i,j-1])
        S1[i,j]=S1[j,i]=crossprod(vs,mu)
      }
    }
  }
  
  return(list(S=S1,a=ahat))
}



