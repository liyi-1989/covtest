library(MASS) # use mvnorm to generate multivariate normal
library(Matrix) # use "dgCMatrix" matrix class, so that to make plot of a matrix like hand written order
library(QRM)
#library(moments) # calculate higher order moments, skewness, kurtosis. not useful
#library(pROC) # make ROC curve, not very useful currently

# Generate V in the factor model
genV=function(p=100,r=0.5,s=0.2,sp,tele=F){
  # p: number of features, assume equal number of latent factors. So V is p by p
  # r: neighbour kernel bump
  # s: tele connection bump scale
  # sp positions of teleconnection, like c(5,90)
  # tele: teleconnection signal. 1: s*(r,1,r) -1: s*(-r,1,-r)
  if(missing(sp)){
    sp=c(0.05*p,0.9*p)
  }
  V=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(abs(i-j)==1){
        V[i,j]=r
      }else if(i-j==0){
        V[i,j]=1
      }
    }
  }
  if(tele==-1){
    V[sp[2]:(sp[2]+2),sp[1]]=c(-r,1,-r)*s
  }else if(tele==1){
    V[sp[2]:(sp[2]+2),sp[1]]=c(r,1,r)*s
  }
  return(V)
}


# print a matrix with heatmap, in a hand written order, left->right, up->down (no need to rotate)
printM=function(M){
  image(as(M, "dgCMatrix"),main=deparse(substitute(M)))
}

# plot the density of a sequence and its normal approx with same mean & sd
plotDen=function(res){
  plot(density(res),col="blue",main=paste("Density of",deparse(substitute(res))))
  x=seq(min(res),max(res),length=100)
  y=dnorm(x,mean=mean(res), sd=sd(res))
  lines(x,y, type="l", lwd=1,col="red")
  abline(v=mean(res),lty=2,col="orange")
  text(mean(res)+2*sd(res),range(density(res)$y)[2]/2,paste("mean:",round(mean(res),4)),col="green")
  # fit Gumbel distribution
  y2=dGumbel(x, mu = mean(res)-0.5772*sd(res)*sqrt(6)/pi, sigma = sd(res)*sqrt(6)/pi, log = FALSE) 
  lines(x,y2, type="l", lwd=1,col="pink")
  #y3=dGumbel(x, mu = -log(sqrt(8*pi))+2*0.5772, sigma = 2*pi/sqrt(6), log = FALSE) 
  y3=dGumbel(x, mu = -log(sqrt(8*pi)), sigma = 2, log = FALSE) 
  lines(x,y3, type="l", lwd=2,col="black")
}

# Calculate "different norms" of the h-nearest-neighbour small block of a matrix M
MnormM=function(M,h=2,type=1,RM=F,Lmax,const){
  # h is the bandwidth
  p=dim(M)[1]
  normM=M
  
  # for(i in h:(p-h+1)){
  #   for(j in h:(p-h+1)){
  for(i in 1:(p)){
    for(j in 1:(p)){    
      il=max(i-h,1)
      iu=min(i+h,p)
      jl=max(j-h,1)
      ju=min(j+h,p)
      n0=(iu-il)*(ju-jl)
      M0=M[il:iu,jl:ju]
      if(RM){
        M0=M0-mean(M0) ########
      }
      
      if(type==1){
        normM[i,j]=base::norm(M0,"2") 
      }else if(type==2){
        normM[i,j]=base::norm(M0,"F")
      }else if(type==3){
        normM[i,j]=mean(M0)
      }else if(type==4){
        normM[i,j]=mean(abs(M0))#sum(M0)/((iu-il+1)*(ju-jl+1))
      }else if(type==5){
        normM[i,j]=log((prod(abs(M0)))^{1/n0})
      }else if(type==6){
        normM[i,j]=sort(M0,decreasing = TRUE)[2]
      }else if(type==7){
        m0=c(M0[2,2],M0[2,1],M0[2,3],M0[1,2],M0[3,2])
        normM[i,j]=mean(m0)
      }else if(type==8){
        m0=c(M0[2,2],M0[2,1],M0[2,3],M0[1,2],M0[3,2])
        normM[i,j]=abs(prod(m0))
      }else if(type==9){
        normM[i,j]=exp(-abs(mean(M0)))
      }else if(type==10){
        normM[i,j]=sd(c(M0))
      }else if(type==11){
        normM[i,j]=skewness(c(M0))
      }else if(type==12){
        normM[i,j]=kurtosis(c(M0))
      }else if(type==13){
        # L=max(Lmax-max(abs(M0)),0.001)
        # thr=max(abs(M0))*L/(L+(Lmax-L)*const)
        # normM[i,j]=M[i,j]*(abs(M[i,j])>thr)
        normM[i,j]=M[i,j]*mean(abs(M0))/Lmax
      }
      
    }
  }
  
  return(normM)
}


# variance of a covariance elements (used in Tony Cai's Adaptive threshouding) 
thetahat=function(X){
  n=dim(X)[1]
  Sigma_hat=cov(X)
  X1=apply(X,2,function(x){return (x-mean(x))}) # remove column mean version
  Sigma_var=((t(X1^2)%*%(X1^2)))/n-2*Sigma_hat*(t(X1)%*%X1)/n+(Sigma_hat)^2
  return(Sigma_var)
}

# adaptive threshoulding estimator

thresholding=function(X,type,const,threshold){
  n=dim(X)[1]
  p=dim(X)[2]
  Shat=cov(X)
  if(type=="hard"){
    ST=Shat*(Shat>threhold)
  }else if(type=="soft"){
    ST=sign(Shat)*pmax(Shat-threshold,0)
  }else if(type=="adaptive"){
    threshold=const*sqrt(log(p)/n)*sqrt(thetahat(X))
    ST=sign(Shat)*pmax(Shat-threshold,0)
  }
  return(ST)
}




# diagoal matrix

mdiag.r=function(n,d){
  M=matrix(0,n,n)
  m=length(d)
  for(i in 1:n){
    for(j in max(1,i-m+1):min(n,i+m-1)){
      M[i,j]=d[abs(i-j)+1]
    }
  }
  return(M)
}

maxd.r=function(x, d) {
  p=nrow(x)
  M=0;
  for(i in 1:(p-d)){
    for(j in min(i+d,p):p){
      M=max(abs(x[i,j]),M);
    }
  }
  return(M);
}
