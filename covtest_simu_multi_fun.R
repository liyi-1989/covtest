source("utils_lscm.R")
library(Rmpi)
library(snow)
library(quadprog)
library(Matrix)

paras=list(p=100,n=1000,a=0.9,r=0.8,sf=1,se=2,i0=NULL,j0=NULL,M=NULL,d=1,method="equalvar",N=100) # one-dimensional model

aas=c(0.7,0.8,0.9)#(0:5)/5
rrs=c(0.7,0.8,0.9)#(4:10)/10
ses=1:3
ns=c(500,1000,2000,3000)#c(2,4,6,8,10,20,30,40,50)*100
ps=c(100,200,300)#c(50,100,150,200,500,1000)

# aas=(0:5)/5
# rrs=(4:10)/10
# ses=1:3
# ns=c(2,4,6,8,10,20,30,40,50)*100
# ps=c(50,100,150,200,500,1000)

  
df=NULL
for(i1 in 1:length(aas)){
  for(i2 in 1:length(rrs)){
    for(i3 in 1:length(ses)){
      for(i4 in 1:length(ns)){
        for(i5 in 1:length(ps)){
            df=rbind(df,c(i1,i2,i3,i4,i5,aas[i1],rrs[i2],ses[i3],ns[i4],ps[i5]))
          
        }
      }
    }
  }
}
colnames(df)=c("id1","id2","id3","id4","i5","a","r","se","n","p")
save(df,file = "./results/results.df.RData")



simu=function(jobid=1,paras,df=NULL){
  write.table(NULL,paste0("./results/working_job_",jobid,".txt"))
  cat("1. Set Parameters ...\n")
  para1=paras
  para1$a=df[jobid,"a"]
  para1$r=df[jobid,"r"]
  para1$se=df[jobid,"se"]
  para1$n=df[jobid,"n"]
  para1$p=df[jobid,"p"]
  para1$i0=round(para1$p*c(0.15,0.23,0.45)); para1$j0=round(para1$p*c(0.6,0.77,0.85))
  
  covmodel=TCM(para1); S0=covmodel$S0; S1=covmodel$S1
  # make a copy
  n=para1$n; p=para1$p; i0=para1$i0; j0=para1$j0; N=para1$N 
  alpha=0.05
  filepath="./results/"
  filename=paste0("job_",jobid)
  
  cat("2. Simulation ...\n")
  R0=as.data.frame(matrix(NA,N,14))
  colnames(R0)=c("a_hat","efratio_hat","s_hat_F","s_bar_F","s_hat_2","s_bar_2","n_detect_hat", "n_detect_bar","max_hat","max_bar","max_hat_i","max_hat_j","max_bar_i","max_bar_j")
  R1=R0
  
  set.seed(1)
  for(k in 1:N){
    #============== H0 ==============
    # data = mvrnorm(n, mu = rep(0,p), Sigma = S0)
    # Shat=cov(data)
    # 
    # hatfit=a_efratio_hat(S=Shat,p=p,d=para1$d) # estimate a, se/sf
    # a_hat=hatfit[1]
    # efratio_hat=hatfit[2]
    # lstar=LS5m(a=a_hat,ef_ratio=ifelse(is.nan(efratio_hat),0,efratio_hat),d=para1$d)$ls5 # estimate lambda star
    # #lstar=ls1
    # Sbar=sigmabar(S=Shat,mu=lstar) # sigmabar
    # 
    # sij_hat=Shat[i0,j0]
    # sij_bar=Sbar[i0,j0]
    # s_hat_F=base::norm(Shat-S0,"F")
    # s_bar_F=base::norm(Sbar-S0,"F")
    # s_hat_2=base::norm(Shat-S0,"2")
    # s_bar_2=base::norm(Sbar-S0,"2")
    # max_hat=max(abs(Shat)*(S0==0))
    # max_bar=maxd.r(abs(Sbar[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(S0==0))
    # mm_hat=abs(Shat[2:(p-1),2:(p-1)])
    # mm_hat_ji=which(mm_hat == maxd.r(mm_hat,4), arr.ind = TRUE)
    # mm_bar=abs(Sbar[2:(p-1),2:(p-1)])
    # mm_bar_ji=which(mm_bar == maxd.r(mm_bar,4), arr.ind = TRUE)
    # 
    # #---------- top 3 positions ------------------
    # Shatoff=myoffDiag(Shat,3)
    # Sbaroff=myoffDiag(Sbar,3)
    # 
    # i0j0hat=which(abs(Shatoff)>tail(sort(abs(Shatoff)),2*4)[1],arr.ind = T)[1:3,]
    # i0j0bar=which(abs(Sbaroff)>tail(sort(abs(Sbaroff)),2*4)[1],arr.ind = T)[1:3,]
    # 
    # i0j0hatSet=apply(i0j0hat,1,function(x){paste0(as.character(x[1]),as.character(x[2]))})
    # i0j0barSet=apply(i0j0bar,1,function(x){paste0(as.character(x[1]),as.character(x[2]))})
    # i0j0Set=paste0(j0,i0)
    # n_detect_hat=sum(is.element(i0j0hatSet,i0j0Set))
    # n_detect_bar=sum(is.element(i0j0barSet,i0j0Set))
    # #-----------------------------------
    # 
    # 
    # R0[k,]=c(a_hat,efratio_hat,s_hat_F,s_bar_F,s_hat_2,s_bar_2, n_detect_hat, n_detect_bar,max_hat,max_bar,mm_hat_ji[1,],mm_bar_ji[1,])
    #============== H1 ==============
    data = mvrnorm(n, mu = rep(0,p), Sigma = S1)
    Shat=cov(data)
    
    hatfit=a_efratio_hat(S=Shat,p=p,d=para1$d) # estimate a, se/sf
    a_hat=hatfit[1]
    efratio_hat=hatfit[2]
    lstar=LS5m(a=a_hat,ef_ratio=ifelse(is.nan(efratio_hat),0,efratio_hat),d=para1$d)$ls5 # estimate lambda star
    #lstar=ls1
    Sbar=sigmabar(S=Shat,mu=lstar) # sigmabar
    
    sij_hat=Shat[i0,j0]
    sij_bar=Sbar[i0,j0]
    s_hat_F=base::norm(Shat-S1,"F")
    s_bar_F=base::norm(Sbar-S1,"F")
    s_hat_2=base::norm(Shat-S1,"2")
    s_bar_2=base::norm(Sbar-S1,"2")
    max_hat=max(abs(Shat)*(S0==0))
    max_bar=maxd.r(abs(Sbar[2:(p-1),2:(p-1)]),4)#max(abs(Sbar)*(S0==0))
    mm_hat=abs(Shat[2:(p-1),2:(p-1)])
    mm_hat_ji=which(mm_hat == maxd.r(mm_hat,4), arr.ind = TRUE)
    mm_bar=abs(Sbar[2:(p-1),2:(p-1)])
    mm_bar_ji=which(mm_bar == maxd.r(mm_bar,4), arr.ind = TRUE)
    
    #---------- top 3 positions ------------------
    Shatoff=myoffDiag(Shat,3)
    Sbaroff=myoffDiag(Sbar,3)
    i0j0Set=paste0(j0,i0) # true set 
    i0j0=data.frame(j0=j0,i0=i0)
    
    i0j0hat=which(abs(Shatoff)>tail(sort(abs(Shatoff)),2*4)[1],arr.ind = T)[1:3,]
    i0j0hatSet=apply(i0j0hat,1,function(x){paste0(as.character(x[1]),as.character(x[2]))})
    n_detect_hat=sum(is.element(i0j0hatSet,i0j0Set))
    
    i0j0bar=which(abs(Sbaroff)>tail(sort(abs(Sbaroff)),2*31)[1],arr.ind = T)[1:30,]
    i0j0bar_value=rep(0,30)
    for(ii in 1:30){i0j0bar_value[ii]=Sbaroff[i0j0bar[ii,1],i0j0bar[ii,2]]}
    i0j0bar=i0j0bar[order(abs(i0j0bar_value),decreasing = T),]
    i0j0barSet=apply(i0j0bar,1,function(x){paste0(as.character(x[1]),as.character(x[2]))})
    
    #sqrt(sum((i0j0[1,]-i0j0[1,])^2))
    n4=rep(0,4)
    ii=0
    while(sum(n4)<3){
      ii=ii+1
      if(is.element(i0j0barSet[ii],i0j0Set[1])){n4[1]=1;next}
      if(is.element(i0j0barSet[ii],i0j0Set[2])){n4[2]=1;next}
      if(is.element(i0j0barSet[ii],i0j0Set[3])){n4[3]=1;next}
      
      d1=sqrt(sum((i0j0bar[ii,]-i0j0[1,])^2))
      d2=sqrt(sum((i0j0bar[ii,]-i0j0[2,])^2))
      d3=sqrt(sum((i0j0bar[ii,]-i0j0[3,])^2))
      
      if(d1>1.42 & d2>1.42 & d3>1.42){n4[4]=n4[4]+1}
      
    }
    
    n_detect_bar=sum(n4[1:3])#sum(is.element(i0j0barSet,i0j0Set))
    #print(n_detect_bar)
    #print(ii)
    #-----------------------------------
    R1[k,]=c(a_hat,efratio_hat,s_hat_F,s_bar_F,s_hat_2,s_bar_2, n_detect_hat, n_detect_bar,max_hat,max_bar,mm_hat_ji[1,],mm_bar_ji[1,])
  }
  
  # sij_0=S0[i0,j0]
  # sij_1=S1[i0,j0]
  # lmax10=R0[,"max_hat"]
  # lmax20=R0[,"max_bar"]
  # lmax11=R1[,"max_hat"]
  # lmax21=R1[,"max_bar"]
  # ############ cut-off ############
  # TS10=n*lmax10^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
  # CUT.simu10=quantile(TS10,1-alpha) # CUT.asym10=-log(8*pi)-2*log(log(1/(1-alpha)))#y_alpha
  # TS20=n*lmax20^2-4*log(p)+log(log(p)) # test statistics of lmax(coherence)
  # CUT.simu20=quantile(TS20,1-alpha) # type I error: mean(TS10>CUT.simu10) # by quantile definition, it is always 0.05
  # ############ type II error: power ############
  # TS11=n*lmax11^2-4*log(p)+log(log(p)) # coherence of S
  # EP11=mean(TS11>CUT.simu10) # Empirical power
  # TS21=n*lmax21^2-4*log(p)+log(log(p)) # coherence of Sbar
  # EP21=mean(TS21>CUT.simu20) # Empirical power

  cat("3. Save results ...\n")
  #save(R0,R1,sij_0,sij_1,EP11,EP21,TS10,TS11,TS20,TS21,para1,file=paste0(filepath,filename,".RData"))
  save(R1,para1,file=paste0(filepath,filename,".RData"))
  
  write.table(NULL,paste0("./results/finishing_job_",jobid,".txt"))
  return(NULL)
}

