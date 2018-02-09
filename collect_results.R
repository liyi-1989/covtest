setwd("./results/small_simu_each_3/")
load("./results.df.RData")


mean0=as.data.frame(matrix(NA,nrow(df),18))
colnames(mean0)=c("a_hat","efratio_hat","sij_hat","sij_bar","s_hat_F","s_bar_F","s_hat_2","s_bar_2",
                 "max_hat","max_bar","max_hat_i","max_hat_j","max_bar_i","max_bar_j","s_ij","a","se","sf")
mean1=mean0
sd0=as.data.frame(matrix(NA,nrow(df),14))
colnames(sd0)=c("a_hat","efratio_hat","sij_hat","sij_bar","s_hat_F","s_bar_F","s_hat_2","s_bar_2",
                  "max_hat","max_bar","max_hat_i","max_hat_j","max_bar_i","max_bar_j")
sd1=sd0
EP1=EP2=rep(NA,nrow(df))

for(jobid in 1:nrow(df)){
  load(paste0("job_",jobid,".RData"))
  mean0[jobid,]=c(apply(R0,2,mean),sij_0,para1$a,para1$se,para1$sf)
  sd0[jobid,]=apply(R0,2,sd)
  mean1[jobid,]=c(apply(R1,2,mean),sij_1,para1$a,para1$se,para1$sf)
  sd1[jobid,]=apply(R1,2,sd)
  EP1[jobid]=EP11
  EP2[jobid]=EP21
}

#XXX=cbind(df,mean0,sd0)

dfep=cbind(df,EP1,EP2,EP2>=EP1)[,-(1:5)]


plot(dfep[,"n"],dfep[,"EP1"])

