#library(tikzDevice)
#library(reshape2)
library(tidyverse)
library(ggplot2)
#library(gridExtra)

setwd("./results/simu_multi_960/")
setwd("./results/simu_multi_300/")
load("./results.df.RData")

N=nrow(df)
mean0=as.data.frame(matrix(NA,N,8))
colnames(mean0)=c("n_detect_hat","n_detect_bar","n","p","a","r","se","sf")
sd1=sd0=mean1=mean0

for(jobid in 1:N){
  load(paste0("job_",jobid,".RData"))
  mean0[jobid,]=c(apply(R0[,c("n_detect_hat","n_detect_bar")],2,mean),para1$n,para1$p,para1$a,para1$r,para1$se,para1$sf)
  sd0[jobid,]=c(apply(R0[,c("n_detect_hat","n_detect_bar")],2,sd),para1$n,para1$p,para1$a,para1$r,para1$se,para1$sf)
  mean1[jobid,]=c(apply(R1[,c("n_detect_hat","n_detect_bar")],2,mean),para1$n,para1$p,para1$a,para1$r,para1$se,para1$sf)
  sd1[jobid,]=c(apply(R1[,c("n_detect_hat","n_detect_bar")],2,sd),para1$n,para1$p,para1$a,para1$r,para1$se,para1$sf)
}

save(df,mean0,mean1,sd0,sd1,file = "multi_300.RData")

load("multi_300.RData")

#-------------------- n ---------------------------
means=subset(mean1,a==0.9 & r==0.9 & se==2 & p==300 & n>500)
df1=data.frame(value=means$n_detect_hat,var=means$n)
df1$estimator="1"
df2=data.frame(value=means$n_detect_bar,var=means$n)
df2$estimator="2"
df=rbind(df1,df2)
ggplot(df,aes(x=as.factor(var),y=value,fill=estimator))+
  geom_bar(stat="identity", position=position_dodge())+
  xlab("n")+ylab("correctly identified pairs")+#ggtitle("Averaged recoveries")+
  theme(legend.position = c(0.85, 0.75))+
  ggtitle(expression("p=300,a="~rho~"=0.9")) +
  theme(plot.title = element_text(hjust = 0.5))+
  ggsave("../../fig/multi_n.pdf")
#-------------------- p ---------------------------
means=subset(mean1,a==0.9 & r==0.9 & se==2 & p>0 & n==1000)
df1=data.frame(value=means$n_detect_hat,var=means$p)
df1$estimator="1"
df2=data.frame(value=means$n_detect_bar,var=means$p)
df2$estimator="2"
df=rbind(df1,df2)
ggplot(df,aes(x=as.factor(var),y=value,fill=estimator))+
  geom_bar(stat="identity", position=position_dodge())+
  xlab("p")+ylab("correctly identified pairs")+#ggtitle("Averaged recoveries")+
  theme(legend.position = c(0.85, 0.75))+
  ggtitle(expression("n=1000,a="~rho~"=0.9")) +
  theme(plot.title = element_text(hjust = 0.5))+
  ggsave("../../fig/multi_p.pdf")
 

#-------------------- table ---------------------------
library(data.table)

meant=subset(mean1, se==1 & p %in% c(100,200,300) & n %in% c(1000,2000,3000) & a %in% c(0.7,0.8,0.9) & r %in% c(0.7,0.8,0.9),
             select = c("n","p","a","r","se","n_detect_hat","n_detect_bar"))
meant1=data.table::dcast(setDT(meant), n+a ~ p+r,value.var = c("n_detect_hat","n_detect_bar"))
meant1=meant1[,c(1,2,rep(c(1,10),9)+rep(0:8,each=2)+2)]
write.csv(meant1,file="multi_1.csv")

meant=subset(mean1, se==2 & p %in% c(100,200,300) & n %in% c(1000,2000,3000) & a %in% c(0.7,0.8,0.9) & r %in% c(0.7,0.8,0.9),
             select = c("n","p","a","r","se","n_detect_hat","n_detect_bar"))
meant1=data.table::dcast(setDT(meant), n+a ~ p+r,value.var = c("n_detect_hat","n_detect_bar"))
meant1=meant1[,c(1,2,rep(c(1,10),9)+rep(0:8,each=2)+2)]
write.csv(meant1,file="multi_2.csv")

meant=subset(mean1, se==3 & p %in% c(100,200,300) & n %in% c(1000,2000,3000) & a %in% c(0.7,0.8,0.9) & r %in% c(0.7,0.8,0.9),
             select = c("n","p","a","r","se","n_detect_hat","n_detect_bar"))
meant1=data.table::dcast(setDT(meant), n+a ~ p+r,value.var = c("n_detect_hat","n_detect_bar"))
meant1=meant1[,c(1,2,rep(c(1,10),9)+rep(0:8,each=2)+2)]
write.csv(meant1,file="multi_3.csv")




