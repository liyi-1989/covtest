source('utils_lscm.R')
library(mvtnorm)
library(tikzDevice)
library(ggplot2)

para1=list(p=100,n=1000,a=0.9,r=0,sf=1,se=2,i0=NULL,j0=NULL,M=NULL,d=1,dd=5,py=NULL,method="equalvar",N=100) # one-dimensional model
para1$i0=round(para1$p*0.23); para1$j0=round(para1$p*0.77)
para1$M=Mp(para1$p)
#para2=para1; para2$d=2; para2$dd=9; para2$py=10 # two-dimensional model 
ll=matrix(c(1,0,0,0,0),5,1)

# start simulation
df=NULL
for(a in c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0)){
  for(se in 1:3){
    print(a)
    para1$a=a; para1$se=se
    covmodel=TCM(para1) # generating true covariance model
    S0=covmodel$S0; S1=covmodel$S1
    ls5fit=LS5(S1,para1)
    
    var1=t(ls5fit$ls5)%*%ls5fit$Dmat%*%ls5fit$ls5
    var0=t(ll)%*%ls5fit$Dmat%*%ll
    
    df=rbind(df,c(a,se,var0,var1,sqrt(var1/var0)))
  }
}

df=as.data.frame(df)
colnames(df)=c("a","se","var0","var1","ratio")

ggplot(df, aes(x=a, y=ratio^2, group=factor(se), colour=factor(se) ) ) + geom_line(size=1)+ theme(legend.title = element_blank())+
  ggtitle(" ")+theme(plot.title = element_text(hjust = 0.5))+xlab("a") + ylab("Reduction Ratio")

ggsave("fig/rr.pdf")

#plot(df$a,df$ratio,col=df$se+1,lty=1)

