library(tikzDevice)
#setwd("./results/small_simu_each_3/")
setwd("./results/large_simu/")
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



dfep=cbind(df,EP1,EP2)[,-(1:5)]
dfep=rbind(cbind(df,EP1,1),cbind(df,EP2,2))
colnames(dfep)=c(colnames(df),"EP","hat1bar2")
dfep=as.data.frame(dfep)

library(reshape2)


dfepw=reshape2::dcast(dfep, se+a+r~p+n+hat1bar2, value.var = "EP")
write.table(dfepw,file = "dfepw.csv",row.names = F)

# Make the power curve table
dfepw=reshape2::dcast(subset(dfep,se %in% c(1,2) & p %in% c(50,100,500) & a %in% c(0.4,0.6,0.8) & r %in% c(0.6,0.7,0.8) & n %in% c(400,1000,2000)), 
                      se+a+r~p+n+hat1bar2, value.var = "EP")
write.table(dfepw,file = "temp.csv",row.names = F)


tikz("../../fig/power_curve_samplesize.tex", width = 6, height = 9)
par(mfrow=c(3,2))
#----------------------------------------------------------------------
# plot power curve with sample size n
dfep1=subset(dfep,p==100 & r==0.6 & se==2 &n>400)
dfep1t1=subset(dfep1,a==0.8&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.8&hat1bar2==2)
plot(dfep1t1$n,dfep1t1$EP,col="orange",type="b",lwd=2,xlab = "sample size", ylab="power", main="power curve")
lines(dfep1t2$n,dfep1t2$EP,col="blue",type="b",lwd=2)
dfep1t1=subset(dfep1,a==0.6&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.6&hat1bar2==2)
lines(dfep1t1$n,dfep1t1$EP,col="orange",type="b",lwd=1)
lines(dfep1t2$n,dfep1t2$EP,col="blue",type="b",lwd=1)
legend("bottom",c("$\\hat{\\Sigma},a=0.8$","$\\bar{\\Sigma},a=0.8$","$\\hat{\\Sigma},a=0.6$","$\\bar{\\Sigma},a=0.6$"),
       col = c("orange","blue","orange","blue"),lwd=c(2,2,1,1),lty=c(1,1,1,1),pch=c(1,1,1,1),cex=0.75,bty = "n")
#dev.off()

# plot power curve with signal strength
#tikz("../../fig/power_curve_signal.tex", width = 3.25, height = 3.25)
dfep1=subset(dfep,p==100 & se==2 &n==1000)
dfep1t1=subset(dfep1,a==0.8&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.8&hat1bar2==2)
plot(dfep1t1$r,dfep1t1$EP,col="orange",type="b",lwd=2,xlab = "signal strength", ylab="power", main="power curve",ylim = c(0,1))
lines(dfep1t2$r,dfep1t2$EP,col="blue",type="b",lwd=2)
dfep1t1=subset(dfep1,a==0.6&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.6&hat1bar2==2)
lines(dfep1t1$r,dfep1t1$EP,col="orange",type="b",lwd=1)
lines(dfep1t2$r,dfep1t2$EP,col="blue",type="b",lwd=1)
legend("topleft",c("$\\hat{\\Sigma},a=0.8$","$\\bar{\\Sigma},a=0.8$","$\\hat{\\Sigma},a=0.6$","$\\bar{\\Sigma},a=0.6$"),
       col = c("orange","blue","orange","blue"),lwd=c(2,2,1,1),lty=c(1,1,1,1),pch=c(1,1,1,1),cex=0.75,bty = "n")
#----------------------------------------------------------------------
# plot power curve with sample size n
dfep1=subset(dfep,p==200 & r==0.6 & se==2 &n>400)
dfep1t1=subset(dfep1,a==0.8&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.8&hat1bar2==2)
plot(dfep1t1$n,dfep1t1$EP,col="orange",type="b",lwd=2,xlab = "sample size", ylab="power", main="power curve")
lines(dfep1t2$n,dfep1t2$EP,col="blue",type="b",lwd=2)
dfep1t1=subset(dfep1,a==0.6&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.6&hat1bar2==2)
lines(dfep1t1$n,dfep1t1$EP,col="orange",type="b",lwd=1)
lines(dfep1t2$n,dfep1t2$EP,col="blue",type="b",lwd=1)
legend("bottom",c("$\\hat{\\Sigma},a=0.8$","$\\bar{\\Sigma},a=0.8$","$\\hat{\\Sigma},a=0.6$","$\\bar{\\Sigma},a=0.6$"),
       col = c("orange","blue","orange","blue"),lwd=c(2,2,1,1),lty=c(1,1,1,1),pch=c(1,1,1,1),cex=0.75,bty = "n")
#dev.off()

# plot power curve with signal strength
#tikz("../../fig/power_curve_signal.tex", width = 3.25, height = 3.25)
dfep1=subset(dfep,p==200 & se==2 &n==1000)
dfep1t1=subset(dfep1,a==0.8&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.8&hat1bar2==2)
plot(dfep1t1$r,dfep1t1$EP,col="orange",type="b",lwd=2,xlab = "signal strength", ylab="power", main="power curve",ylim = c(0,1))
lines(dfep1t2$r,dfep1t2$EP,col="blue",type="b",lwd=2)
dfep1t1=subset(dfep1,a==0.6&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.6&hat1bar2==2)
lines(dfep1t1$r,dfep1t1$EP,col="orange",type="b",lwd=1)
lines(dfep1t2$r,dfep1t2$EP,col="blue",type="b",lwd=1)
legend("topleft",c("$\\hat{\\Sigma},a=0.8$","$\\bar{\\Sigma},a=0.8$","$\\hat{\\Sigma},a=0.6$","$\\bar{\\Sigma},a=0.6$"),
       col = c("orange","blue","orange","blue"),lwd=c(2,2,1,1),lty=c(1,1,1,1),pch=c(1,1,1,1),cex=0.75,bty = "n")
#----------------------------------------------------------------------
# plot power curve with sample size n
dfep1=subset(dfep,p==500 & r==0.6 & se==2 &n>400)
dfep1t1=subset(dfep1,a==0.8&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.8&hat1bar2==2)
plot(dfep1t1$n,dfep1t1$EP,col="orange",type="b",lwd=2,xlab = "sample size", ylab="power", main="power curve")
lines(dfep1t2$n,dfep1t2$EP,col="blue",type="b",lwd=2)
dfep1t1=subset(dfep1,a==0.6&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.6&hat1bar2==2)
lines(dfep1t1$n,dfep1t1$EP,col="orange",type="b",lwd=1)
lines(dfep1t2$n,dfep1t2$EP,col="blue",type="b",lwd=1)
legend("bottom",c("$\\hat{\\Sigma},a=0.8$","$\\bar{\\Sigma},a=0.8$","$\\hat{\\Sigma},a=0.6$","$\\bar{\\Sigma},a=0.6$"),
       col = c("orange","blue","orange","blue"),lwd=c(2,2,1,1),lty=c(1,1,1,1),pch=c(1,1,1,1),cex=0.75,bty = "n")
#dev.off()

# plot power curve with signal strength
#tikz("../../fig/power_curve_signal.tex", width = 3.25, height = 3.25)
dfep1=subset(dfep,p==500 & se==2 &n==1000)
dfep1t1=subset(dfep1,a==0.8&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.8&hat1bar2==2)
plot(dfep1t1$r,dfep1t1$EP,col="orange",type="b",lwd=2,xlab = "signal strength", ylab="power", main="power curve",ylim = c(0,1))
lines(dfep1t2$r,dfep1t2$EP,col="blue",type="b",lwd=2)
dfep1t1=subset(dfep1,a==0.6&hat1bar2==1)
dfep1t2=subset(dfep1,a==0.6&hat1bar2==2)
lines(dfep1t1$r,dfep1t1$EP,col="orange",type="b",lwd=1)
lines(dfep1t2$r,dfep1t2$EP,col="blue",type="b",lwd=1)
legend("topleft",c("$\\hat{\\Sigma},a=0.8$","$\\bar{\\Sigma},a=0.8$","$\\hat{\\Sigma},a=0.6$","$\\bar{\\Sigma},a=0.6$"),
       col = c("orange","blue","orange","blue"),lwd=c(2,2,1,1),lty=c(1,1,1,1),pch=c(1,1,1,1),cex=0.75,bty = "n")

#----------------------------------------------------------------------
dev.off()






