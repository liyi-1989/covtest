library(rcd)
library(igraph)
library("maps")
library("geosphere")
source('../clinet/utils_network.R')
library(scatterplot3d)
library(ggplot2)
library(glmnet)
library(e1071)
library(kernlab)
library(maptools)
library(fields)
data(wrld_simpl)
#https://www.esrl.noaa.gov/psd/data/gridded/data.cmap.html
#load("./dat/precip_mon_mean_mon_mean_removed_sub.RData")
load("./dat/precip_mon_mean_mon_mean_not_removed_sub.RData")
#load("../clinet/data/air_mon_mean_mon_mean_removed_sub.RData")
#load("../clinet/data/air_mon_mean_mon_mean_not_removed_sub.RData")
#load("../clinet/data/cor_rcd_matrix.RData")
source('utils_lscm.R')
library(mvtnorm)
library(tikzDevice)

idlon=1:72#1:144#1:72#(1:36)*2
idlat=7:32#16:63#7:32#(2:17)*2
NLON=length(idlon); NLAT=length(idlat)
LON=LON[idlon]
LAT=LAT[idlat]
X=X[idlon,idlat,]

dfv=NULL
count=0
for(i in 1:NLON){ # stack by row
  for(j in 1:NLAT){
    count=count+1
    dfv=rbind(dfv,c(count,i,j,LON[i],LAT[j],X[i,j,1]))#372+360+1
  }
}
colnames(dfv)=c("vertex","idlon","idlat","lon","lat","x")
dfv=data.frame(dfv)
p=dim(dfv)[1]
X1=NULL # Long Data (2D)
for(i in 1:NLON){
  for(j in 1:NLAT){
    X1=cbind(X1,X[i,j,])
  }
} # X1 is the final data matrix to work on. Vertex data frame is in dfv. Edge data frame need analysis with correlation
#plot_lonlat_df(dfv,vcol="x",region="world",CEX=2)
# sum(abs(X1)>1e10)
# which(abs(X)>1e3,arr.ind = T)
# plot(apply(X1,2,function(x){sum(abs(x)>1e3)==0}))
# X1=apply(X1,2,function(x){if(any(abs(x)>1e3)){x[abs(x)>1e3]=mean(x[abs(x)<1e3]);return(x)}else{return(x)}})
# X1[abs(X1)>1e3]=NA

#X1n=scale(X1)

t0=proc.time()
# Step 1: Estimate model parameters 
#S=cor(X1,use="pairwise.complete.obs")#+1e-5*diag(1:dim(X1)[2])
S=cor(X1)
hatfit=a_efratio_hat(S=S,p=nrow(S),d=2,la=1)
a_hat=hatfit[-length(hatfit)]
efratio_hat=hatfit[length(hatfit)]
#lstar=LS5m(a=a_hat,ef_ratio=efratio_hat,d=2,dd=9,sf=1,r=1)$ls5 # length(a)=1

# Step 2: Reconstruct the model with the estimated parameters (empirical Shat will be useless)
para1=list(n=dim(X1)[1],p=dim(X1)[2],a=a_hat,r=0.8,sf=1,se=efratio_hat,i0=NULL,j0=NULL,M=NULL,d=1,dd=5,py=NLAT,method="equalvar",N=100) 
para1$i0=round(para1$p*0.23); para1$j0=round(para1$p*0.77)
para2=para1; para2$d=2; para2$dd=9; 
covmodel=TCM(para2) # generating true covariance model
lstar=LS5(covmodel$S1,para2,Amat=NULL)$ls5
lstar
# Step 3: use model optimal lambda to find modified estimator (need empirical Shat and optimal lambda)
Sbar=sigmabar(S=S,mu=lstar,bandwidth=NLAT,dd=9,py=NLAT)

t1=proc.time()-t0

# par(mfrow=c(1,1))
# printM(S*(S>0.8))
# printM(Sbar*(Sbar>0.8))#Sbar[232,441] # 100 vs 207 225
# dfv[100,]
# dfv[207,]
# plot(X[6,10,],X[12,9,])
# 
# plot(X1[,232],X1[,441]) # Indian vs Brazil part
# plot(X1[,40],X1[,137]) # El nino part

# plot arc map
thres=5
S0=nn_ind(NLON,NLAT,thres,thres) # nearest neighbour indicator matrix
par(mfrow=c(1,1))
#net1=graph_from_adjacency_matrix((abs(S)>0.8)*(!S0),mode = "undirected")
net1=graph_from_adjacency_matrix(((abs(S)*(!S0))>0.83) & ((abs(S)*(!S0))<0.832),mode = "undirected") #0.525precip 0.515/150 # 0.8485 to avoid tele
dfe1=as_edgelist(net1)
plot_arc(dfv,dfe1,cap="Precipitation")

net1=graph_from_adjacency_matrix(((abs(S)*(!S0))>0.8485) ,mode = "undirected") #0.525precip 0.515/150 # 0.8485 to avoid tele
dfe1=as_edgelist(net1)
plot_arc(dfv,dfe1,cap="Precipitation")

#dfv[dfe1[,1],5]< -10 # 79 261 (14th row)

#net2=graph_from_adjacency_matrix((abs(Sbar)>0.8)*(!S0),mode = "undirected")
### paper  plot ******
net2=graph_from_adjacency_matrix((abs(Sbar)*(!S0))>1.1,mode = "undirected") #0.675/100precip 0.61305/150; 1.056
dfe2=as_edgelist(net2)
plot_arc(dfv,dfe2,cap="Precipitation")
nrow(dfe2)


net3=graph_from_adjacency_matrix((abs(S)*(!S0)-abs(Sbar)*(!S0))>0.08,mode = "undirected")
net3=graph_from_adjacency_matrix((S-Sbar)>0.2 & S>0.8 & !S0,mode = "undirected")
dfe3=as_edgelist(net3)
plot_arc(dfv,dfe3,cap="Precipitation")

#net3=graph_from_adjacency_matrix((abs(S)*(!S0)-abs(Sbar)*(!S0))>0.08,mode = "undirected")
### paper plot ******
net3=graph_from_adjacency_matrix((S-Sbar)>0.49 & S>0.4 & !S0,mode = "undirected")
dfe3=as_edgelist(net3)
plot_arc(dfv,dfe3,cap="Precipitation")

for(i in 1:nrow(dfe3)){
  cat(i,S[dfe3[i,1],dfe3[i,2]],Sbar[dfe3[i,1],dfe3[i,2]],"\n")
}

for(i in 1:nrow(dfe1)){
  cat(i,S[dfe1[i,1],dfe1[i,2]],Sbar[dfe1[i,1],dfe1[i,2]],"\n")
}

############### rank S as a long vector ###############
#### list of S Sbar S0
LS=matrix(0,sum(S0==0),3)
count=0
for(i in 1:p){
  for(j in i:p){
    if(!S0[i,j]){
      count=count+1
      LS[count,]=c(S[i,j],abs(S[i,j])-abs(Sbar[i,j]),Sbar[i,j])
    }
  }
}
#i00=688;j00=1398; i00=899; j00=1133; i00=1473; j00=1479
i00=660; j00=1077
i00=1078; j00=1398
q1=ecdf(LS[,1])
q2=ecdf(LS[,3])
q1(S[i00,j00])
q2(Sbar[i00,j00])

sum(LS[,1]>S[i00,j00])
sum(LS[,3]>Sbar[i00,j00])

tikz("./fig/real_plots.tex", width = 6.5, height = 3.25)
par(mfrow=c(1,2))
xx=X1[,232]; yy=X1[,441]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$5^\\circ$S, $35^\\circ$W",
     ylab="$5^\\circ$N, $95^\\circ$E",pch=19,main="Temperature")
abline(lm(xx~yy), col="blue",lty=2)
#text(c(-1,1),paste0("$\\rho=$",round(cor(xx,yy),2)),col="blue")
xx=X1[,40]; yy=X1[,137]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$5^\\circ$S, $155^\\circ$W",
     ylab="$5^\\circ$N, $95^\\circ$W",pch=19,main="Temperature")
abline(lm(xx~yy), col="blue",lty=2)
#text(c(-1,1),paste0("$\\rho=$",round(cor(xx,yy),2)),col="blue")
dev.off()



tikz("./fig/real_plots2.tex", width = 6.5, height = 6.5)

par(mfrow=c(2,2))
ii=79;jj=261
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$58.75^\\circ$S, $163.75^\\circ$W",
     ylab="$58.75^\\circ$S, $128.75^\\circ$W",pch=19,main="Precipitation")
abline(lm(xx~yy), col="blue",lty=2)

ii=79;jj=262
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$58.75^\\circ$S, $163.75^\\circ$W",
     ylab="$53.75^\\circ$S, $128.75^\\circ$W",pch=19,main="Precipitation")
abline(lm(xx~yy), col="blue",lty=2)


ii=79+NLAT; jj=261
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$58.75^\\circ$S, $158.75^\\circ$W",
     ylab="$58.75^\\circ$S, $128.75^\\circ$W",pch=19,main="Precipitation")
abline(lm(xx~yy), col="blue",lty=2)

ii=79+NLAT; jj=262
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$58.75^\\circ$S, $158.75^\\circ$W",
     ylab="$53.75^\\circ$S, $128.75^\\circ$W",pch=19,main="Precipitation")
abline(lm(xx~yy), col="blue",lty=2)

dev.off()


#### screening features ####
dfs=NULL
N=dim(S)[1]
for(i in 1:(N-2)){
  #print(i)
  for(j in (i+1):(N-1)){
    if((i%%NLON>=2)){
      if((abs(S[i,j])-mean(abs(S[i,j+1]),abs(S[i,j-1])))>0.3){
        dfs=rbind(dfs,c(i,j,S[i,j],S[i,j-1],S[i,j+1]))
      }
    }
  }
}

# Temp, removed, 5*5 deg, model length(a)=2, 2d model

i00=173; j00=7070
i00=57; j00=4559
i00=2544; j00=4028
i00=981; j00=1136
i00=871; j00=1398
i00=1078; j00=1398
i00=660; j00=1077

S[i00+c(-1,0,1),j00+c(-1,0,1)]
S[i00+c(-NLAT,0,NLAT),j00+c(-NLAT,0,NLAT)]
plot(X1[,i00],X1[,j00+NLAT])#,xlim=c(-10,10),ylim=c(-10,10)
plot(rank(X1[,i00]),rank(X1[,j00-1]))

############### paper plots ###############
###### real_plot3 ######
tikz("./fig/real_plots3.tex", width = 6.5, height = 6.5)
i00=1078; j00=1398; 
#i00=688; j00=1398
par(mfrow=c(2,2))
ii=i00;jj=j00
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$3.75^\\circ$S, $26.25^\\circ$W",
     ylab="$36.25^\\circ$N, $86.25^\\circ$E",pch=19,main="Precipitation")
#abline(lm(xx~yy), col="blue",lty=2)
text(2,8,paste0("cor=",round(cor(xx,yy),2)),col="blue")

ii=i00+1;jj=j00
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$1.25^\\circ$N, $26.25^\\circ$W",
     ylab="$36.25^\\circ$N, $86.25^\\circ$E",pch=19,main="Precipitation")
#abline(lm(xx~yy), col="blue",lty=2)
text(2,8,paste0("cor=",round(cor(xx,yy),2)),col="blue")

ii=i00;jj=j00+NLAT
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$3.75^\\circ$S, $26.25^\\circ$W",
     ylab="$36.25^\\circ$N, $91.25^\\circ$E",pch=19,main="Precipitation")
#abline(lm(xx~yy), col="blue",lty=2)
text(2,4,paste0("cor=",round(cor(xx,yy),2)),col="blue")

ii=i00;jj=j00-NLAT
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$3.75^\\circ$S, $26.25^\\circ$W",
     ylab="$36.25^\\circ$N, $81.25^\\circ$E",pch=19,main="Precipitation")
#abline(lm(xx~yy), col="blue",lty=2)
text(2,2,paste0("cor=",round(cor(xx,yy),2)),col="blue")

dev.off()


###### real_plot4 ######
tikz("./fig/real_plots4.tex", width = 6.5, height = 6.5)
i00=660; j00=1077; #i00=688; j00=1398
par(mfrow=c(2,2))
ii=i00;jj=j00
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$13.75^\\circ$S, $53.75^\\circ$W",
     ylab="$8.75^\\circ$S, $26.25^\\circ$E",pch=19,main="Precipitation")
#abline(lm(xx~yy), col="blue",lty=2)
text(2,8,paste0("cor=",round(cor(xx,yy),2)),col="blue")

ii=i00+1;jj=j00
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$8.75^\\circ$S, $53.75^\\circ$W",
     ylab="$8.75^\\circ$S, $26.25^\\circ$E",pch=19,main="Precipitation")
#abline(lm(xx~yy), col="blue",lty=2)
text(2,8,paste0("cor=",round(cor(xx,yy),2)),col="blue")

ii=i00;jj=j00+NLAT
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$13.75^\\circ$S, $53.75^\\circ$W",
     ylab="$8.75^\\circ$S, $31.25^\\circ$E",pch=19,main="Precipitation")
#abline(lm(xx~yy), col="blue",lty=2)
text(2,8,paste0("cor=",round(cor(xx,yy),2)),col="blue")

ii=i00;jj=j00-NLAT
xx=X1[,ii]; yy=X1[,jj]
plot(xx,yy,col=rgb(1,0.5,0,alpha=0.5),xlab="$13.75^\\circ$S, $53.75^\\circ$W",
     ylab="$8.75^\\circ$S, $21.25^\\circ$E",pch=19,main="Precipitation")
#abline(lm(xx~yy), col="blue",lty=2)
text(2,8,paste0("cor=",round(cor(xx,yy),2)),col="blue")

dev.off()


############ Add: find null cutoff ############

nk=1000
lmax=rep(0,nk)

for(k in 1:nk){
  #t0=proc.time()
  print(k)
  data = mvrnorm(dim(X)[3], mu = rep(0,p), Sigma = covmodel$S0)
  Shat=cor(data)
  lmax[k]=maxd.r(abs(Shat)*(!S0),6)
  rm(data)
  #proc.time()-t0
}

save(lmax,file="results/real_H0_lmax.RData")
load("results/real_H0_lmax.RData")

quantile(lmax,0.999) #0.4935794 
max(abs(LS[,1])) #0.889298
sum(abs(LS[,1])>quantile(lmax,0.999)) #55387

