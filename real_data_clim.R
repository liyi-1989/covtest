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
load("../clinet/data/air_mon_mean_mon_mean_removed_sub.RData")
load("../clinet/data/air_mon_mean_mon_mean_not_removed_sub.RData")
load("../clinet/data/cor_rcd_matrix.RData")
source('utils_lscm.R')
library(mvtnorm)


idlon=(1:18)*4
idlat=(1:18)*2
NLON=NLAT=18
LON=LON[idlon]
LAT=LAT[idlat]
X=X[idlon,idlat,]

dfv=NULL
count=0
for(i in 1:NLON){
  for(j in 1:NLAT){
    count=count+1
    dfv=rbind(dfv,c(count,i,j,LON[i],LAT[j],X[i,j,372+360+1]))
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
plot_lonlat_df(dfv,vcol="x",region="world",CEX=2)


#X1n=scale(X1)
S=cor(X1)
hatfit=a_efratio_hat(S=S,p=nrow(S),d=2)
a_hat=hatfit[1]
efratio_hat=hatfit[2]
lstar=LS5m(a=a_hat,ef_ratio=efratio_hat,d=2)$ls5 
Sbar=sigmabar(S=S,mu=lstar) 

par(mfrow=c(1,1))
printM(S*(S>0.8))
printM(Sbar*(Sbar>0.75))# 100 vs 207 225
dfv[100,]
dfv[207,]
plot(X[6,10,],X[12,9,])
