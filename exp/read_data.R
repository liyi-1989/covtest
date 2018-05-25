library(rhdf5)
#https://www.esrl.noaa.gov/psd/data/gridded/data.cmap.html
fileloc="C:/Users/liyi/Desktop/xw/precip.mon.mean.nc"
fileloc="./dat/precip.mon.mean.nc"
h5ls(fileloc)

LON=h5read(fileloc,"lon")
LAT=h5read(fileloc,"lat")
TIME=h5read(fileloc,"time")
X=h5read(fileloc,"precip") 
NLON=length(LON)
NLAT=length(LAT)

LAT=rev(LAT)
LON=c(LON[(NLON/2+1):NLON]-360,LON[1:(NLON/2)])
X=X[, length(LAT):1, ]
X=X[c((NLON/2+1):NLON,1:(NLON/2)),,-dim(X)[3]]

# remove monthly mean
X_mon_mean=X[,,1:12]
for(i in 0:11){
  X0=X[,,(1:dim(X)[3])%%12==i]
  X_mon_mean[,,i+1]=apply(X0,c(1,2),mean) # [12,1,2,3,...,11]
}

for(i in 1:dim(X)[3]){
  X[,,i]=X[,,i]-X_mon_mean[,,(i%%12)+1]
}

############### 1.1 Get subsample ################
idsellon=2*(1:(NLON/2))-1
idsellat=2*(1:(NLAT/2))-1
LON=LON[idsellon]
LAT=LAT[idsellat]
NLON=length(LON)
NLAT=length(LAT)
X=X[idsellon,idsellat,]
save(LON,LAT,NLON,NLAT,X,file="./dat/precip_mon_mean_mon_mean_removed_sub.RData")
save(LON,LAT,NLON,NLAT,X,file="./dat/precip_mon_mean_mon_mean_not_removed_sub.RData")








idlon=(1:18)*8
idlat=(1:18)*4
NLON=NLAT=18
LON=LON[idlon]
LAT=LAT[idlat]
X=X[idlon,idlat,]




