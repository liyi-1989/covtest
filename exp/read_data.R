library(rhdf5)
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
X=X[c((NLON/2+1):NLON,1:(NLON/2)),,]


idlon=(1:18)*8
idlat=(1:18)*4
NLON=NLAT=18
LON=LON[idlon]
LAT=LAT[idlat]
X=X[idlon,idlat,]




