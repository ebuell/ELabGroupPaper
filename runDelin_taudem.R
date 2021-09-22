#Run R delineation for !!!!!!!!!!!one specific DEM and DEM process!!!!!!!!!!!!!!
#Save the outputs from the delineation (slope, TIC, TIV, and SCA)
#to google drive
#code sourced in part from:
# https://hydrology.usu.edu/taudem/taudem5/TauDEMRScript.txt
# must install tuadem and mpich (https://www.mpich.org/downloads/) first

#rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,GSODR,curl,elevatr,raster,soilDB,rgdal,googlesheets4,shapefiles,sp,rgeos)
pacman::p_load(reticulate,sf,beepr)
library(beepr);library('RColorBrewer');library(corrplot)
options("download.file.extra"="-L -k")
options("download.file.method"="curl")

# Elyce's directory
basedir="G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Spatial Data/DEMs"
diroptions = c("1as","1_3as","1m","1m")
fileoptions = c("/1asproj.tif","/1_3asproj_trim.tif","/1mLOC_proj&trim.tif","/1mDC_proj&trim.tif")
#current clipped boundaries for 1/3as: xlim = -73,-73.41; ylim = 43.9,44.35
#current clipped boundaries for 1mDC: xlim = -73.4, -73.26; ylim = 43.92,44.09
#current clipped boundaries for 1mLOC: xlim = -73.26, -73.05; ylim = 44.06,44.28
#for some reason this does not work for 1/3as Dinf DC

#decide which DEM
demorigin = readline('1 for usgs, 2 for LIDAR ')
if(demorigin==1){
  aschoice = readline('1 for 1arcsec, 2 for 1/3arcsec, 3 for 1m ')
  if(aschoice==1){DEMchoice = 1;name = "1arcsec"}
  if(aschoice==2){DEMchoice = 2;name = "1_3arcsec"}
  if(aschoice==3){DEMchoice = 3;name = "1meter"}
}
if(demorigin==2){
  aschoice = readline('1 for 2010, 2 for 2018 ')
  if(aschoice==1){DEMchoice = 3;name = "2010LIDAR"}
  if(aschoice==2){DEMchoice = 4;name = "2018LIDAR"}
}

#decide which watershed
watershednum = readline('1 for Little Otter Creek and 2 for Dead Creek')
if(watershednum==1){
  watershedstring = "/Shapes/"
  watershedname = "LOC"
}else{
  watershedstring = "/Shapes-Deadcr/"
  watershedname = "DC"  
  if(aschoice==3){DEMchoice = 4}
}


d8ordinf = readline('1 for D8, 2 for Dinf ')
if(d8ordinf==1){d8ordinfstr = "d8"}
if(d8ordinf==2){d8ordinfstr = "dinf"}

#Get and process DEM desired
folder = paste(diroptions[DEMchoice],sep = "")
url=paste0(basedir,"/",folder,fileoptions[DEMchoice])
dem_utm=raster(url)
#Projection change
proj4_utm = paste0("+proj=utm +zone=", 18, " +datum=WGS84 +units=m +no_defs")
proj4_ll = "+proj=longlat"; crs_ll=CRS(proj4_ll); crs_utm=CRS(proj4_utm)
#dem_utm <- projectRaster(dem_utm,crs = crs_ll)

#get and plot watershed shape
url=paste0("G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Spatial Data/Assorted Shape Files",watershedstring,"subs1.shp")
roughboundary = read.shp(url)
roughboundary = SpatialPoints(cbind(as.numeric(unlist(roughboundary$shp[[1]]$points[,1])),as.numeric(unlist(roughboundary$shp[[1]]$points[,2]))),proj4string=crs_utm)
roughboundary = spTransform(roughboundary,crs_ll)



plot(dem_utm)
points(roughboundary)

###End DEM manipulations

#Import outlet point
url=paste0("G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Spatial Data/Assorted Shape Files/",watershedstring,"outlets1.shp")
outlet = readOGR(url)
outlet = spTransform(outlet,crs_ll) #check that we are in the correct coordinate system
#chose 1 outlet
ind = which.min(extract(dem_utm,outlet))
outlet = SpatialPoints(matrix(c(outlet$Long_[ind],outlet$Lat[ind]),1,2),proj4string = crs_ll)
points(outlet,col = 'blue',pch = 16)

#pit remove
setwd("C:/Users/ebuell/Documents/RDump_TauDEM")
writeRaster(dem_utm,"logan",format = "GTiff",overwrite = TRUE)
system("mpiexec -n 8 pitremove -z logan.tif -fel loganfel.tif")
dem_utm_pitrm=raster("loganfel.tif")
plot(dem_utm_pitrm)


#slope and contributing area calculation
if(d8ordinf=="1"){
  #D8 flow directions
  system("mpiexec -n 8 D8Flowdir -p loganp.tif -sd8 logans.tif -fel loganfel.tif",show.output.on.console=F,invisible=F)
  beep('coin')
  p=raster("loganp.tif")
  plot(p)
  slp=raster("logans.tif")
  plot(slp)
  # Contributing area
  system("mpiexec -n 8 AreaD8 -p loganp.tif -ad8 logana.tif")
  ca=raster("logana.tif")
  #plot(ca)
  sca = ca*res(dem_utm)[1]
  plot(log(sca))
}
if(d8ordinf=="2"){
  # DInf flow directions
  system("mpiexec -n 8 DinfFlowdir -ang loganang.tif -slp loganslp.tif -fel loganfel.tif",show.output.on.console=F,invisible=F)
  beep('coin')
  ang=raster("loganang.tif")
  #plot(ang)
  slp=raster("loganslp.tif")
  #plot(slp)
  # Contributing area
  system("mpiexec -n 8 AreaDinf -ang loganang.tif -sca logansca.tif")
  sca=raster("logansca.tif")
  plot(log(sca))
}
beep('coin')


#define threshhold for picking outlet point
threshhold = ceiling(max(values(sca),na.rm = TRUE)/40)#/res(dem_utm)[1])
if(d8ordinf=="1"){system(paste0("mpiexec -n 8 Threshold -ssa loganp.tif -src logansrc.tif -thresh ",threshhold))
}else{system(paste0("mpiexec -n 8 Threshold -ssa loganang.tif -src logansrc.tif -thresh ",threshhold))}
src=raster("logansrc.tif")
plot(src)
points(outlet)

#write outlet
shapefile(outlet, "approxoutlet.shp",overwrite = TRUE)

# Move Outlet
if(d8ordinf=="1"){system("mpiexec -n 8 moveoutletstostreams -p loganp.tif -src logansrc.tif -o approxoutlet.shp -om Outlet.shp")
}else{system("mpiexec -n 8 moveoutletstostreams -p loganang.tif -src logansrc.tif -o approxoutlet.shp -om Outlet.shp")}
outpt=read.shp("Outlet.shp")

# Contributing area upstream of outlet
if(d8ordinf=="1"){system("mpiexec -n 8 Aread8 -p loganp.tif -o Outlet.shp -ad8 loganssa.tif")
}else{system("mpiexec -n 8 Aread8 -p loganang.tif -o Outlet.shp -ad8 loganssa.tif")}

ssa=raster("loganssa.tif")
plot(ssa) 
points(roughboundary,cex=.1)
if(length(values(ssa)[!is.na(values(ssa))])<100){beep('mario');stop()}

#clip the catchment area and slp to the watershed
sca = mask(sca,ssa)
slp = mask(slp,ssa)
# riverrast = mask(src,ssa); values(riverrast)[which(values(riverrast)==0)] = NA
# if(DEMchoice==3|DEMchoice==4){
#   if(DEMchoice==3){thinby = 3}
#   if(DEMchoice==4){thinby = 5}
#   ind = which(values(riverrast)==1)
#   for(i in 1:length(ind)){
#     if(ind[i]%%thinby==0){ind[i] = NA}
#   }
#   ind = ind[!is.na(ind)]
#   values(riverrast)[ind] = NA
# }
#river = SpatialPointsDataFrame(rasterToPoints(riverrast),data.frame(hold = matrix(1,length(rasterToPoints(riverrast)[,1]),1)))
#plot(river,pch = 15)
if(d8ordinf==2){p = ang}
p = mask(p,ssa)


#calculated tiv
TI = log((sca+1)/(slp+0.00001))
plot(TI)

#split into TICs
pacman::p_load(classInt)
nTIclass=10 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI); v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
#plot(TIC)#,col = rainbow(10))

#trim all rasters for export
p = trim(p,values=NA)
plot(p)
beep('coin')
sca = trim(sca,values=NA)
plot(sca)
beep('coin')
slp = trim(slp,values=NA)
plot(slp)
beep('coin')
TI = trim(TI,values=NA)
plot(TI)
beep('coin')
TIC = trim(TIC,values=NA)
plot(TIC)
beep('coin')


#stop()
#### UNITS of things!!! #####
#Slope: length/length
#SCA: meters
#TIC: unitless
#TI: something complicated that relates to slope and SCA units
#beep("mario")
#stop()
writeRaster(TIC,paste0("G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Spatial Data/Spatially Derived/TICs/TIC_",name,d8ordinfstr,"_",watershedname),format = "GTiff",overwrite = TRUE)
writeRaster(TI,paste0("G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Spatial Data/Spatially Derived/TIVs/TIV_",name,d8ordinfstr,"_",watershedname),format = "GTiff",overwrite = TRUE)
writeRaster(slp,paste0("G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Spatial Data/Spatially Derived/Slopes/Slp_",name,d8ordinfstr,"_",watershedname),format = "GTiff",overwrite = TRUE)
writeRaster(sca,paste0("G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Spatial Data/Spatially Derived/SCAs/sca_",name,d8ordinfstr,"_",watershedname),format = "GTiff",overwrite = TRUE)
writeRaster(p,paste0("G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Spatial Data/Spatially Derived/Aspects/aspect_",name,d8ordinfstr,"_",watershedname),format = "GTiff",overwrite = TRUE)


print(diroptions[DEMchoice])
print(d8ordinfstr)
print(watershedname)
beep("mario")
stop()
#Higher res
DEMchoice = 4; d8ordinf = 2
#Lower res
DEMchoice2 = 3; d8ordinf2 = 2
#Compare SCAs to other SCAs
d8ordinfstr3 = c("_D8","_Dinf")
name3 = c("1arcsec","1_3arcsec","LIDAR_2010","LIDAR_2018")
name4 = c("1arcsec","1_3arcsec","2010LIDAR","2018LIDAR")
url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/SCA_",name3[DEMchoice],"/SCA_",name3[DEMchoice],d8ordinfstr3[d8ordinf],"/SCA_",name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")
url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/SCA_",name3[DEMchoice2],"/SCA_",name3[DEMchoice2],d8ordinfstr3[d8ordinf2],"/SCA_",name4[DEMchoice2],d8ordinfstr3[d8ordinf2],".tif")
url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[DEMchoice],"/TIClass_",name3[DEMchoice],d8ordinfstr3[d8ordinf],"/TIClass_",name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")
url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[DEMchoice2],"/TIClass_",name3[DEMchoice2],d8ordinfstr3[d8ordinf2],"/TIClass_",name4[DEMchoice2],d8ordinfstr3[d8ordinf2],".tif")
rast=raster(url); rast2 = raster(url2)

#Comepare DEMs to other DEMs
folder = paste(diroptions[DEMchoice],sep = "")
url=paste0(basedir,"/",folder,fileoptions[DEMchoice])
rast=raster(url); if(DEMchoice==3||DEMchoice==4){values(rast) = values(rast)*0.3048}
folder = paste(diroptions[DEMchoice2],sep = "")
url=paste0(basedir,"/",folder,fileoptions[DEMchoice2])
rast2 = raster(url); if(DEMchoice2==3||DEMchoice2==4){values(rast2) = values(rast2)*0.3048}

#Extract rast2 values for each raster in rast
x = SpatialPoints(rast)
rast2altered = raster::extract(rast2,x)
rast2alteredrast = rast
values(rast2alteredrast) = rast2altered
plot((rast2alteredrast-rast),ylim = c(4117000,4118800),xlim = c(549050,551600))
points(outlet)
#zoom(rast2alteredrast-rast)
hist((rast2alteredrast-rast))
sd(values(rast2alteredrast-rast),na.rm = TRUE)

corrtf = TRUE
if(corrtf){
#Paper figure: dems to eachother
  par(mar = c(4,1,1,1))
  par(mfrow = c(1,3))
  construction = readOGR("G:/My Drive/Which DEM Paper/SWATHelperFiles/Construction_exclusion.shp")
  constrast = 
for(j in 1:3){
whichfigcor = j #1 for TIC, 2 for Slp, 3 for ln(SCA)
corspatial = NULL
url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[4],"/TIClass_",name3[4],d8ordinfstr3[1],"/TIClass_",name4[4],d8ordinfstr3[1],".tif")
rast=raster(url); rast = mask(rast,construction,inverse = TRUE) 
x = SpatialPoints(rast)
for (i in 1:8){
  #lower resolution is DEMchoice2
  if(i == 1){DEMchoice = 4; d8ordinf = 2}
  if(i == 2){DEMchoice = 4; d8ordinf = 1}
  if(i == 3){DEMchoice = 3; d8ordinf = 2}
  if(i == 4){DEMchoice = 3; d8ordinf = 1}
  if(i == 5){DEMchoice = 2; d8ordinf = 2}
  if(i == 6){DEMchoice = 2; d8ordinf = 1}
  if(i == 7){DEMchoice = 1; d8ordinf = 2}
  if(i == 8){DEMchoice = 1; d8ordinf = 1}
  if(whichfigcor==3){url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[DEMchoice],"/TIClass_",name3[DEMchoice],d8ordinfstr3[d8ordinf],"/TIClass_",name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")}
  if(whichfigcor==1){url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/Slope_",name3[DEMchoice],"/Slope_",name3[DEMchoice],d8ordinfstr3[d8ordinf],"/Slope_",name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")}
  if(whichfigcor==2){url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/SCA_",name3[DEMchoice],"/SCA_",name3[DEMchoice],d8ordinfstr3[d8ordinf],"/SCA_",name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")}
  rast2 = raster(url2); rast2 = mask(rast2,construction,inverse = TRUE)
  hold = raster::extract(rast2,x)
  if(whichfigcor==3){hold = log(hold)}
  if(i == 1){corspatial$LI18Dinf = hold}
  if(i == 2){corspatial$LI18D8 = hold}
  if(i == 3){corspatial$LI10Dinf = hold}
  if(i == 4){corspatial$LI10D8 = hold}
  if(i == 5){corspatial$as13Dinf = hold}
  if(i == 6){corspatial$as13D8 = hold}
  if(i == 7){corspatial$as1Dinf = hold}
  if(i == 8){corspatial$as1D8 = hold}
  print(i)
}
beep("mario")
par(mar = c(4,1,1,1))
cor4plotspatial = matrix(NA,8,8); pval4plotspatial = matrix(NA,8,8)
cor4plotspatial[1,2] = cor(corspatial$LI18Dinf,corspatial$LI18D8,use = "complete.obs"); cor4plotspatial[1,3] = cor(corspatial$LI18Dinf,corspatial$LI10Dinf,use = "complete.obs"); cor4plotspatial[1,4] = cor(corspatial$LI18Dinf,corspatial$LI10D8,use = "complete.obs"); 
cor4plotspatial[1,5] = cor(corspatial$LI18Dinf,corspatial$as13Dinf,use = "complete.obs"); cor4plotspatial[1,6] = cor(corspatial$LI18Dinf,corspatial$as13D8,use = "complete.obs"); cor4plotspatial[1,7] = cor(corspatial$LI18Dinf,corspatial$as1Dinf,use = "complete.obs"); cor4plotspatial[1,8] = cor(corspatial$LI18Dinf,corspatial$as1D8,use = "complete.obs"); 
cor4plotspatial[2,3] = cor(corspatial$LI18D8,corspatial$LI10Dinf,use = "complete.obs"); cor4plotspatial[2,4] = cor(corspatial$LI18D8,corspatial$LI10D8,use = "complete.obs"); cor4plotspatial[2,5] = cor(corspatial$LI18D8,corspatial$as13Dinf,use = "complete.obs");
cor4plotspatial[2,6] = cor(corspatial$LI18D8,corspatial$as13D8,use = "complete.obs"); cor4plotspatial[2,7] = cor(corspatial$LI18D8,corspatial$as1Dinf,use = "complete.obs"); cor4plotspatial[2,8] = cor(corspatial$LI18D8,corspatial$as1D8,use = "complete.obs");
cor4plotspatial[3,4] = cor(corspatial$LI10Dinf,corspatial$LI10D8,use = "complete.obs"); cor4plotspatial[3,5] = cor(corspatial$LI10Dinf,corspatial$as13Dinf,use = "complete.obs"); cor4plotspatial[3,6] = cor(corspatial$LI10Dinf,corspatial$as13D8,use = "complete.obs"); cor4plotspatial[3,7] = cor(corspatial$LI10Dinf,corspatial$as1Dinf,use = "complete.obs"); cor4plotspatial[3,8] = cor(corspatial$LI10Dinf,corspatial$as1D8,use = "complete.obs");
cor4plotspatial[4,5] = cor(corspatial$LI10D8,corspatial$as13Dinf,use = "complete.obs"); cor4plotspatial[4,6] = cor(corspatial$LI10D8,corspatial$as13D8,use = "complete.obs"); cor4plotspatial[4,7] = cor(corspatial$LI10D8,corspatial$as1Dinf,use = "complete.obs"); cor4plotspatial[4,8] = cor(corspatial$LI10D8,corspatial$as1D8,use = "complete.obs");
cor4plotspatial[5,6] = cor(corspatial$as13Dinf,corspatial$as13D8,use = "complete.obs"); cor4plotspatial[5,7] = cor(corspatial$as13Dinf,corspatial$as1Dinf,use = "complete.obs"); cor4plotspatial[5,8] = cor(corspatial$as13Dinf,corspatial$as1D8,use = "complete.obs");
cor4plotspatial[6,7] = cor(corspatial$as13D8,corspatial$as1Dinf,use = "complete.obs"); cor4plotspatial[6,8] = cor(corspatial$as13D8,corspatial$as1D8,use = "complete.obs");
cor4plotspatial[7,8] = cor(corspatial$as1Dinf,corspatial$as1D8,use = "complete.obs");
# pval4plotspatial[1,2] = cor.test(corspatial$LI18Dinf,corspatial$LI18D8,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[1,3] = cor.test(corspatial$LI18Dinf,corspatial$LI10Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[1,4] = cor.test(corspatial$LI18Dinf,corspatial$LI10D8,method = "pearson",use = "complete.obs")$p.value; 
# pval4plotspatial[1,5] = cor.test(corspatial$LI18Dinf,corspatial$as13Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[1,6] = cor.test(corspatial$LI18Dinf,corspatial$as13D8,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[1,7] = cor.test(corspatial$LI18Dinf,corspatial$as1Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[1,8] = cor.test(corspatial$LI18Dinf,corspatial$as1D8,method = "pearson",use = "complete.obs")$p.value; 
# pval4plotspatial[2,3] = cor.test(corspatial$LI18D8,corspatial$LI10Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[2,4] = cor.test(corspatial$LI18D8,corspatial$LI10D8,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[2,5] = cor.test(corspatial$LI18D8,corspatial$as13Dinf,method = "pearson",use = "complete.obs")$p.value;
# pval4plotspatial[2,6] = cor.test(corspatial$LI18D8,corspatial$as13D8,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[2,7] = cor.test(corspatial$LI18D8,corspatial$as1Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[2,8] = cor.test(corspatial$LI18D8,corspatial$as1D8,method = "pearson",use = "complete.obs")$p.value;
# pval4plotspatial[3,4] = cor.test(corspatial$LI10Dinf,corspatial$LI10D8,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[3,5] = cor.test(corspatial$LI10Dinf,corspatial$as13Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[3,6] = cor.test(corspatial$LI10Dinf,corspatial$as13D8,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[3,7] = cor.test(corspatial$LI10Dinf,corspatial$as1Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[3,8] = cor.test(corspatial$LI10Dinf,corspatial$as1D8,method = "pearson",use = "complete.obs")$p.value;
# pval4plotspatial[4,5] = cor.test(corspatial$LI10D8,corspatial$as13Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[4,6] = cor.test(corspatial$LI10D8,corspatial$as13D8,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[4,7] = cor.test(corspatial$LI10D8,corspatial$as1Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[4,8] = cor.test(corspatial$LI10D8,corspatial$as1D8,method = "pearson",use = "complete.obs")$p.value;
# pval4plotspatial[5,6] = cor.test(corspatial$as13Dinf,corspatial$as13D8,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[5,7] = cor.test(corspatial$as13Dinf,corspatial$as1Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[5,8] = cor.test(corspatial$as13Dinf,corspatial$as1D8,method = "pearson",use = "complete.obs")$p.value;
# pval4plotspatial[6,7] = cor.test(corspatial$as13D8,corspatial$as1Dinf,method = "pearson",use = "complete.obs")$p.value; pval4plotspatial[6,8] = cor.test(corspatial$as13D8,corspatial$as1D8,method = "pearson",use = "complete.obs")$p.value;
# pval4plotspatial[7,8] = cor.test(corspatial$as1Dinf,corspatial$as1D8,method = "pearson",use = "complete.obs")$p.value;
colnames(cor4plotspatial) = c("2018 LI Dinf","2018 LI D8","2010 LI Dinf", "2010 LI D8","1/3 as Dinf","1/3 as D8","1 as Dinf","1 as D8")
rownames(cor4plotspatial) = c("2018 LI Dinf","2018 LI D8","2010 LI Dinf", "2010 LI D8","1/3 as Dinf","1/3 as D8","1 as Dinf","1 as D8")
titleforcorplot = c("(a) Slope (m/m)","(b) ln(SCA)","(c) TIC")
corrplot(cor4plotspatial, tl.col = "black",type = "upper",title = titleforcorplot[whichfigcor],diag = FALSE,mar=c(0,0,2.5,0))
}}

#Comepare DEMs to other DEMs/ TICs to TICs
whichgraphic = 2 #=1 when you want to compare raw DEMs; =2 when you want to compare D8 TICs; =3 when you want to compare Dinf TICs
if(whichgraphic==1){d8ordinf = NULL}
if(whichgraphic==2){d8ordinf = 1}
if(whichgraphic==3){d8ordinf = 2}
a = NULL; b = NULL; c = NULL; d = NULL; e = NULL; f = NULL;
for (i in 1:6){
  #lower resolution is DEMchoice2
  if(i == 1){DEMchoice = 2; DEMchoice2 = 1}
  if(i == 2){DEMchoice = 3; DEMchoice2 = 1}
  if(i == 3){DEMchoice = 4; DEMchoice2 = 1}
  if(i == 4){DEMchoice = 3; DEMchoice2 = 2}
  if(i == 5){DEMchoice = 4; DEMchoice2 = 2}
  if(i == 6){DEMchoice = 4; DEMchoice2 = 3}
  if(DEMchoice<DEMchoice2){stop()}
  if(whichgraphic == 1){
    folder = paste(diroptions[DEMchoice],sep = "")
    url=paste0(basedir,"/",folder,fileoptions[DEMchoice])
    rast=raster(url); if(DEMchoice==3||DEMchoice==4){values(rast) = values(rast)*0.3048}
    folder = paste(diroptions[DEMchoice2],sep = "")
    url=paste0(basedir,"/",folder,fileoptions[DEMchoice2])
    rast2 = raster(url); if(DEMchoice2==3||DEMchoice2==4){values(rast2) = values(rast2)*0.3048}
  }else{
    url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[DEMchoice],"/TIClass_",name3[DEMchoice],d8ordinfstr3[d8ordinf],"/TIClass_",name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")
    url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[DEMchoice2],"/TIClass_",name3[DEMchoice2],d8ordinfstr3[d8ordinf],"/TIClass_",name4[DEMchoice2],d8ordinfstr3[d8ordinf],".tif")
    rast=raster(url); rast2 = raster(url2);
  }
  x = SpatialPoints(rast)
  rast2altered = raster::extract(rast2,x)
  rast2alteredrast = rast
  values(rast2alteredrast) = rast2altered
  if(i == 1){a = rast2alteredrast-rast;values(a) = round(values(a),1)}
  if(i == 2){b = rast2alteredrast-rast;values(b) = round(values(b),1)}
  if(i == 3){c = rast2alteredrast-rast;values(c) = round(values(c),1)}
  if(i == 4){d = rast2alteredrast-rast;values(d) = round(values(d),1)}
  if(i == 5){e = rast2alteredrast-rast;values(e) = round(values(e),1)}
  if(i == 6){f = rast2alteredrast-rast;values(f) = round(values(f),1)}
  print(i)
}
beep("mario")
as1 = NULL; as13 = NULL; LI10 = NULL; LI18 = NULL;
for (i in 1:4){
  url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/SCA_",name3[i],"/SCA_",name3[i],d8ordinfstr3[1],"/SCA_",name4[i],d8ordinfstr3[1],".tif")
  rast = raster(url); values(rast)[!is.na(values(rast))] = 0
  if(i ==1){as1 = rasterToPolygons(rast, dissolve=TRUE,na.rm = TRUE)}
  if(i ==2){rast = aggregate(rast,2);as13 = rasterToPolygons(rast, dissolve=TRUE,na.rm = TRUE)}
  if(i ==3){rast = aggregate(rast,10);LI10 = rasterToPolygons(rast, dissolve=TRUE,na.rm = TRUE)}
  if(i ==4){rast = aggregate(rast,20);LI18 = rasterToPolygons(rast, dissolve=TRUE,na.rm = TRUE)}
  print(i)
}
beep("coin")
exty = c(4117200,4118590); extx =  c(549150,551470)
color = brewer.pal(5,"Set1")
color = color[2:5]
colorforTIC = matrix("grey",10,1)
colorforTIC = c(brewer.pal(9,"Oranges")[9:1],"grey100",brewer.pal(9,"BuGn"))
colorforDEM = matrix("grey",17,1)
colorforDEM = paste0(colorforDEM,c(ceiling(seq(0,100,length.out = 8)),ceiling(seq(86,0,length.out = 7))))

if(whichgraphic==1){
  values(a)[which(values(a)>8)] = 8;values(a)[which(values(a)< -8)] = -8;
  values(b)[which(values(b)>8)] = 8;values(b)[which(values(b)< -8)] = -8;
  values(c)[which(values(c)>8)] = 8;values(c)[which(values(c)< -8)] = -8;
  values(d)[which(values(d)>8)] = 8;values(d)[which(values(d)< -8)] = -8;
  values(e)[which(values(e)>8)] = 8;values(e)[which(values(e)< -8)] = -8;
  values(f)[which(values(f)>8)] = 8;values(f)[which(values(f)< -8)] = -8;
  zlimit = c(-8,8)
  naturallog = TRUE; if(naturallog){a = log(abs(a));b = log(abs(b));c = log(abs(c));d = log(abs(d));e = log(abs(e));f = log(abs(f));zlimit = c(0,2.08)}
  d8 = FALSE; if(d8){d8ordinfstr="d8"}else{d8ordinfstr="dinf"}
  river1as = readOGR(paste0("G:/My Drive/Which DEM Paper/SWATHelperFiles/Rivers/1arcsec/",d8ordinfstr,"/hold.shp"))
  river13as = readOGR(paste0("G:/My Drive/Which DEM Paper/SWATHelperFiles/Rivers/1_3arcsec/",d8ordinfstr,"/hold.shp"))
  river2010LI = readOGR(paste0("G:/My Drive/Which DEM Paper/SWATHelperFiles/Rivers/LIDAR_2010/",d8ordinfstr,"/hold.shp"))
  river2018LI = readOGR(paste0("G:/My Drive/Which DEM Paper/SWATHelperFiles/Rivers/LIDAR_2018/",d8ordinfstr,"/hold.shp"))
  river2018LIdinf = readOGR(paste0("G:/My Drive/Which DEM Paper/SWATHelperFiles/Rivers/LIDAR_2018/dinf/hold.shp"))
#plot(as1)
#plot(as13)
#plot(LI10)
#plot(LI18)
m <- do.call(bind, list(as1,as13)); m = gUnaryUnion(m); a = mask(a,m)
m <- do.call(bind, list(as1,LI10)); m = gUnaryUnion(m); b = mask(b,m)
m <- do.call(bind, list(as1,LI18)); m = gUnaryUnion(m); c = mask(c,m)
m <- do.call(bind, list(as13,LI10)); m = gUnaryUnion(m); d = mask(d,m)
m <- do.call(bind, list(as13,LI18)); m = gUnaryUnion(m); e = mask(e,m)
m <- do.call(bind, list(LI10,LI18)); m = gUnaryUnion(m); f = mask(f,m)

par(mfrow = c(4,4))
par(mar = c(2.5,1,1,4))

#1
plot(0,0,col = "white", xaxt = "n", yaxt = "n")
text(0,0,label = "2018 LI",cex = 3)
#2
plot(f, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000)
mtext("m", side=4,cex = .6,line = 3,las=2)
plot(river2010LI,pch = 18,add = TRUE, col = color[2],cex = .35)
plot(river2018LI,pch = 18,add = TRUE, col = color[1],cex = .25)
lines(LI18, col = color[1])
lines(LI10, col = color[2])
points(outlet,pch = 16,cex = 1,col = 'darkgreen')

#3
plot(e, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000)
mtext("m", side=4,cex = .6,line = 3,las=2)
plot(river13as,pch = 18,add = TRUE, col = color[3],cex = .35)
plot(river2018LI,pch = 18,add = TRUE, col = color[1],cex = .25)
lines(LI18, col = color[1])
lines(as13,col = color[3])
points(outlet,pch = 16,cex = 1,col = 'darkgreen')

#4
plot(c, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000)
mtext("m", side=4,cex = .6,line = 3,las=2)
plot(river1as,pch = 18,add = TRUE, col = color[4],cex = .35)
plot(river2018LI,pch = 18,add = TRUE, col = color[1],cex = .25)
lines(as1,col = color[4])
lines(LI18, col = color[1])
points(outlet,pch = 16,cex = 1,col = 'darkgreen')

#5
hist(f, yaxt = "n",main = "",xlim = c(-7.5,7.5))

#6
plot(0,0,col = "white", xaxt = "n", yaxt = "n")
text(0,0,label = "2010 LI",cex = 3)

#7
plot(d, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000)
mtext("m", side=4,cex = .6,line = 3,las=2)
plot(river13as,pch = 18,add = TRUE, col = color[3],cex = .35)
plot(river2010LI,pch = 18,add = TRUE, col = color[2],cex = .25)
lines(as13,col = color[3])
lines(LI10, col = color[2])
points(outlet,pch = 16,cex = 1,col = 'darkgreen')

#8
plot(b, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000)
mtext("m", side=4,cex = .6,line = 3,las=2)
plot(river1as,pch = 15,add = TRUE, col = color[4],cex = .45)
plot(river2010LI,pch = 18,add = TRUE, col = color[2],cex = .25)
lines(as1,col = color[4])
lines(LI10, col = color[2])
points(outlet,pch = 16,cex = 1,col = 'darkgreen')

#9
hist(e, yaxt = "n",main = "",xlim = c(-7.5,7.5))

#10 
hist(d, yaxt = "n",main = "",xlim = c(-7.5,7.5))

#11
plot(0,0,col = "white", xaxt = "n", yaxt = "n")
text(0,0,label = "1/3 as",cex = 3)

#12
plot(a, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000)
mtext("m", side=4,cex = .6,line = 3,las=2)
plot(river1as,pch = 15,add = TRUE, col = color[4],cex = .45)
plot(river13as,pch = 18,add = TRUE, col = color[3],cex = .25)
lines(as1,col = color[4])
lines(as13,col = color[3])
points(outlet,pch = 16,cex = 1,col = 'darkgreen')

#13
hist(c, yaxt = "n",main = "",xlim = c(-7.5,7.5))
#14
hist(b, yaxt = "n",main = "",xlim = c(-7.5,7.5))
#15
hist(a, yaxt = "n",main = "",xlim = c(-7.5,7.5))
#16
plot(0,0,col = "white", xaxt = "n", yaxt = "n")
text(0,0,label = "1 as",cex = 3)


# #plot zoomed in version of lidar comparison
# par(mfrow = c(1,2))
# par(mar = c(1,1,2,5))
# plot(f, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,main = "2010LI-2018LI",zlim = c(-10,10),col = colorforTIC)
# lines(LI18,col = color[1])
# lines(LI10, col = color[2])
# points(outlet,pch = 16,cex = 1,col = 'darkgreen')
# plot(f, xaxt = "n", yaxt = "n",ylim = c(4117500,4118550),xlim = c(549000,550450),main = "2010LI-2018LI zoomed")
# lines(LI18,col = color[1])
# lines(LI10, col = color[2])
# points(outlet,pch = 16,cex = 1,col = 'darkgreen')
# 
# par(mfrow = c(1,1))
# newf = f; values(newf)[which(abs(values(newf))<=.25)] = 0; values(newf)[which(abs(values(newf))>.25)] = .25
# plot(newf,col = colorRampPalette(c("white","red"))(2),main = "Values above 25 cm difference",xlim = extx,ylim = exty,yaxt = "n",xaxt = "n")
# lines(LI18,col = color[1])
# lines(LI10, col = color[2])
}
beep("mario")
if(whichgraphic==3|whichgraphic==2){
  par(mfrow = c(4,4))
  par(mar = c(3,1,2,3.5))
  #exty = c(4117000,4118800); extx =  c(549050,551600)
  #1
  plot(0,0,col = "white", xaxt = "n", yaxt = "n")
  text(0,0,label = "2018 LI",cex = 3)
  #2
  plot(f, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-10,10),col = colorforTIC)
  points(outlet,pch = 16,cex = 1,col = 'darkgreen')
  #3
  plot(e, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-10,10),col = colorforTIC)
  points(outlet,pch = 16,cex = 1,col = 'darkgreen')
  #4
  plot(c, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-10,10),col = colorforTIC)
  points(outlet,pch = 16,cex = 1,col = 'darkgreen')
  #5
  hist(f, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(f),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-10,10))
  #6
  plot(0,0,col = "white", xaxt = "n", yaxt = "n")
  text(0,0,label = "2010 LI",cex = 3)
  #7
  plot(d, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-10,10),col = colorforTIC)
  points(outlet,pch = 16,cex = 1,col = 'darkgreen')
  #8
  plot(b, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-10,10),col = colorforTIC)
  points(outlet,pch = 16,cex = 1,col = 'darkgreen')
  #9
  hist(e, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(e),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-10,10))
  #10 
  hist(d, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(d),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-10,10))
  #11
  plot(0,0,col = "white", xaxt = "n", yaxt = "n")
  text(0,0,label = "1/3 as",cex = 3)
  #12
  plot(a, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-10,10),col = colorforTIC)
  points(outlet,pch = 16,cex = 1,col = 'darkgreen')
  lines(as1,col = color[4])
  lines(as13,col = color[3])
  lines(LI10, col = color[2])
  lines(LI18,col = color[1])
  #13
  hist(c, yaxt = "n",main = paste("\n SD:", sprintf("%.2f",round(sd(values(c),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-10,10))
  #14
  hist(b, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(b),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-10,10))
  #15
  hist(a, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(a),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-10,10))
  #16
  plot(0,0,col = "white", xaxt = "n", yaxt = "n")
  text(0,0,label = "1 as",cex = 3)
}
beep("mario")
# if(whichgraphic==2){
#   par(mfrow = c(2,2))
#   par(mar = c(3,1,2,3))
#   #exty = c(4117000,4118800); extx =  c(549050,551600)
#   #1
#   plot(0,0,col = "white", xaxt = "n", yaxt = "n")
#   text(0,0,label = "2018 LI",cex = 3)
#   #2
#   plot(f, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-10,10),col = gray.colors(40, start = 0.01, end = 1, gamma = 2.2, alpha = NULL))
#   points(outlet,pch = 16,cex = 1,col = 'darkgreen')
#   #3
#   hist(f, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(f),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-10,10))
#   #4
#   plot(0,0,col = "white", xaxt = "n", yaxt = "n")
#   text(0,0,label = "2010 LI",cex = 3)
# }

#Comepare TICs to other TICs for D8 and Dinf varying ones
if(FALSE){
a = NULL; b = NULL; c = NULL; d = NULL; e = NULL; f = NULL;
for (i in 1:4){
  #lower resolution is DEMchoice2
  if(i == 1){DEMchoice = 1; d8ordinf = 1; d8ordin2 = 2}
  if(i == 2){DEMchoice = 2; d8ordinf = 1; d8ordin2 = 2}
  if(i == 3){DEMchoice = 3; d8ordinf = 1; d8ordin2 = 2}
  if(i == 4){DEMchoice = 4; d8ordinf = 1; d8ordin2 = 2}
  url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[DEMchoice],"/TIClass_",name3[DEMchoice],d8ordinfstr3[d8ordinf],"/TIClass_",name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")
  url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[DEMchoice],"/TIClass_",name3[DEMchoice],d8ordinfstr3[d8ordinf2],"/TIClass_",name4[DEMchoice],d8ordinfstr3[d8ordinf2],".tif")
  rast=raster(url); rast2 = raster(url2);
  if(i == 1){a = rast2-rast}
  if(i == 2){b = rast2-rast}
  if(i == 3){c = rast2-rast}
  if(i == 4){d = rast2-rast}
  print(i)
}


par(mfrow = c(4,2))
par(mar = c(2,2,5,2))
exty = c(4117000,4118800); extx =  c(549050,551600)
#1
hist(d,main = paste("\n SD:",round(sd(values(a),na.rm = TRUE),2)), yaxt = "n")
#2
plot(d, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,main = "\n TIC Differences 2018 LI",col = colorforTIC,zlim = c(-10,10))
points(outlet,pch = 16,col = 'darkgreen',cex = 1)

#3
hist(c,main = paste("\n SD:",round(sd(values(b),na.rm = TRUE),2)), yaxt = "n")
#4
plot(c, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,main = "\n TIC Differences 2010 LI",col = colorforTIC,zlim = c(-10,10))
points(outlet,pch = 16,col = 'darkgreen',cex = 1)

#5
hist(b,main = paste("\n SD:",round(sd(values(c),na.rm = TRUE),2)), yaxt = "n")
#6
plot(b, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,main = "\n TIC Differences 1/3 as",col = colorforTIC,zlim = c(-10,10))
points(outlet,pch = 16,col = 'darkgreen',cex = 1)

#7
hist(a,main = paste("\n SD:",round(sd(values(d),na.rm = TRUE),2)), yaxt = "n")
#8
plot(a, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,main = "\n TIC Differences 1 as",col = colorforTIC,zlim = c(-10,10))
points(outlet,pch = 16,col = 'darkgreen',cex = 1)
lines(as1,col = color[4])
lines(as13,col = color[3])
lines(LI10, col = color[2])
lines(LI18,col = color[1])
}
beep()

