# The purpose of this code is to extract calculated spatial properties for 
#rasters (generated via the "runDelin_taudem.R") for each location of soil sample

pacman::p_load(devtools)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(rgdal,googlesheets4,raster)
pacman::p_load(elevatr,raster,soilDB,rgdal,readxl)

#Options for file names
sheetname = c("TIC", "TIV", "Slp", "SCA")
sheetname2 = c("TIC", "TIV", "Slp", "sca")
resolution = c("1arcsec","1_3arcsec","1meter")
watershed = c("_LOC","_DC")
methodoptions = c("d8","dinf")

#Load Sheet
#sheetid="1Ny241_GTJ_8S1wvQHF2G8Ev9I1l9p4QTIlPa_a4d9GA"
#TIC=read_sheet(sheetid,sheet = "TIC") #unitless
#TIV=read_sheet(sheetid,sheet = "TIV") #unitless
#SCA=read_sheet(sheetid,sheet = "SCA") #meters
#Slope=read_sheet(sheetid,sheet = "Slope") #length/length
soils <- read_excel("2020 Soil Health Raw Data_50 samples (003).xlsx")

TIC = data.frame(soilsamp = soils$SampleNumber, lat = soils$Latitude,long = soils$Longitude)
TIV = TIC; Slp = TIC; SCA = TIC

for (i in 1:length(sheetname)){
  for (j in 1:length(resolution)){
    for (k in 1:length(methodoptions)){
      #Get extraction raster
      url= paste(sheetname[i],"/",sheetname2[i],"_",resolution[j],methodoptions[k],watershed[1],".tif",sep = "")
      rast=raster(url)
      url= paste(sheetname[i],"/",sheetname2[i],"_",resolution[j],methodoptions[k],watershed[2],".tif",sep = "")
      rast2=raster(url)
      origin(rast2) = origin(rast)
      rast = merge(rast,rast2)
      plot(rast)
      
      proj4string(rast)
      proj4_utm = proj4string(rast)
      proj4_ll = "+proj=longlat"
      
      # Now we will build our proj4strings which define our “Coordinate 
      # Reference Systems” or CRS in future geographic manipulations. 
      crs_ll=CRS(proj4_ll)
      crs_utm=CRS(proj4_utm)
      
      
      #extraction
      SP_ll=SpatialPoints(cbind(TIC$long[!is.na(TIC$long)],TIC$lat[!is.na(TIC$lat)]),proj4string =crs_ll)
      if(sheetname[i] == "TIC"){
        TIC[!is.na(TIC$lat),paste(sheetname[i],"_",resolution[j],methodoptions[k],sep = "")]=raster::extract(rast,SP_ll)
      }else if(sheetname[i] == "Slp"){
        Slp[!is.na(TIC$lat),paste(sheetname[i],"_",resolution[j],methodoptions[k],sep = "")]=raster::extract(rast,SP_ll)
      }else if (sheetname[i] == "SCA"){
        SCA[!is.na(TIC$lat),paste(sheetname[i],"_",resolution[j],methodoptions[k],sep = "")]=raster::extract(rast,SP_ll)
      }else {
        TIV[!is.na(TIC$lat),paste(sheetname[i],"_",resolution[j],methodoptions[k],sep = "")]=raster::extract(rast,SP_ll)
      }
      
      
      # #Comparison
       # SP_ll=SpatialPoints(matrix(c(ExtractionSheet$Long,ExtractionSheet$Lat), 
       #                            ncol = 2, byrow = FALSE),proj4string =crs_ll)
       # SP_utm <- spTransform(SP_ll, crs_utm)
       # plot(rast)
       # points(SP_utm,pch = 24, cex=2, col="blue", bg="red", lwd=2)   #
       # raster::extract(rast,SP_utm)
       # ExtractionSheet[paste(Extraction)]=raster::extract(rast,SP_utm)
      print(i)
      print(j)
      print(k)
    }
  }
}

#at the moment I am copy and pasting things into google docs, not the most elegant but works in ARC
