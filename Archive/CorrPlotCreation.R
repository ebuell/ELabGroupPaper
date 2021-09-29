# The purpose of this code is largely to generate visualizations for correlations
#between measured and calculated spatial and physical properties

library('googlesheets4')
library('zoo')
library('dplyr')
library('corrplot')
library('plot.matrix')
library('grid')
library(readxl)

rm(list=ls())

#### read in data ####
spatial_sheetid="1Ny241_GTJ_8S1wvQHF2G8Ev9I1l9p4QTIlPa_a4d9GA"
SCA <- read_sheet(spatial_sheetid, sheet = "SCA")
Slope <- read_sheet(spatial_sheetid, sheet = "Slp")
TIC <- read_sheet(spatial_sheetid, sheet = "TIC")
TIV <- read_sheet(spatial_sheetid, sheet = "TIV")
soilsreadin = read_excel("G:/.shortcut-targets-by-id/15fVCtxOyan-kvhW9iOPkLv5ono_BCHpk/Lab/ELGroupPaperFall2021/Soil Physical Properties/2020 Soil Health Raw Data_50 samples (003).xlsx")
soilswcat = data.frame(ID = soilsreadin$SampleNumber,lat = soilsreadin$Latitude,long = soilsreadin$Longitude,
                   name = soilsreadin$SoilName,tilage = soilsreadin$Tillage_1to4,crplastyr = soilsreadin$CropLastYr,
                   crpthsyr = soilsreadin$CropThisYr,texture = soilsreadin$soil_texture_class,
                   sand = soilsreadin$soil_texture_sand,silt = soilsreadin$soil_texture_silt,
                   clay = soilsreadin$soil_texture_clay,predwater = soilsreadin$pred_water_capacity,
                   predwaterrating = soilsreadin$pred_water_capacity_rating,stability = soilsreadin$aggregate_stability,
                   stabilityrating = soilsreadin$aggregate_stability_rating,OM = soilsreadin$organic_matter,
                   OMrating = soilsreadin$organic_matter_rating,TC = soilsreadin$total_c,
                   activecarbon = soilsreadin$active_carbon,activecarbonrating = soilsreadin$active_carbon_rating,TN = soilsreadin$total_n,
                   soilprotein = soilsreadin$pred_soil_protein,soilproteinrating = soilsreadin$pred_soil_protein_rating,
                   respiration = soilsreadin$respiration,respirationrating = soilsreadin$respiration_rating,
                   ph = soilsreadin$ph,phrating = soilsreadin$ph_rating,phos = soilsreadin$p,phosreading = soilsreadin$p_rating,
                   potas = soilsreadin$k,potasrating = soilsreadin$k_rating,mg = soilsreadin$mg,fe = soilsreadin$fe,
                   mn = soilsreadin$mn,zn = soilsreadin$zn,minorelrating = soilsreadin$minor_elements_rating,
                   minoroverallscore = soilsreadin$`overall score`,al = soilsreadin$`Mod. Morgan Al. ppm`,
                   ca = soilsreadin$`Mod. Morgan Ca. ppm`,cu = soilsreadin$`Mod. Morgan Cu. ppm`,sulf = soilsreadin$`Mod. Morgan S. ppm`,
                   b = soilsreadin$`Mod. Morgan B. ppm`)
reduced = TRUE
if(reduced){soils = soilswcat[,c(1:3,9:12,14,16,21:22,24,26,28,30,37)]
}else{soils = soilswcat[,c(1:3,9:42)]}
#####

#### Clean Spatial Data ####
for(i in 4:length(names(Slope))){
  Slope[,i] = as.numeric(unlist(Slope[,i]))
  SCA[,i] = as.numeric(unlist(SCA[,i]))
  TIV[,i] = as.numeric(unlist(TIV[,i]))
  TIC[,i] = as.numeric(unlist(TIC[,i]))
}


### create correlation matrix and visualize correlations ###
cor4plot = data.frame(slp_1asd8 = NA,lnsca_1asd8 = NA,tiv_1asd8 = NA,tic_1asd8 = NA,
                      slp_1asdinf = NA,lnsca_1asdinf = NA,tiv_1asdinf = NA,tic_1asdinf = NA,
                      slp_1_3asd8 = NA,lnsca_1_3asd8 = NA,tiv_1_3asd8 = NA,tic_1_3asd8 = NA,
                      slp_1_3asdinf = NA,lnsca_1_3asdinf = NA,tiv_1_3asdinf = NA,tic_1_3asdinf = NA,
                      slp_1md8 = NA,lnsca_1md8 = NA,tiv_1md8 = NA,tic_1md8 = NA,
                      slp_1mdinf = matrix(NA,length(names(soils))-3,1),lnsca_1mdinf = NA,tiv_1mdinf = NA,tic_1mdinf = NA)
rownames(cor4plot) = names(soils)[4:length(names(soils))]
par(mfrow = c(6,4))
par(mar = c(4,4,2,1))
for(i in 4:length(names(soils))){
    for(j in 4:(length(names(Slope)))){
      cor4plot[i-3,1+(j-4)*4] = cor(soils[,i],Slope[,j],use = "complete.obs")
      plot(soils[!is.na(Slope[,j]),i],pull(Slope[!is.na(Slope[,j]),j]),xlab = names(soils)[i],ylab = names(Slope)[j],main = paste("Cor:",round(cor4plot[i-3,1+(j-4)*4],2)))
      cor4plot[i-3,2+(j-4)*4] = cor(soils[,i],log(SCA[,j]+0.000001),use = "complete.obs")
      plot(soils[!is.na(SCA[,j]),i],log(pull(SCA[!is.na(SCA[,j]),j])),xlab = names(soils)[i],ylab = paste("ln",names(SCA)[j]),main = paste("Cor:",round(cor4plot[i-3,2+(j-4)*4],2)))
      cor4plot[i-3,3+(j-4)*4] = cor(soils[,i],TIC[,j],use = "complete.obs")
      plot(soils[!is.na(TIV[,j]),i],pull(TIV[!is.na(TIV[,j]),j]),xlab = names(soils)[i],ylab = names(TIV)[j],main = paste("Cor:",round(cor4plot[i-3,3+(j-4)*4],2)))
      cor4plot[i-3,4+(j-4)*4] = cor(soils[,i],TIV[,j],use = "complete.obs")
      plot(soils[!is.na(TIC[,j]),i],pull(TIC[!is.na(TIC[,j]),j]),xlab = names(soils)[i],ylab = names(TIC)[j],main = paste("Cor:",round(cor4plot[i-3,4+(j-4)*4],2)))
      #print(paste("i=",i,"j=",j,"i ind=",i-3,"sCA j ind=",2+(j-4)*4))
    }
}

##### create corplot
#reorder cor4plot
cor4plot = cor4plot[,c(21:24,17:20,13:16,9:12,5:8,1:4)]

par(mfrow = c(1,1))
corrplot(as.matrix(cor4plot),tl.col = "black",mar=c(0,0,1,0),tl.cex=.75)
segments(.54,.5,.54,40); segments(4.54,.5,4.54,38,lty = 4); segments(8.54,.5,8.54,40); segments(12.54,.5,12.54,38,lty = 4); segments(16.54,.5,16.54,40); 
segments(20.54,.5,20.54,38,lty = 4); segments(24.54,.5,24.54,40);
if(!reduced){
  segments(.5,15.5,24.5,15.5); segments(.5,17.5,24.5,17.5,lty = 4); segments(.5,19.5,24.5,19.5,lty = 4);
  segments(.5,21.5,24.5,21.5,lty = 4); segments(.5,22.5,24.5,22.5); segments(.5,27.5,24.5,27.5);
  segments(.5,29.5,24.5,29.5); segments(.5,31.5,24.5,31.5);
}
