StudyArea<-extent(1.75,7.6,41.8,44.5)###Définition des contours de la zone d'étude
StudyArea2<-extent(0.34,12.83,40.66,45.21)###Définition des contours de la zone d'étude
StudyAreaAlberte<-extent(2.08333,7.58332975,42.16666667,45.00)
RasterAlberte<-raster(nrows=25, ncols=66, xmn=2.08330,xmx=7.58, ymn=42.16, ymx=45, resolution=0.08333333, crs=4326 )

###Importer les données la reconstuction climatique Guiot###
load("C:/Users/Bernigaud Nicolas/Documents/Guiot_data/reclim2_co2.RData")

###TEMPERATURES###
###Créer des rasters des anomalies###
ANO_T<-matrix(nrow = 81, ncol=14)
ANO_T[,13:14]<-refco2[,3:4]
Points_Numbers<-rownames(refco2)
Col_names<-c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "LAT", "LONG")
Centuries<-c("1300BP", "1400BP", "1500BP", "1600BP", "1700BP", "1800BP", "1900BP", "2000BP", "2100BP", "2200BP", "2300BP", "2400BP", "2500BP")
ANO_T_Ar<-array(data = ANO_T, dim =c(81,14,13), dimnames = list(Points_Numbers,Col_names,Centuries))

###boucle pour remplir l'array des valeurs d'anomalie de TEMPERATURE de reco2 entre 1300 BP et 2500 BP###
for (i in 1:13){
for (j in 1:12) {
    ANO_T_Ar[,j,i]<-reco2[i+13,,j+19]
}
}
###Convertir pour chaque siècle les matrices de l'array en SpatialPointDataFrame###
for (i in 1:13){
  assign(paste("ANO_T",i, sep = "_"), as.data.frame(ANO_T_Ar[,,i]))
  assign(paste("ANO_T", i, "SPDF", sep="_" ), SpatialPointsDataFrame(coords = get(paste("ANO_T",i,sep="_"))[,13:14], data= get(paste("ANO_T",i, sep="_"))[,1:12], proj4string = CRS("+init=epsg:4326")))
}
###Interpolation des SpatialPointsDataframe pour chaque siècle et rasterisation###
library(rgdal)
library(sp)
require(fields)
Coastline <- readOGR(dsn = 'C://Users/Bernigaud Nicolas/Documents/Article_SMA/LPJmL',
               layer = 'Europe_coastline_poly2')
CoastlineC<-crop(Coastline, extent(refco2[,3:4]))###nécessaire de couper ce fichiers très lourd, sinon les calculs sont très longs###
MyRaster3<-raster(ext=extent(StudyArea2), resolution=0.1)
Coord<-coordinates(ANO_T_1_SPDF)
Coord2<-refco2[,3:4]
for (i in 1:13){
for (j in 1:12){
  assign(paste("ANO_T",i, j, "Interpol", sep="_"),interpolate(MyRaster3, Tps(Coord2, get(paste("ANO_T", i, sep="_"))[,j])))
  assign(paste("ANO_T",i, j, "Interpol", "C", sep="_"), mask(get(paste("ANO","T",i, j, "Interpol", sep="_")), CoastlineC))
}
}
###Créer des rasterbrick pour chaque siècle###
for (i in 1:13){
assign(paste( "ANO_T_", i, sep=""),brick(
         list(
           get(paste("ANO_T_",i, "_1_Interpol_C", sep ="")),
           get(paste("ANO_T_",i, "_2_Interpol_C", sep ="")),
           get(paste("ANO_T_",i, "_3_Interpol_C", sep="")),
           get(paste("ANO_T_",i, "_4_Interpol_C", sep="")),
           get(paste("ANO_T_",i, "_5_Interpol_C", sep="")),
           get(paste("ANO_T_",i, "_6_Interpol_C", sep="")),
           get(paste("ANO_T_",i, "_7_Interpol_C", sep ="")),
           get(paste("ANO_T_",i, "_8_Interpol_C", sep ="")),
           get(paste("ANO_T_",i, "_9_Interpol_C", sep="")),
           get(paste("ANO_T_",i, "_10_Interpol_C", sep="")),
           get(paste("ANO_T_",i, "_11_Interpol_C", sep="")),
           get(paste("ANO_T_",i, "_12_Interpol_C", sep=""))
         )
       )
)
}
###Créer le raster de REF actuel des TEMPERATURES###
###TEMPERATURES MENSUELLES (= NORTAV) ####
Tempe_REF_mensuel<-read.table(file="https://www.data.gouv.fr/s/resources/indices-mensuels-de-temperature-et-nombre-de-jours-de-temperature-issus-du-modele-aladin-pour-la-periode-de-reference/20150930-200335/Tempe_REF_mensuel.txt", sep=";")
colnames(Tempe_REF_mensuel)<-c("Point", "Latitude", "Longitude", "Contexte", "Integration", "Mois", "NORTAV","NORTNAV", "NORTXAV", "NORTRAV", "NORTXQ90", "NORTXQ10","NORTNQ10","NORTNQ90","NORTXND","NORTNND","NORTNHT", "NORTXHWD", "NORTNCWD", "NORTNFD", "NORTXFD", "NORSD", "NORTR", "NORHDD", "NORCDD")
NORTAV<-matrix(nrow = 8602,ncol=14)
colnames(NORTAV)<-c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "LAT", "LONG")
seq1<-seq(1,103224, by=12)
for (k in 1:12){
  NORTAV[,k]<-Tempe_REF_mensuel[seq1+k-1,7]
}
NORTAV[,13]<-Tempe_REF_mensuel[seq1,2]
NORTAV[,14]<-Tempe_REF_mensuel[seq1,3]
NORTAV<-as.data.frame(NORTAV)
Ext_Nortav<-extent(min(NORTAV[,14]),max(NORTAV[,14]),min(NORTAV[,13]),max(NORTAV[,13]))
MyRaster<-raster(ext=Ext_Nortav, resolution=0.1)
REF_T<-rasterize(NORTAV[,14:13], MyRaster, NORTAV[1:12] )
REF_T<-crop(REF_T, StudyAreaAlberte)
REF_T<-resample(REF_T, RasterAlberte)
crs(REF_T)<-CRS("+init=epsg:4326")

###Ecrire le raster de REF des T###
writeRaster(REF_T, "C:/Users/Bernigaud Nicolas/Documents/Article_SMA/LPJmL/1300_2500BP_V2/REF/ALADIN_TEMP_REF.asc", bylayer=T)

###Resampler les raster d'anomalies d'après le modéle de REF###
for (i in 1:13){
assign(paste("ANO_T_",i, "_RS", sep = ""), resample(get(paste("ANO_T_", i, sep="")), REF_T))
}

###Additionner les anomalies de TEMPERATURE au modèle actuel###
for (i in 1:13){
  assign((paste("T_", 12 + i,"00BP", sep="")), get(paste("ANO_T_",i, "_RS", sep = "")) + REF_T)
}

###Ecrire les rasters###
for (i in 1:13){
  writeRaster(get(paste("T_", 12 + i, "00BP", sep="")), filename = paste0("Temp_",12 + i, "00BP.asc"), bylayer=T, overwrite=TRUE, crs=4326)
}

####PRECIPITATIONS####

###Création d'un array avec les valeurs d'anomalies par siècles###
ANO_P<-matrix(nrow = 81, ncol=14)
col_names<-c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12", "LAT", "LONG")
Centuries<-c("1300BP", "1400BP", "1500BP", "1600BP", "1700BP", "1800BP", "1900BP", "2000BP", "2100BP", "2200BP", "2300BP", "2400BP", "2500BP")
ANO_P[,13:14]<-refco2[,3:4]
ANO_P_Ar<-array(data = ANO_P, dim =c(81,14,13), dimnames = list(rownames(refco2),col_names,Centuries))

###boucle pour remplir l'array des valeurs d'anomalie de PRECIPITATIONS de reco2 entre 1300 BP et 2500 BP###
for (i in 1:13){
  for (j in 1:12) {
    ANO_P_Ar[,j,i]<-reco2[i+13,,j+31]
  }
}
###Convertir pour chaque siècle les matrices de l'array en SpatialPointDataFrame###
for (i in 1:13){
  assign(paste("ANO_P",i, sep = "_"), as.data.frame(ANO_P_Ar[,,i]))
  assign(paste("ANO_P", i, "SPDF", sep="_" ), SpatialPointsDataFrame(coords = get(paste("ANO_P",i,sep="_"))[,13:14], data= get(paste("ANO_P",i, sep="_"))[,1:12], proj4string = CRS("+init=epsg:4326")))
}
###Interpolation des SpatialPointsDataframe pour chaque siècle et rasterisation###
require(fields)
for (i in 1:13){
  for (j in 1:12){
    assign(paste("ANO_P",i, j, "Interpol", sep="_"),interpolate(MyRaster3, Tps(Coord, get(paste("ANO_P", i, sep="_"))[,j])))
    assign(paste("ANO_P",i, j, "Interpol", "C", sep="_"), mask(get(paste("ANO","P",i, j, "Interpol", sep="_")), CoastlineC))
  }
}
###Création de rasterbricks pour chaque siècle###
for (i in 1:13){
  assign(paste( "ANO_P_", i, sep=""),brick(
    list(
      get(paste("ANO_P_",i, "_1_Interpol_C", sep ="")),
      get(paste("ANO_P_",i, "_2_Interpol_C", sep ="")),
      get(paste("ANO_P_",i, "_3_Interpol_C", sep="")),
      get(paste("ANO_P_",i, "_4_Interpol_C", sep="")),
      get(paste("ANO_P_",i, "_5_Interpol_C", sep="")),
      get(paste("ANO_P_",i, "_6_Interpol_C", sep="")),
      get(paste("ANO_P_",i, "_7_Interpol_C", sep ="")),
      get(paste("ANO_P_",i, "_8_Interpol_C", sep ="")),
      get(paste("ANO_P_",i, "_9_Interpol_C", sep="")),
      get(paste("ANO_P_",i, "_10_Interpol_C", sep="")),
      get(paste("ANO_P_",i, "_11_Interpol_C", sep="")),
      get(paste("ANO_P_",i, "_12_Interpol_C", sep=""))
    )
  )
  )
}
###Importer et reformater les données de références sur la quantité de PRECIPITATIONS MENSUELLES en mm (= NORRR)####
Precip_REF_mensuel<-read.table(file="https://www.data.gouv.fr/s/resources/indices-mensuels-de-precipitations-et-nombre-de-jours-de-precipitations-issus-du-modele-aladin-pour-la-periode-de-reference/20150930-201506/Precip_REF_mensuel.txt", sep=";")
colnames(Precip_REF_mensuel)<-c("Point", "Latitude", "Longitude", "Contexte", "Integration", "Mois", "NORPAV", "NORPINT", "NORRR", "NORPFL90", "NORRR1MM", "NORPXCWD", "NORPN20MM", "NORPXCDD")
NORRR<-matrix(nrow = 8602,ncol=14)
colnames(NORRR)<-c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "LAT", "LONG")
seq1<-seq(1,103224, by=12)
for (i in 1:12){
  NORRR[,i]<-Precip_REF_mensuel[seq1+i-1,9]
}
NORRR[,13]<-Precip_REF_mensuel[seq1,2]
NORRR[,14]<-Precip_REF_mensuel[seq1,3]
NORRR_DF<-as.data.frame(NORRR)

###Créer le raster de REF actuel des PRECIPITATIONS###
Ext_NORRR<-extent(min(NORRR[,14]),max(NORRR[,14]),min(NORRR[,13]),max(NORRR[,13]))
MyRaster<-raster(ext=Ext_NORRR, resolution=0.1)
REF_P<-rasterize(NORRR_DF[,14:13], MyRaster, NORRR_DF[1:12] )
REF_P<-crop(REF_P, StudyAreaAlberte)
REF_P<-resample(REF_P, RasterAlberte)
crs(REF_T)<-CRS("+init=epsg:4326")

###Ecrire le raster de REF des P###
writeRaster(REF_P, "C:/Users/Bernigaud Nicolas/Documents/Article_SMA/LPJmL/1300_2500BP_V2/REF/ALADIN_PREC_REF.asc", bylayer=T)

###Resampler les raster d'anomalies d'après le modéle de REF###
for (i in 1:13){
  assign(paste("ANO_P_",i, "_RS", sep = ""), resample(get(paste("ANO_P_", i, sep="")), REF_P))
}

###Additionner les anomalies de PRECIPITATIONS au modèle actuel###
for (i in 1:13){
  assign((paste("P_", 12 + i,"00BP", sep="")), get(paste("ANO_P_",i, "_RS", sep = "")) + REF_P)
}

###Ecrire les rasters de PRECIPITATIONS###
for (i in 1:13){
  writeRaster(get(paste("P_", 12 + i, "00BP", sep="")), filename = paste0("Precip_",12 + i, "00BP.asc"), bylayer=T, proj4string = CRS("+init=epsg:4326"), overwrite=TRUE)
}

####SUNSHINE####

###Création d'un array avec les valeurs d'anomalies par siècles###
ANO_S<-matrix(nrow = 81, ncol=14)
colnames(ANO_S)<-c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "LAT", "LONG")
rownames(ANO_S)<-rownames(refco2)
Centuries<-c("1300BP", "1400BP", "1500BP", "1600BP", "1700BP", "1800BP", "1900BP", "2000BP", "2100BP", "2200BP", "2300BP", "2400BP", "2500BP")
ANO_S[,13:14]<-refco2[,3:4]
ANO_S_Ar<-array(data = ANO_S, dim =c(81,14,13), dimnames = list(rownames(ANO_S),colnames(ANO_S),Centuries))

###boucle pour remplir l'array des valeurs d'anomalie de SUNSHINE de reco2 entre 1300 BP et 2500 BP###
for (i in 1:13){
  for (j in 1:12) {
    ANO_S_Ar[,j,i]<-reco2[i+13,,j+43]
  }
}
###Convertir pour chaque siècle les matrices de l'array en SpatialPointDataFrame###
for (i in 1:13){
  assign(paste("ANO_S",i, sep = "_"), as.data.frame(ANO_S_Ar[,,i]))
  assign(paste("ANO_S", i, "SPDF", sep="_" ), SpatialPointsDataFrame(coords = get(paste("ANO_S",i,sep="_"))[,13:14], data= get(paste("ANO_S",i, sep="_"))[,1:12], proj4string = CRS("+init=epsg:4326")))
}
###Interpolation des SpatialPointsDataframe pour chaque siècle et rasterisation###
require(fields)
for (i in 1:13){
  for (j in 1:12){
    assign(paste("ANO_S",i, j, "Interpol", sep="_"),interpolate(MyRaster3, Tps(Coord, get(paste("ANO_S", i, sep="_"))[,j])))
    assign(paste("ANO_S",i, j, "Interpol", "C", sep="_"), mask(get(paste("ANO","S",i, j, "Interpol", sep="_")), CoastlineC))
  }
}
###Créer des rasterbrick pour chaque siècle###
for (i in 1:13){
  assign(paste( "ANO_S_", i, sep=""),brick(
    list(
      get(paste("ANO_S_",i, "_1_Interpol_C", sep ="")),
      get(paste("ANO_S_",i, "_2_Interpol_C", sep ="")),
      get(paste("ANO_S_",i, "_3_Interpol_C", sep="")),
      get(paste("ANO_S_",i, "_4_Interpol_C", sep="")),
      get(paste("ANO_S_",i, "_5_Interpol_C", sep="")),
      get(paste("ANO_S_",i, "_6_Interpol_C", sep="")),
      get(paste("ANO_S_",i, "_7_Interpol_C", sep ="")),
      get(paste("ANO_S_",i, "_8_Interpol_C", sep ="")),
      get(paste("ANO_S_",i, "_9_Interpol_C", sep="")),
      get(paste("ANO_S_",i, "_10_Interpol_C", sep="")),
      get(paste("ANO_S_",i, "_11_Interpol_C", sep="")),
      get(paste("ANO_S_",i, "_12_Interpol_C", sep=""))
    )
  )
  )
}

###Télécharger les données à 10' de résolution au lien suivant: https://crudata.uea.ac.uk/cru/data/hrg/tmc/ ###
CRU<-read.table(file="C:/Users/Bernigaud Nicolas/Documents/rasters_alberte_2020/grid_10min_sunp.dat/grid_10min_sunp.dat", sep="")
colnames(CRU)<-c("LAT_", "LONG_", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
CRU[c("LAT", "LONG")]<-NA
CRU[,15]<-CRU[,1]
CRU[,16]<-CRU[,2]
CRU<-CRU[,c(-1, -2)]
CRU<-as.data.frame(CRU)
CRU_SPDF<-SpatialPointsDataFrame(coords = CRU[,14:13],data = CRU[,1:12], proj4string = CRS('+init=epsg:4326'))
CRU_SPDF_CROP<-crop(CRU_SPDF, extent(StudyAreaAlberte + 2))
crs(CRU_SPDF_CROP)<-CRS('+init=epsg:4326')
SUN_REF<-as.data.frame(CRU_SPDF_CROP)

###Créer le raster de REF actuel du SUNSHINE###
MyRaster2<-raster(ext=extent(StudyAreaAlberte + 2), resolution=0.2)
REF_S<-rasterize(SUN_REF[,13:14], MyRaster2, SUN_REF[1:12] )
REF_S<-resample(REF_S, RasterAlberte)
REF_S<-mask(REF_S, REF_T)
crs(REF_S)<-CRS("+init=epsg:4326")

###Ecrire le raster de REF du SUNSHINE T###
writeRaster(REF_S, "C:/Users/Bernigaud Nicolas/Documents/Article_SMA/LPJmL/1300_2500BP_V2/REF/CRU_SUN_REF.asc", bylayer=T)

###Resampler les raster d'anomalies d'après le modéle de REF###
for (i in 1:13){
  assign(paste("ANO_S_",i, "_RS", sep = ""), resample(get(paste("ANO_S_", i, sep="")), REF_S))
}

###Additionner les anomalies de SUNSHINE au modèle actuel###
for (i in 1:13){
  assign((paste("S_", 12 + i,"00BP", sep="")), get(paste("ANO_S_",i, "_RS", sep = "")) + REF_S)
}

###Transformer les fichiers SUNSHINE en CLOUDINESS###
for (i in 1:13){
  assign((paste("C_", 12 + i,"00BP", sep="")), 100 - get(paste("S_", 12 + i,"00BP", sep="")))
}

###Ecrire les rasters de CLOUDINESS###
for (i in 1:13){
  writeRaster(get(paste("C_", 12 + i, "00BP", sep="")), filename = paste0("Cloud_",12 + i, "00BP.asc"), bylayer=T, overwrite=TRUE)
}

###NOMBRE_DE_JOURS_PLUIES###
NORRR1MM<-matrix(nrow = 8602,ncol=14)
colnames(NORRR1MM)<-c("NJP1", "NJP2", "NJP3", "NJP4", "NJP5", "NJP6", "NJP7", "NJP8", "NJP9", "NJP10", "NJP11", "NJP12", "LAT", "LONG")
seq1<-seq(1,103224, by=12)
for (i in 1:12){
  NORRR1MM[,i]<-Precip_REF_mensuel[seq1+i-1,11]
}
NORRR1MM[,13]<-Precip_REF_mensuel[seq1,2]
NORRR1MM[,14]<-Precip_REF_mensuel[seq1,3]
NORRR1MM_DF<-as.data.frame(NORRR1MM)

###Créer le raster de REF actuel du nombre de jours de pluies###
REF_NJP<-rasterize(NORRR1MM_DF[,14:13], MyRaster, NORRR1MM_DF[1:12] )
REF_NJP<-crop(REF_NJP, StudyAreaAlberte)
REF_NJP<-resample(REF_NJP, RasterAlberte)

###Ecrire le raster de réf###
writeRaster(REF_NJP, filename = "REF_NJDP.asc", bylayer=T, overwrite=TRUE)
