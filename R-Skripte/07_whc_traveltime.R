library(terra)
library(dplyr)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Predictors/, Predictors_raw/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  usable field capacity 1:50,000 (HLNUG 2025b) and travel time to cities
#         2015 (Weiss et al. 2018), Predictors_raw/;
#         Predictors/DGM_hesse.tif (04_DGM_und_abhaengige_Variablen.R)
# Output: Predictors/whc_hesse.tif, Predictors/travel_time_hesse.tif

base_raster <- rast("Predictors/DGM_hesse.tif")

#whc
##############################
whc <- vect("Predictors_raw/water_holding_capacity/Nutzbare_Feldkapazitaet100cm_50000__epsg25832/Nutzbare_Feldkapazitaet100cm_50000.shp") %>% project("EPSG:4326")
whc$BEZEICHNER[whc$BEZEICHNER == "Abbau- und Aufschüttungsflächen"] <- "Siedlung und Verkehr"
whc$BEZEICHNER[whc$BEZEICHNER == "Flächen für Siedlung, Industrie und Verkehr"] <- "Siedlung und Verkehr"
whc$BEZEICHNER[whc$BEZEICHNER == "Siedlungs-Kernflächen"] <- "Siedlung und Verkehr"

whc_ras <- rasterize(whc, base_raster, field = "BEZEICHNER", touches = TRUE)

whc_cropped <- crop(whc_ras, base_raster)
whc_hesse <- mask(whc_cropped, base_raster)

levels(whc_hesse)
writeRaster(whc_hesse, "Predictors/whc_hesse.tif", overwrite = TRUE)

#travel time
##############################
travel_time <- rast("Predictors_raw/Travel_time/201501_Global_Travel_Time_to_Cities_2015.tif") %>% crop(base_raster)
travel_time_hesse <- mask(travel_time, base_raster)
plot(travel_time_hesse)

writeRaster(travel_time_hesse, "Predictors/travel_time_hesse.tif", overwrite = TRUE)
