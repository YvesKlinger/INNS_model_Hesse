library(terra)
library(sf)
library(dplyr)
library(whitebox)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Predictors/, Predictors_raw/, Shape_Hesse/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  DEM 1000 m (BKG 2021), Predictors_raw/DGM/; heat index shapefile
#         (Wan et al. 2015 / HLNUG), Predictors_raw/MITHBEL20J/;
#         Shape_Hesse/Hessen.shp; Predictors/CHELSA/bio1.tif (03_CHELSA.R, target grid)
# Output: Predictors/{DGM,aspect,twi,heat_index}_hesse.tif
#         (intermediate files in Predictors_raw/processing/)

for (d in c("DGM", "aspect", "twi", "slope")) dir.create(file.path("Predictors_raw/processing", d), recursive = TRUE, showWarnings = FALSE)

DGM_ori <- rast("Predictors_raw/DGM/dgm_original_data/dgm1000/dgm1000_utm32s.asc")
hessen_shp <- vect("Shape_Hesse/Hessen.shp")
DGM <- crop(DGM_ori, hessen_shp)

writeRaster(DGM,"Predictors_raw/processing/DGM/dgm1000_cropped.tif",overwrite = TRUE)

#aspect
#####################
wbt_aspect(dem = "Predictors_raw/processing/DGM/dgm1000_cropped.tif", output = "Predictors_raw/processing/aspect/aspect.tif")

aspect <- rast("Predictors_raw/processing/aspect/aspect.tif")

#slope & twi
##################################
flow_acc <- wbt_d_inf_flow_accumulation(input = "Predictors_raw/processing/DGM/dgm1000_cropped.tif",
                            output = "Predictors_raw/processing/twi/flow_acc.tif",
                            out_type = "Specific Contributing Area",
                            log = FALSE)

wbt_slope(dem = "Predictors_raw/processing/DGM/dgm1000_cropped.tif",
          output = "Predictors_raw/processing/slope/slope.tif",
          units = "degrees")

wbt_wetness_index(sca = "Predictors_raw/processing/twi/flow_acc.tif",
                  slope = "Predictors_raw/processing/slope/slope.tif",
                  output = "Predictors_raw/processing/twi/TWI.tif")

twi <- rast("Predictors_raw/processing/twi/TWI.tif")

#wgs & hessen shape
######################
value_list <- list(DGM = DGM, aspect = aspect, twi = twi)
chelsa <- rast("Predictors/CHELSA/bio1.tif")
hessen_wgs <- project(hessen_shp,"EPSG:4326")
for (value_name in names(value_list)) {

  value <- value_list[[value_name]]

  value_wgs <- project(value, "EPSG:4326")
  value_resampled <- resample(value_wgs, chelsa, method = "bilinear")
  value_cropped <- crop(value_resampled, hessen_wgs)
  value_hesse <- mask(value_cropped, hessen_wgs)

  writeRaster(value_hesse, paste0("Predictors/", value_name, "_hesse.tif"), overwrite = TRUE)
}


#Hitzeindex
#################################
dgm <- rast("Predictors/DGM_hesse.tif")

HI_shp <- st_read("Predictors_raw/MITHBEL20J/MITHBEL20J.shp") %>%
  st_transform("EPSG:4326")

HI_rasterized <- rasterize(HI_shp, dgm, field = "KLASSE", touches = TRUE)
HI_cropped <- crop(HI_rasterized, dgm)
HI_hesse <- mask(HI_cropped, dgm)
plot(HI_hesse)
levels(HI_hesse)

writeRaster(HI_hesse, "Predictors/heat_index_hesse.tif", overwrite = TRUE)
