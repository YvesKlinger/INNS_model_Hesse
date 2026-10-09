library(raster)
library(terra)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Predictors/, Predictors_raw/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  raw CHELSA v2.1 bioclim layers (1981-2010), Predictors_raw/CHELSA/
# Output: Predictors/CHELSA/bio<1-19>.tif (scaled to degC / mm), read by
#         08_correlation_matrix_complete_model_area.R (all 19) and
#         10_SDM_Code_Model_final.R (bio1, 4, 15, 16, 17, 18)

dir.create("Predictors/CHELSA", recursive = TRUE, showWarnings = FALSE)

results_list <- list()
bioclim_list <- 1:19

for (bioclim in bioclim_list) {

  bioclim_path <- paste0("Predictors_raw/CHELSA/CHELSA_bio",bioclim,"_1981-2010_V.2.1.tif")
  raster_layer <- raster(bioclim_path)
  names(raster_layer) <- paste0("bio", bioclim)
  results_list[[bioclim]] <- raster_layer
}

chelsa_bioclim <- stack(results_list)

#scale rasters
#################
bioclim_scale_list_1 <- list(1,5,6,8,9,10,11)

for (bioclim in bioclim_scale_list_1) {

  bioclim_1 <- chelsa_bioclim[[bioclim]]
  bioclim_1_scaled <- bioclim_1 * 0.1 -273.15
  names(bioclim_1_scaled) <- paste0("bio", bioclim,"_scaled")
  writeRaster(bioclim_1_scaled, paste0("Predictors/CHELSA/bio",bioclim,".tif"), overwrite = FALSE)
}

bioclim_scale_list_2 <- list(2,3,4,7,12,13,14,15,16,17,18,19)

for (bioclim in bioclim_scale_list_2) {

  bioclim_2 <- chelsa_bioclim[[bioclim]]
  bioclim_2_scaled <- bioclim_2 * 0.1
  names(bioclim_2_scaled) <- paste0("bio", bioclim,"_scaled")
  writeRaster(bioclim_2_scaled, paste0("Predictors/CHELSA/bio",bioclim,".tif"), overwrite = FALSE)
}
