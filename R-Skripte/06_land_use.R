library(terra)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Predictors/, Predictors_raw/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  CORINE Land Cover 2018 (EEA 2019), cropped to Hesse in ArcGIS
#         (all attribute columns except "value" deleted), Predictors_raw/land_use/CORINE/;
#         Predictors/DGM_hesse.tif (04_DGM_und_abhaengige_Variablen.R)
# Output: Predictors/CORINE_hesse.tif

base_raster <- rast("Predictors/DGM_hesse.tif")

#CORINE
###############################
#delete all columns exept value in ArcGIS and crop to Hesse
corine_cropped <- rast("Predictors_raw/land_use/CORINE/corine_crop.tif")
reclass_matrix <- matrix(c(
  1, 11,
  2, 11,
  3, 12,
  4, 12,
  5, 12,
  6, 12,
  7, 13,
  8, 13,
  9, 13,
  10, 14,
  11, 14,
  12, 21,
  15, 22,
  16, 22,
  18, 23,
  20, 24,
  21, 24,
  23, 31,
  24, 31,
  25, 31,
  26, 32,
  27, 32,
  29, 32,
  31, 33,
  35, 41,
  36, 41,
  40, 51,
  41, 51
), ncol = 2, byrow = TRUE)

corine_reclassified <- classify(corine_cropped, rcl = reclass_matrix)

corine_wgs <- project (corine_reclassified,"EPSG:4326", method = "near")
corine_crop2 <- crop(corine_wgs, ext(base_raster))
corine_crop_resampled <- resample(corine_crop2, base_raster, method = "near")
corine_hesse <- mask(corine_crop_resampled, base_raster)
plot(corine_hesse)
writeRaster(corine_hesse, "Predictors/CORINE_hesse.tif", overwrite = TRUE)
