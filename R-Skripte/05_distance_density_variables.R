library(terra)
library(dplyr)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Predictors/, Predictors_raw/, Shape_Hesse/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  ATKIS Basis-DLM Hesse (HVBG 2025), German road/railroad/settlement layers
#         (BKG 2024), water bodies (BfG 2024; HLNUG 2025a), Predictors_raw/;
#         Shape_Hesse/Hessen.shp; Predictors/DGM_hesse.tif (04_DGM_und_abhaengige_Variablen.R)
# Output: Predictors/{roads,railroads,flowing_water,standing_water}_distance_hesse.tif,
#         Predictors/railroads_density_hesse.tif, Predictors/urban_builup_density_hesse.tif

base_shp <- vect("Shape_Hesse/Hessen.shp") %>% project("EPSG:4326")
base_raster <- rast("Predictors/DGM_hesse.tif")

road_hess <- vect("Predictors_raw/land_use/ATKIS/ATKIS Basis DLM Hessen (shape)/dlm_42002_ax_strasse.shp") %>% project("EPSG:4326")
road_ger <- vect("Predictors_raw/Distance_density_variables/roads/ver01_l_buff.shp") %>% project("EPSG:4326")
road_diff <- erase(road_ger, base_shp)
road_full <- rbind(road_diff, road_hess)

railroad_hess <- vect("Predictors_raw/land_use/ATKIS/ATKIS Basis DLM Hessen (shape)/dlm_42014_ax_bahnstrecke.shp") %>% project("EPSG:4326")
railroad_ger <- vect("Predictors_raw/Distance_density_variables/railroads/ver03_l_buff.shp") %>% project("EPSG:4326")
railroad_diff <- erase(railroad_ger, base_shp)
railroad_full <- rbind(railroad_diff, railroad_hess)

urban_buildup_hess <- vect("Predictors_raw/Distance_density_variables/urban_buildup/Objektartengruppe_Siedlung_merge_diss.shp") %>% project("EPSG:4326")

flowing_water_hess <- vect("Predictors_raw/Distance_density_variables/flowing_water/Gewaessernetz_DLM25__epsg25832/Gewaessernetz_DLM25.shp") %>% project("EPSG:4326")
flowing_water_ger <- vect("Predictors_raw/Distance_density_variables/flowing_water/riverWaterBody_DE_buff.shp") %>% project("EPSG:4326")
flowing_water_diff <- erase(flowing_water_ger, base_shp)
flowing_water_full <- rbind(flowing_water_diff, flowing_water_hess)

standing_water_hess <- vect("Predictors_raw/land_use/ATKIS/ATKIS Basis DLM Hessen (shape)/dlm_44006_ax_stehndgewssr.shp") %>% project("EPSG:4326")
standing_water_ger <- vect("Predictors_raw/Distance_density_variables/standing_water/lakeWaterBody_DE_buff.shp") %>% project("EPSG:4326")
standing_water_diff <- erase(standing_water_ger, base_shp)
standing_water_full <- rbind(standing_water_diff, standing_water_hess)
plot(standing_water_full)

#distance to
########################

list_objects <- list(roads = road_full, railroads = railroad_full, flowing_water = flowing_water_full, standing_water = standing_water_full)

for (object_name in names(list_objects)) {

object <- list_objects[[object_name]]
dist_raster <- distance(base_raster, object)
dist_raster_hesse <- mask(dist_raster, base_raster)

writeRaster(dist_raster_hesse, paste0("Predictors/",object_name,"_distance_hesse.tif"),overwrite=TRUE)
}

#density
##########################

railroad_lines <- as.lines(railroad_hess)
railroad_length <- rasterizeGeom(railroad_lines, base_raster, fun="length", "m")
railroad_length_hesse <- mask(railroad_length, base_raster)
writeRaster(railroad_length_hesse, "Predictors/railroads_density_hesse.tif", overwrite = TRUE)

urban_buildup_dens <-rasterize(urban_buildup_hess, base_raster,update = TRUE, fun='sum',cover=TRUE, background = 0)
urban_buildup_dens_hesse <- mask(urban_buildup_dens, base_raster)
writeRaster(urban_buildup_dens_hesse, "Predictors/urban_builup_density_hesse.tif", overwrite = TRUE)
