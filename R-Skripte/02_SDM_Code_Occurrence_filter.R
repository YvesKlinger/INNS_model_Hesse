#Setup
#################################################################
library(dismo)
library(raster)
library(sf)
library(sp)
library(dplyr)
library(CoordinateCleaner)
library(data.table)
library(terra)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Occurrences/, Predictors/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

#load needed data & prepare
#################################################################

#load predictors needed
results_list <- list()
bioclim_list <- c(1, 4, 15, 16, 17, 18)

for (bioclim in bioclim_list) {
  bioclim_path <- paste0("Predictors/CHELSA/bio", bioclim, ".tif")  
  raster_layer <- rast(bioclim_path)
  names(raster_layer) <- paste0("bio", bioclim)
  results_list[[paste0("bio", bioclim)]] <- raster_layer
}

predictors_global <- rast(results_list)

#extent
raster_hesse <- raster("Predictors/DGM_hesse.tif")
raster_global <- predictors_global[[1]]

###load gbif data

occ_raster_ <- list(local = raster_hesse, global = raster_global, hlnug =raster_hesse, combined = raster_hesse)

occ_table_global <- fread("Occurrences/Occurrence_data_GBIF/global/0021796-250717081556266/occurrence.txt", 
                          select = c("gbifID", "license", "institutionCode", "basisOfRecord", "year", "month", "day", 
                                     "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", 
                                     "genus", "specificEpithet", "speciesKey")) %>%  filter(!is.na(speciesKey)) %>% filter(!speciesKey %in% c(3170247, 2869311))
occ_table_global$speciesKey <- as.character(occ_table_global$speciesKey)
table(occ_table_global$speciesKey)

occ_table_local <- fread("Occurrences/Occurrence_data_GBIF/local/0021801-250717081556266/occurrence.txt", 
                         select = c("gbifID", "license", "institutionCode", "basisOfRecord", "year", "month", "day", 
                                    "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", 
                                    "genus", "specificEpithet", "speciesKey")) %>%  filter(!is.na(speciesKey))  %>% filter(!speciesKey %in% c(3170247, 2869311))
occ_table_local$speciesKey <- as.character(occ_table_local$speciesKey)
table(occ_table_local$speciesKey)

#load hlnug data
shp_Ailanthus <- st_read("Occurrences/Occurrence_data_HLNUG/unfiltered/Chinesischer_Götterbaum/Chinesischer_Götterbaum.shp")
shp_Heracleum <- st_read("Occurrences/Occurrence_data_HLNUG/unfiltered/Drüsiges_Springkraut/Drüsiges_Springkraut.shp")
shp_Impatiens <- st_read("Occurrences/Occurrence_data_HLNUG/unfiltered/Riesen_Bärenklau/Riesen_Bärenklau.shp")
shp_list <- list(shp_Ailanthus, shp_Heracleum, shp_Impatiens)
shp_all <- do.call(rbind, shp_list)
shp_all_wgs <- st_transform(shp_all, CRS("+proj=longlat +datum=WGS84"))
coords <- st_coordinates(shp_all_wgs) %>% as.data.frame()
shp_all_wgs$decimalLongitude <- coords$X  
shp_all_wgs$decimalLatitude <- coords$Y 
table(shp_all_wgs$ID_ART)
occ_table_hlnug <- shp_all_wgs %>%
  rename(speciesKey = ART_WISS) %>%
  filter(speciesKey != "Lysichiton americanus") %>%
  mutate(speciesKey = recode(speciesKey,
                             "Ailanthus altissima" = 3190653,
                             "Heracleum mantegazzianum" = 3034825,
                             "Impatiens glandulifera" = 2891770)) %>%
  st_drop_geometry()

#combine both
hlnug_all <- as.data.frame(occ_table_hlnug)
hlnug_all$speciesKey <- as.character(hlnug_all$speciesKey)
occ_table_combined <- bind_rows(hlnug_all,occ_table_local); table(occ_table_combined$speciesKey)


#create lists
occ_table_ <- list(local = occ_table_local, global = occ_table_global, hlnug = occ_table_hlnug, combined = occ_table_combined)

list_occ_type <- list("global","local","hlnug","combined")


#Filtering loop
#################################################################

for (occ_name in list_occ_type) {
  results_list <- list()
  species_list <- unique(occ_table_[[occ_name]]$speciesKey); species_list
  ras_results_list <- list()
  all_data_list <- list()
  flag_list <- list()
  
  for (species in species_list) {
    occ <- occ_table_[[occ_name]] %>% filter(speciesKey == species)
    
    #filter year and uncertainty
    if (occ_name %in% c("combined", "hlnug")) {
      occ <- occ %>% 
        filter(is.na(JAHR) | JAHR >= 1990) %>% 
        filter(is.na(TOLERANZ) | !(TOLERANZ %in% c("Raster Quadrant", "ungenau"))) 
    }
    
    occ_val <- cc_val(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "clean", verbose = TRUE )
    #filter standard occ problems with coordinate cleaner
    occ_clean <- clean_coordinates(x = occ_val,
                                   lon = "decimalLongitude",
                                   lat = "decimalLatitude",
                                   species = "speciesKey",
                                   value = "clean",
                                   tests = c("centroids", "duplicates", "equal", "gbif", "outliers", "institutions", "zeros", "seas"),
                                   verbose = TRUE
    )
    
    occ_coordinates <- occ_clean %>%
      dplyr::select(decimalLongitude, decimalLatitude)
    
    occ_coordinates_sp <- SpatialPoints(occ_coordinates, proj4string = CRS("+proj=longlat +datum=WGS84")) 
    
    if (occ_name == "global") {  
      hesse_extent <- ext(raster_hesse)
      hesse_extent_poly <- as.polygons(hesse_extent)
      occ_coordinates_vect <- vect(occ_coordinates_sp)
      occ_coordinates_sp <- erase(occ_coordinates_vect, hesse_extent_poly)
      occ_coordinates_sp <- as(occ_coordinates_sp, "Spatial")
    }
    
    occ_grid <- gridSample(occ_coordinates_sp, occ_raster_[[occ_name]], n = 1) %>% as.data.frame 
    
    if (occ_name == "global") {  
      occ_grid <- occ_grid %>%
        rename(decimalLongitude = x, decimalLatitude = y)
    }
    
    occ_reconnected <- occ_grid %>%
      inner_join(occ_clean, by = c("decimalLongitude", "decimalLatitude"))
    
    results_list[[species]] <- occ_reconnected
    
    if (occ_name %in% c("global", "local")) {
      file_path <- paste0("Occurrences/Occurrence_data_GBIF/",occ_name,"/filtered/", species, ".csv")
    }
    
    if (occ_name %in% c("combined", "hlnug")) {
      file_path <- paste0("Occurrences/Occurrence_data_",occ_name,"/filtered/", species, ".csv")
    }
    
    write.table(occ_reconnected, file_path, sep = ";", row.names = FALSE)
    
    species <- as.character(species)
    assign(species, read.table(file_path, sep = ";"))
    
    all_data_list[[species]] <- read.table(file_path, sep = ";", header = TRUE)
  }
  
  combined_species <- bind_rows(all_data_list)
  print(table(combined_species$speciesKey))
  
  if (occ_name %in% c("global", "local")) {
    write.table(combined_species, paste0("Occurrences/Occurrence_data_GBIF/",occ_name,"/filtered/all_species_gbif_",occ_name,".csv"), sep = ";", row.names = FALSE)
  }
  
  if (occ_name %in% c("combined", "hlnug")) {
    write.table(combined_species, paste0("Occurrences/Occurrence_data_",occ_name,"/filtered/all_species_",occ_name,".csv"), sep = ";", row.names = FALSE)
  }
  
}


