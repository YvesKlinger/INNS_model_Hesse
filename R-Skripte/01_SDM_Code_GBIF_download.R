library(rgbif)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Occurrences/, Predictors/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

species_list <- c("Ailanthus altissima", "Heracleum mantegazzianum", "Impatiens glandulifera")

species_keys <- sapply(species_list, function(sp) name_backbone(sp)$usageKey)

download_key_global <-occ_download(
  pred_in("taxonKey", species_keys),
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_in("basisOfRecord",c("OBSERVATION","MACHINE_OBSERVATION", "HUMAN_OBSERVATION","OCCURRENCE")),
  pred_gte("year", 1990),
  pred_or(  
    pred_lt("coordinateUncertaintyInMeters",1000),
    pred_isnull("coordinateUncertaintyInMeters")),
  format = "DWCA"
)

occ_download_meta(download_key_global)

download_key_local <-occ_download(
  pred_in("taxonKey", species_keys),
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_in("basisOfRecord",c("OBSERVATION","MACHINE_OBSERVATION", "HUMAN_OBSERVATION","OCCURRENCE")),
  pred("gadm","DEU.7_1"),
  pred_gte("year", 1990),
  pred_or(  
    pred_lt("coordinateUncertaintyInMeters",1000),
    pred_isnull("coordinateUncertaintyInMeters")),
  format = "DWCA"
)

occ_download_meta(download_key_local)
