library(dplyr)
library(data.table)
library(sf)
library(terra)
library(CoordinateCleaner)
library(ggplot2)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Occurrences/, Predictors/, Shape_Hesse/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Number of cleaned but unthinned occurrence records per 1-km cell in Hesse per
# species and data source (CS = GBIF, StAg = HLNUG). The records are cleaned exactly as
# in 02_SDM_Code_Occurrence_filter.R (year, coordinate uncertainty,
# CoordinateCleaner); only the spatial thinning to one record per 1-km cell is
# skipped.
#
# Input:  Occurrences/Occurrence_data_GBIF/local/0021801-250717081556266/occurrence.txt (GBIF download, Hesse)
#         Occurrences/Occurrence_data_HLNUG/unfiltered/<name>/<name>.shp (HLNUG, on request)
#         Shape_Hesse/Hessen.shp, Predictors/DGM_hesse.tif (1-km model grid)
# Output: Model_results/Figures/Occurrence_maps_cells.{png,pdf}   (records per 1-km cell; Suppl. App. 1.2)
#         Model_results/Figures/Occurrence_maps_counts.csv        (records and occupied cells per panel)

out_dir <- "Model_results/Figures/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_names <- c("3190653" = "A. altissima",
                   "3034825" = "H. mantegazzianum",
                   "2891770" = "I. glandulifera")
source_names <- c(gbif = "Citizen science (CS)", hlnug = "State agency (StAg)")


#load data (as in 02_SDM_Code_Occurrence_filter.R)
##########################################################################
##########################################################################

occ_table_local <- fread("Occurrences/Occurrence_data_GBIF/local/0021801-250717081556266/occurrence.txt",
                         select = c("gbifID", "license", "institutionCode", "basisOfRecord", "year", "month", "day",
                                    "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
                                    "genus", "specificEpithet", "speciesKey")) %>%  filter(!is.na(speciesKey))  %>% filter(!speciesKey %in% c(3170247, 2869311))
occ_table_local$speciesKey <- as.character(occ_table_local$speciesKey)

shp_list <- lapply(c("Chinesischer_Götterbaum", "Drüsiges_Springkraut", "Riesen_Bärenklau"),
                   function(n) st_read(paste0("Occurrences/Occurrence_data_HLNUG/unfiltered/", n, "/", n, ".shp"), quiet = TRUE))
shp_all_wgs <- st_transform(do.call(rbind, shp_list), 4326)
coords <- st_coordinates(shp_all_wgs) %>% as.data.frame()
shp_all_wgs$decimalLongitude <- coords$X
shp_all_wgs$decimalLatitude <- coords$Y
occ_table_hlnug <- shp_all_wgs %>%
  rename(speciesKey = ART_WISS) %>%
  filter(speciesKey != "Lysichiton americanus") %>%
  mutate(speciesKey = recode(speciesKey,
                             "Ailanthus altissima" = 3190653,
                             "Heracleum mantegazzianum" = 3034825,
                             "Impatiens glandulifera" = 2891770)) %>%
  st_drop_geometry() %>%
  as.data.frame()
occ_table_hlnug$speciesKey <- as.character(occ_table_hlnug$speciesKey)

occ_table_ <- list(gbif = occ_table_local, hlnug = occ_table_hlnug)


#cleaning without thinning (filters as in 02_SDM_Code_Occurrence_filter.R)
##########################################################################
##########################################################################

clean_list <- list()

for (occ_name in names(occ_table_)) {
  for (species in names(species_names)) {
    occ <- occ_table_[[occ_name]] %>% filter(speciesKey == species)

    #filter year and uncertainty (GBIF: already applied in the download query)
    if (occ_name == "hlnug") {
      occ <- occ %>%
        filter(is.na(JAHR) | JAHR >= 1990) %>%
        filter(is.na(TOLERANZ) | !(TOLERANZ %in% c("Raster Quadrant", "ungenau")))
    }

    occ_val <- cc_val(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "clean", verbose = TRUE )
    occ_clean <- clean_coordinates(x = occ_val,
                                   lon = "decimalLongitude",
                                   lat = "decimalLatitude",
                                   species = "speciesKey",
                                   value = "clean",
                                   tests = c("centroids", "duplicates", "equal", "gbif", "outliers", "institutions", "zeros", "seas"),
                                   verbose = TRUE
    )

    clean_list[[paste(occ_name, species)]] <- data.frame(
      source = occ_name, species = species,
      decimalLongitude = occ_clean$decimalLongitude, decimalLatitude = occ_clean$decimalLatitude)
  }
}

occ_clean_all <- bind_rows(clean_list) %>%
  mutate(source  = factor(source, levels = names(source_names), labels = source_names),
         species = factor(species, levels = names(species_names), labels = species_names))


#map basics
##########################################################################
##########################################################################

map_crs <- 25832  # ETRS89 / UTM 32N for display
hesse <- st_read("Shape_Hesse/Hessen.shp", quiet = TRUE) %>% st_transform(map_crs)
grid_1km <- rast("Predictors/DGM_hesse.tif")

# number of records and of occupied 1-km cells per panel
occ_clean_all$cell <- cellFromXY(grid_1km, as.matrix(occ_clean_all[, c("decimalLongitude", "decimalLatitude")]))
counts <- occ_clean_all %>%
  group_by(source, species) %>%
  summarise(n_records = n(), n_cells = n_distinct(cell, na.rm = TRUE), .groups = "drop")
write.csv(counts, paste0(out_dir, "Occurrence_maps_counts.csv"), row.names = FALSE)
print(counts)

panel_labels <- counts %>%
  mutate(label = paste0("n = ", format(n_records, big.mark = ",", trim = TRUE)))

map_theme <- theme_void(base_size = 10) +
  theme(
    strip.text.x = element_text(face = "italic", size = 10, margin = margin(b = 4)),
    strip.text.y = element_text(size = 10, angle = 270, margin = margin(l = 4)),
    panel.spacing = unit(6, "pt"),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", colour = NA)
  )

label_layer <- geom_text(data = panel_labels, aes(label = label), x = -Inf, y = Inf,
                         hjust = -0.1, vjust = 1.5, size = 3, colour = "grey20", inherit.aes = FALSE)


#Records per 1-km cell
##########################################################################
##########################################################################

cell_list <- list()
for (s in levels(occ_clean_all$source)) {
  for (sp in levels(occ_clean_all$species)) {
    d <- occ_clean_all %>% filter(source == s, species == sp)
    r <- rasterize(vect(d, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326"),
                   grid_1km, fun = "count")
    names(r) <- "n"
    cells <- st_as_sf(as.polygons(r, dissolve = FALSE, na.rm = TRUE))
    cells$source <- s
    cells$species <- sp
    cell_list[[paste(s, sp)]] <- cells
  }
}
cells_all <- bind_rows(cell_list) %>%
  mutate(source  = factor(source, levels = source_names),
         species = factor(species, levels = species_names))

p_cells <- ggplot() +
  geom_sf(data = hesse, fill = "white", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = cells_all, aes(fill = n), colour = NA) +
  scale_fill_gradient(low = "#6BAED6", high = "#08306B", trans = "log10",
                      breaks = c(1, 3, 10, 30, 100),
                      name = "Records per 1-km cell") +
  label_layer +
  facet_grid(source ~ species) +
  coord_sf(crs = map_crs, datum = NA) +
  map_theme +
  guides(fill = guide_colourbar(barwidth = unit(4, "cm"), barheight = unit(0.3, "cm"), title.vjust = 0.8))

ggsave(paste0(out_dir, "Occurrence_maps_cells.png"), p_cells, width = 18, height = 14, units = "cm", dpi = 300, bg = "white")
ggsave(paste0(out_dir, "Occurrence_maps_cells.pdf"), p_cells, width = 18, height = 14, units = "cm")
