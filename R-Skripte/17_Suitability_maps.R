# ------------------------ Libraries ------------------------
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(stringr)
library(stars)
library(cowplot)   # plot_grid(), get_legend(), ggdraw()
library(grid)      # unit()

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Model_results/, Shape_Hesse/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  Model_results/Maps/inside_aoa/pred_ras_mean_aoa_local_<origin>_<species>_mean.tif
#         (mean suitability masked by each model's own AOA; 10_SDM_Code_Model_final.R),
#         Shape_Hesse/Hessen.shp
# Output: Model_results/Maps/Map_Plots/final_suitability_plot.png (Figure 5),
#         suitability and difference rasters in Model_results/Maps/Map_Plots/Maps_for_Plots/

# ------------------------ Inputs ------------------------
hesse <- vect("Shape_Hesse/Hessen.shp")
# Keep a single sf copy; don't overwrite inside the loop
hesse_sf <- st_as_sf(hesse)

path <- "Model_results/Maps/inside_aoa"

species_list <- c("2891770", "3190653", "3034825")
species_names <- c(
  "2891770" = "I. glandulifera",
  "3190653" = "A. altissima",
  "3034825" = "H. mantegazzianum"
)

# ------------------------ Files & metadata ------------------------
map_files <- list.files(path, pattern = "\\.tif$", full.names = TRUE)

parse_meta <- function(fp){
  fn <- basename(fp)
  origin <- str_extract(fn, "local_[^_]+")
  sp <- str_extract(fn, "\\d{7}")
  data.frame(file=fp, origin=origin, species=sp, stringsAsFactors=FALSE)
}
meta_df <- do.call(rbind, lapply(map_files, parse_meta))

# ------------------------ Theme (TIGHTER) ------------------------
map_theme <- theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.title     = element_text(size = 11),
    plot.tag         = element_text(face = "plain", size = 13),
    plot.tag.position = c(0, 1),          # top-left inside panel
    plot.margin      = margin(1,1,1,1)    # tighter margins
  )

# ------------------------ Storage ------------------------
list_plot_com   <- list()
list_plot_gbif  <- list()
list_plot_hlnug <- list()

# ------------------------ Loop over species ------------------------
for(sp_id in species_list){
  row <- meta_df %>% filter(species == sp_id)
  f_comb  <- row$file[row$origin=="local_combined"]
  f_gbif  <- row$file[row$origin=="local_gbif"]
  f_hlnug <- row$file[row$origin=="local_hlnug"]
  
  # --- Load rasters ---
  r_comb  <- rast(f_comb)
  r_gbif  <- rast(f_gbif)
  r_hlnug <- rast(f_hlnug)
  
  # --- Transform Hesse to target CRS (once per species) ---
  target_crs <- crs(r_comb)
  hesse_sf_sp <- st_transform(hesse_sf, st_crs(target_crs))
  
  # --- Align CRSs & grids to r_comb ---
  if (crs(r_gbif) != crs(r_comb))  r_gbif  <- terra::project(r_gbif,  r_comb, method = "bilinear")
  if (!terra::compareGeom(r_gbif, r_comb, stopOnError = FALSE))
    r_gbif <- terra::resample(r_gbif, r_comb, method = "bilinear")
  
  if (crs(r_hlnug) != crs(r_comb)) r_hlnug <- terra::project(r_hlnug, r_comb, method = "bilinear")
  if (!terra::compareGeom(r_hlnug, r_comb, stopOnError = FALSE))
    r_hlnug <- terra::resample(r_hlnug, r_comb, method = "bilinear")
  
  # --- Mask to Hesse (clean edges) ---
  r_comb  <- terra::mask(r_comb,  vect(hesse_sf_sp))
  r_gbif  <- terra::mask(r_gbif,  vect(hesse_sf_sp))
  r_hlnug <- terra::mask(r_hlnug, vect(hesse_sf_sp))
  
  # --- Differences ---
  diff_gbif  <- r_gbif  - r_comb
  diff_hlnug <- r_hlnug - r_comb
  
  # --- Convert to stars for consistent plotting ---
  comb_st       <- st_as_stars(r_comb)
  diff_gbif_st  <- st_as_stars(diff_gbif)
  diff_hlnug_st <- st_as_stars(diff_hlnug)
  
  # --- Symmetric limits for the diff color bar (per-species) ---
  max_abs <- max(
    abs(c(as.vector(terra::values(diff_gbif)),
          as.vector(terra::values(diff_hlnug)))),
    na.rm = TRUE
  )
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  # --- Plot: COM (light grey NA background inside Hesse) ---
  p_comb <- ggplot() +
    geom_sf(data = hesse_sf_sp, fill = "lightgrey", color = NA) +
    geom_stars(data = comb_st, na.rm = FALSE) +
    scale_fill_viridis_c(
      name = "Suitability",
      option = "D",
      na.value = NA,
      guide = guide_colorbar(
        direction = "vertical",
        barheight = unit(10, "lines"),
        barwidth  = unit(2,  "lines"),
        ticks = FALSE,
        label.position = "right",
        title.position = "top",
        title.hjust = 0,
        title.vjust = 2
      )
    ) +
    geom_sf(data = hesse_sf_sp, fill = NA, color = "grey40", linewidth = 0.3) +
    coord_sf(crs = st_crs(comb_st), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    map_theme
  list_plot_com[[sp_id]] <- p_comb
  
  # --- Shared diff scale (VERTICAL, bigger, same size as COM legend) ---
  diff_scale <- scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0, limits = c(-max_abs, max_abs),
    oob = scales::squish, na.value = NA,
    name = NULL,
    breaks = c(-max_abs, max_abs),
    labels = c("less suitable", "more suitable"),
    guide = guide_colorbar(
      direction = "vertical",
      barheight = unit(10, "lines"),
      barwidth  = unit(2,  "lines"),
      ticks = FALSE,
      label.position = "right"
    )
  )
  
  # --- Plot: CS (GBIF vs COM) with light grey NA interior ---
  p_gbif <- ggplot() +
    geom_sf(data = hesse_sf_sp, fill = "lightgrey", color = NA) +
    geom_stars(data = diff_gbif_st, na.rm = FALSE) +
    diff_scale +
    geom_sf(data = hesse_sf_sp, fill = NA, color = "grey40", linewidth = 0.3) +
    coord_sf(crs = st_crs(comb_st), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    map_theme +
    theme(plot.title = element_blank())
  list_plot_gbif[[sp_id]] <- p_gbif
  
  # --- Plot: StAg (HLNUG vs COM) with light grey NA interior ---
  p_hlnug <- ggplot() +
    geom_sf(data = hesse_sf_sp, fill = "lightgrey", color = NA) +
    geom_stars(data = diff_hlnug_st, na.rm = FALSE) +
    diff_scale +
    geom_sf(data = hesse_sf_sp, fill = NA, color = "grey40", linewidth = 0.3) +
    coord_sf(crs = st_crs(comb_st), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    map_theme +
    theme(plot.title = element_blank())
  list_plot_hlnug[[sp_id]] <- p_hlnug
  
  # ------------------------ Save rasters for ArcGIS ------------------------
  ras_out_dir <- "Model_results/Maps/Map_Plots/Maps_for_Plots"
  if (!dir.exists(ras_out_dir)) dir.create(ras_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # helper to sanitize file names
  safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  
  sp_label <- species_names[[sp_id]]
  
  # filenames (suitability rasters by origin)
  f_comb_tif   <- file.path(ras_out_dir, paste0(safe_name(sp_label), "_COM.tif"))
  f_gbif_tif   <- file.path(ras_out_dir, paste0(safe_name(sp_label), "_CS.tif"))
  f_hlnug_tif  <- file.path(ras_out_dir, paste0(safe_name(sp_label), "_StAg.tif"))
  
  # filenames (difference rasters vs COM)
  f_diff_cs_tif   <- file.path(ras_out_dir, paste0(safe_name(sp_label), "_CS_vs_COM.tif"))
  f_diff_stag_tif <- file.path(ras_out_dir, paste0(safe_name(sp_label), "_StAg_vs_COM.tif"))
  
  # write options for ArcGIS-friendly GeoTIFF
  gtiff_opts <- c("COMPRESS=LZW", "TILED=YES")  # small files, better I/O
  na_flag <- -9999                               # explicit NoData value (ArcGIS reads this well)
  
  # Save suitability rasters (continuous float)
  terra::writeRaster(r_comb,  f_comb_tif,
                     filetype = "GTiff", gdal = gtiff_opts,
                     datatype = "FLT4S", NAflag = na_flag, overwrite = TRUE)
  
  terra::writeRaster(r_gbif,  f_gbif_tif,
                     filetype = "GTiff", gdal = gtiff_opts,
                     datatype = "FLT4S", NAflag = na_flag, overwrite = TRUE)
  
  terra::writeRaster(r_hlnug, f_hlnug_tif,
                     filetype = "GTiff", gdal = gtiff_opts,
                     datatype = "FLT4S", NAflag = na_flag, overwrite = TRUE)
  
  # Save difference rasters (GBIF–COM, StAg–COM)
  terra::writeRaster(diff_gbif,  f_diff_cs_tif,
                     filetype = "GTiff", gdal = gtiff_opts,
                     datatype = "FLT4S", NAflag = na_flag, overwrite = TRUE)
  
  terra::writeRaster(diff_hlnug, f_diff_stag_tif,
                     filetype = "GTiff", gdal = gtiff_opts,
                     datatype = "FLT4S", NAflag = na_flag, overwrite = TRUE)
}

# ------------------------ Assemble Panels ------------------------

safe_get_plot <- function(lst, key){
  if(is.null(lst) || !(key %in% names(lst)) || is.null(lst[[key]])){
    ggplot() + theme_void() + labs(title = paste0("missing: ", key))
  } else lst[[key]]
}

species_order <- c("3190653","3034825","2891770")

row_com   <- lapply(species_order, function(s) safe_get_plot(list_plot_com, s))
row_gbif  <- lapply(species_order, function(s) safe_get_plot(list_plot_gbif, s))
row_hlnug <- lapply(species_order, function(s) safe_get_plot(list_plot_hlnug, s))

# Remove legends from maps
row_com_no_legend   <- lapply(row_com,   function(p) p + theme(legend.position = "none"))
row_gbif_no_legend  <- lapply(row_gbif,  function(p) p + theme(legend.position = "none"))
row_hlnug_no_legend <- lapply(row_hlnug, function(p) p + theme(legend.position = "none"))

# ------------------------ Assign panel tags a) ... i) ------------------------
panel_tags <- paste0(letters[1:9], ")")
all_plots <- c(row_com_no_legend, row_gbif_no_legend, row_hlnug_no_legend)
for (i in seq_along(all_plots)) {
  all_plots[[i]] <- all_plots[[i]] +
    labs(tag = panel_tags[i]) +
    theme(plot.tag.position = c(0, 1))
}
row_com_no_legend   <- all_plots[1:3]
row_gbif_no_legend  <- all_plots[4:6]
row_hlnug_no_legend <- all_plots[7:9]

# ------------------------ Row labels ------------------------
label_plot <- function(txt){
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = txt, size = 6, fontface = "plain") +
    theme(plot.margin = margin(0,0,0,0))
}

# Build a row: left label + three species maps
build_row_grob <- function(label, maps_no_legend){
  maps_row <- plot_grid(plotlist = maps_no_legend, ncol = 3,
                        align = "hv", axis = "tblr")
  plot_grid(label_plot(label), maps_row, ncol = 2,
            rel_widths = c(0.06, 0.94),  # a tiny bit more left padding to avoid cut-off
            align = "h")
}

row1_grob <- build_row_grob("COM",  row_com_no_legend)
row2_grob <- build_row_grob("CS",   row_gbif_no_legend)
row3_grob <- build_row_grob("StAg", row_hlnug_no_legend)

# ---------- Strong vertical compression, but safe from clipping ----------
maps_grid <- plot_grid(
  row1_grob,
  row2_grob,
  row3_grob,
  ncol = 1,
  rel_heights = c(1, 0.85, 0.85),
  align = "v"
)
# keep tight spacing but avoid cut-off: moderate negative, not extreme
maps_grid <- ggdraw(maps_grid) + theme(plot.margin = margin(-12, 0, -12, 0, "pt"))

# ------------------------ Legends ------------------------
legend_theme <- theme_void() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.justification = c(0, 0.5),
    legend.text = element_text(size = 10.5),
    legend.title = element_text(size = 10.5, face = "plain"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0)
  )

df_dummy_c <- data.frame(x = 1:2, y = 1:2, z = c(0, 1))
p_dummy_c <- ggplot(df_dummy_c, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(
    name = "Suitability",
    option = "D",
    na.value = "white",
    guide = guide_colorbar(
      direction = "vertical",
      barheight = unit(8,  "lines"),
      barwidth  = unit(1.6,"lines"),
      ticks = FALSE,
      label.position = "right",
      title.position = "top",
      title.hjust = 0,
      title.vjust = 2
    )
  ) +
  legend_theme
com_legend_grob <- get_legend(p_dummy_c)

df_dummy_d <- data.frame(x = 1:2, y = 1:2, z = c(-1, 1))
p_diff_vertical <- ggplot(df_dummy_d, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0, limits = c(-1, 1),
    name = NULL,
    breaks = c(-1, 1),
    labels = c("less suitable", "more suitable"),
    guide = guide_colorbar(
      direction = "vertical",
      barheight = unit(8,  "lines"),
      barwidth  = unit(1.6,"lines"),
      ticks = FALSE,
      label.position = "right"
    )
  ) +
  legend_theme
diff_legend_grob <- get_legend(p_diff_vertical)

legend_stack <- plot_grid(
  ggdraw() + draw_grob(com_legend_grob),
  ggdraw() + draw_grob(diff_legend_grob),
  ncol = 1,
  rel_heights = c(1, 1.1),
  align = "v",
  axis = "l"
)

legend_col <- ggdraw() +
  draw_plot(
    legend_stack,
    x = 0,
    y = 0.28,
    width = 1,
    height = 0.72,
    hjust = 0,
    vjust = 0
  ) +
  theme(plot.margin = margin(0, 0, 0, 0))

# ------------------------ Italic species headers above columns ------------------------
header_plot <- function(txt){
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = txt,
             fontface = "italic", size = 5) +
    theme(plot.margin = margin(0,0,0,0))
}

header_row <- plot_grid(
  header_plot("A. altissima"),
  header_plot("H. mantegazzianum"),
  header_plot("I. glandulifera"),
  ncol = 3, align = "hv"
)

# add a left spacer (for the COM/CS/StAg label column)
left_spacer <- ggplot() + theme_void()
header_with_pad <- plot_grid(left_spacer, header_row, ncol = 2,
                             rel_widths = c(0.06, 0.94), align = "h")

# ------------------------ Final composition ------------------------
# Increase header-to-map spacing slightly (match origin label distance)
maps_with_header <- plot_grid(
  header_with_pad,
  maps_grid,
  ncol = 1,
  rel_heights = c(0.10, 0.90),   # was 0.07/0.93 → now a bit more breathing room
  align = "v"
)

final_all_core <- plot_grid(
  maps_with_header,
  legend_col,
  ncol = 2,
  rel_widths = c(1, 0.13),
  align = "h"
)

# Outer margin to prevent clipping anywhere
final_all <- ggdraw(final_all_core) +
  theme(plot.margin = margin(14, 16, 16, 16, "pt"))

plot(final_all)

# ------------------------ Display / Save ------------------------
out_path <- "Model_results/Maps/Map_Plots/final_suitability_plot.png"

ggsave(
  filename = out_path,
  plot = final_all,
  width = 14,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)
