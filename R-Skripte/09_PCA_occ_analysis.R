#Setup
##########################################################################
##########################################################################
library(terra)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(ecospat)
library(patchwork)
library(ggrepel)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Occurrences/, Predictors/, PCA_results/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

#load data & prepare
##########################################################################
##########################################################################

# load local predictors
##########################################################################
main_directory <- file.path(getwd(),"Predictors/")
raster_files <- list.files(path = main_directory, pattern = "_hesse.*\\.tif$", recursive = TRUE, full.names = TRUE)
raster_stack <- rast(raster_files)
raster_names <- tools::file_path_sans_ext(basename(raster_files))
raster_names_clean <- gsub("_hesse$", "", raster_names)
names(raster_stack) <- raster_names_clean
predictors_local <- raster_stack[[!duplicated(names(raster_stack))]] %>% na.omit() 

df_predictors_local_num <- terra::extract(predictors_local, terra::cells(predictors_local)) %>%
  dplyr::select(any_of(c("aspect","DGM","flowing_water_distance",
                         "railroads_density","railroads_distance",
                         "roads_distance","standing_water_distance",
                         "travel_time","twi","urban_builup_density"))) 

#stack species-specific
predictors_local_2891770 <- predictors_local[[names(predictors_local) != "heat_index"]]
predictors_local_3034825 <- predictors_local[[names(predictors_local) != "heat_index"]]
predictors_local_3190653 <- predictors_local[[names(predictors_local) != "whc"]]

# load local occurrences
##########################################################################
occ_local_gbif <- read.csv2("Occurrences/Occurrence_data_GBIF/local/filtered/all_species_gbif_local.csv") %>% 
  mutate(presence = "present") %>% 
  mutate(weight = 1)

occ_local_hlnug <- read.csv2("Occurrences/Occurrence_data_HLNUG/filtered/all_species_hlnug.csv") %>% 
  mutate(presence = "present") %>% 
  mutate(weight = 1) 

loadings <- list()
p <- list()
pca <- list()
xlims <- list()
ylims <- list()
D_list<- list()
eq_test_list <- list()
species_list <- list("2891770","3190653","3034825")
species_names <- c("2891770" = "I. glandulifera",
                   "3190653" = "A. altissima",
                   "3034825" = "H. mantegazzianum")

for (species in species_list) {
  species_loop    <- species
  predictors_loop <- get(paste0("predictors_local_", species_loop))
  
  #Prepare & tag each origin, then filter by species
  occ_gbif  <- occ_local_gbif  %>% 
    filter(speciesKey == species_loop) %>%
    mutate(origin = "GBIF",
           decimalLongitude = as.numeric(decimalLongitude),
           decimalLatitude  = as.numeric(decimalLatitude)
    )
  
  occ_hlnug <- occ_local_hlnug %>% 
    filter(speciesKey == species_loop) %>%
    mutate(origin = "HLNUG",decimalLongitude = as.numeric(decimalLongitude),
           decimalLatitude  = as.numeric(decimalLatitude)
    )
  
  #same sampling
  NoP <- min(nrow(occ_gbif), nrow(occ_hlnug))

  # bind all
  occ_all <- bind_rows(occ_gbif, occ_hlnug) %>% 
    dplyr::select(speciesKey, decimalLongitude, decimalLatitude, origin)
  
  #Extract predictors at all points
  occ_all[c("decimalLongitude","decimalLatitude")] <-
    lapply(occ_all[c("decimalLongitude","decimalLatitude")], as.numeric)
  
  occ_env <- terra::extract(
    predictors_loop,
    vect(occ_all, geom = c("decimalLongitude","decimalLatitude"), 
         crs  = crs(predictors_loop))
  ) %>%
    bind_cols(occ_all) %>%
    drop_na() 
  
  #Select numeric predictors 
  env_vars <- occ_env %>%
    dplyr::select(any_of(c("aspect","DGM","flowing_water_distance",
                    "railroads_density","railroads_distance",
                    "roads_distance","standing_water_distance",
                    "travel_time","twi","urban_builup_density"))) %>% 
    rename(
      Aspect = aspect,
      DEM = DGM,
      Distance_Flowing_Water = flowing_water_distance,
      Railroad_Density = railroads_density,
      Distance_Railroads = railroads_distance,
      Distance_Roads = roads_distance,
      Distance_Standing_Water = standing_water_distance,
      Travel_Time_to_Cities = travel_time,
      TWI = twi,
      Urban_Density = urban_builup_density
    )
  
  
  #Run PCA
  ##########################################################################
  ##########################################################################
  #GBIF & HLNUG together in one PCA space
  pca_all <- prcomp(env_vars, scale. = TRUE, center = TRUE)
   
  pca[[species]] <- pca_all
  
  var_explained <- (pca_all$sdev)^2 / sum((pca_all$sdev)^2)
  var_explained[1:2]
  var1_2 <- var_explained[1:2] * 100
  cumulative <- sum(var1_2); print(paste0(species, ":", cumulative))

  #plot PCA_all
  ###################################################
  # Build scores_df (with origin)
  scores_df <- as.data.frame(pca_all$x[,1:2]) %>%
    rename(PC1 = PC1, PC2 = PC2) %>%
    bind_cols(occ_env %>% dplyr::select(origin))

  #calculate centroid for each origin
  scores_df <- scores_df %>%
    group_by(origin) %>%
    mutate(
      centroid_PC1     = mean(PC1),
      centroid_PC2     = mean(PC2),
      dist_to_centroid = sqrt((PC1 - centroid_PC1)^2 +
                                (PC2 - centroid_PC2)^2)
    ) %>%
    ungroup()


  #Build loadings_df
  loadings_df <- as.data.frame(pca_all$rotation[,1:2]) %>%
    rownames_to_column("variable") %>%
    rename(PC1 = PC1, PC2 = PC2)
  
  loadings[[species]] <- loadings_df

  write.csv(loadings_df,paste0("PCA_results/PCA_loadings/PCA_loadings_",species,".csv"), row.names = FALSE)


  #Compute arrow multiplier
  arrow_mult <- min(
    (max(scores_df$PC1) - min(scores_df$PC1)) /
      (max(loadings_df$PC1) - min(loadings_df$PC1)),
    (max(scores_df$PC2) - min(scores_df$PC2)) /
      (max(loadings_df$PC2) - min(loadings_df$PC2))
  ) * 0.7

  # Compute mean positions for each origin
  means_df <- scores_df %>%
    group_by(origin) %>%
    summarise(mean_PC1 = mean(PC1),
              mean_PC2 = mean(PC2))

  #find overlapping occs and sign origin as "overlapping"
  scores_df_wOL <- scores_df %>%
    group_by(PC1, PC2) %>%
    mutate(
      n_gh = sum(origin %in% c("GBIF", "HLNUG"))
    ) %>%
    ungroup() %>%
    mutate(
      overlap_status = case_when(
        origin %in% c("GBIF", "HLNUG") & n_gh > 1 ~ "overlapping",
        TRUE                                     ~ as.character(origin)
      )
    )

  #filter outliers & build hulls
  filtered_df <- scores_df %>%
    group_by(origin) %>%
    filter(dist_to_centroid <= quantile(dist_to_centroid, 0.95)) %>%
    ungroup()

  hulls_df <- filtered_df %>%
    group_by(origin) %>%
    slice(chull(PC1, PC2)) %>%
    ungroup()
  
  # Plot
  ##########################################################################
  ##########################################################################
  #set colors
  color_origins <- c(
    "GBIF" = "darkgreen",
    "HLNUG" = "purple",
    "overlapping" = "brown"
  )

  #boundaries based on hulls
  xlims[[species]] <- range(hulls_df$PC1)
  ylims[[species]] <- range(hulls_df$PC2)
  
  species_table <- data.frame(
    letter  = c("(a)", "(b)", "(c)"),
    id      = c("3190653", "3034825", "2891770"),   # a) A. altissima, b) Heracleum, c) Impatiens
    species = c("A. altissima", "H. mantegazzianum", "I. glandulifera"),
    stringsAsFactors = FALSE
  )
  
  plot_letter <- species_table$letter[species_table$id == species]
  
  
  # Plot
  p[[species]] <- ggplot()  +
    geom_point(
      data = scores_df_wOL,
      aes(x = PC1, y = PC2, color = overlap_status),
      alpha = 0.6,
      size = 1.8,
      show.legend = FALSE
    ) +
    geom_polygon(
      data = hulls_df,
      aes(x = PC1, y = PC2, color = origin, fill = origin),
      alpha = 0.1,
      linewidth = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = color_origins) +
    scale_fill_manual(values = color_origins) +
    geom_point(
      data = means_df,
      aes(x = mean_PC1, y = mean_PC2, fill = origin),
      shape = 25,
      size = 3.5,
      color = "black",
      show.legend = FALSE
    ) +
    # arrow segments
    geom_segment(
      data = loadings_df,
      aes(x = 0, y = 0,
          xend = PC1 * arrow_mult,
          yend = PC2 * arrow_mult),
      arrow = arrow(length = unit(0.3, "cm")),
      color = "grey15"
    ) +
    # arrow text
    geom_text_repel(
      data = loadings_df,
      aes(x = PC1 * arrow_mult,
          y = PC2 * arrow_mult,
          label = variable),
      size = 5.5,
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = Inf,
      segment.size = 0.2,
      seed = 42,
      fontface = "bold"     
    ) +
    labs(
      title = bquote(.(plot_letter) ~ italic(.(species_names[[species]]))), 
      x     = paste0("PC1 (", round(summary(pca_all)$importance[2,1]*100,1), "%)"),
      y     = paste0("PC2 (", round(summary(pca_all)$importance[2,2]*100,1), "%)"),
      color = "Origin",
      fill  = "Origin"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title   = element_text(size = 16),
      plot.title   = element_text(size = 20),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white")
    )
  

  print(p[[species]])

  ggsave(filename = paste0("PCA_results/PCA_plots/PCA_gbif_vs_hlnug_",species,".png"), plot = p[[species]], width = 10, height = 10, dpi = 400, bg = "white")

  # Niche Overlap calculation after Broennimann
  ##########################################################################
  ##########################################################################
  set.seed(798454087)

  coords <- terra::xyFromCell(predictors_local, 1:ncell(predictors_local))
  predictors_local_coords <- cbind(coords, predictors_local)

  occ_gbif_vec <- terra::vect(occ_gbif, geom = c("decimalLongitude","decimalLatitude"), crs = crs(predictors_local))
  occ_hlnug_vec <- terra::vect(occ_hlnug, geom = c("decimalLongitude","decimalLatitude"), crs = crs(predictors_local))

  coords_gbif <- terra::crds(occ_gbif_vec)
  cells_gbif <- terra::cellFromXY(predictors_local, coords_gbif)

  coords_hlnug <- terra::crds(occ_hlnug_vec)
  cells_hlnug <- terra::cellFromXY(predictors_local, coords_hlnug)

  df_predictors_local <- as.data.frame(predictors_local, xy = TRUE, na.rm = FALSE)
  df_predictors_local$origin <- "none"
  df_predictors_local$origin[cells_gbif] <- "GBIF"
  df_predictors_local$origin[cells_hlnug] <- "HLNUG"
  df_predictors_local$origin[intersect(cells_gbif, cells_hlnug)] <- "BOTH"

  df_predictors_local_num <- df_predictors_local %>%
  dplyr::select(any_of(c("aspect","DGM","flowing_water_distance",
                         "railroads_density","railroads_distance",
                         "roads_distance","standing_water_distance",
                         "travel_time","twi","urban_builup_density", "origin")))
  df_predictors_local_num_clean <- df_predictors_local_num[complete.cases(df_predictors_local_num) &
                                                 !apply(is.infinite(as.matrix(df_predictors_local_num)), 1, any), ]

  df_predictors_pca <- df_predictors_local_num_clean %>%
    dplyr::select(-origin)
  
  pca_all_hesse <- prcomp(df_predictors_pca, center = FALSE, scale. = TRUE)

  scores_proc <- as.data.frame(pca_all_hesse$x[, 1:2]) %>%
    bind_cols(df_predictors_local_num_clean %>% dplyr::select(origin))
  
  scores_all_bg_hesse <- scores_proc %>%
    dplyr::select(-origin)
  scores_hlnug_bg_hesse <- scores_proc %>%
    filter(origin %in% c("HLNUG", "BOTH")) %>%
    dplyr::select(-origin)
  scores_gbif_bg_hesse <- scores_proc %>%
    filter(origin %in% c("GBIF", "BOTH")) %>%
    dplyr::select(-origin)
  
    bg <- "hesse"
    glob <- get(paste0("scores_all_bg_",bg))
    sp_gbif <- get(paste0("scores_gbif_bg_",bg))
    sp_hlnug <- get(paste0("scores_hlnug_bg_",bg))
    pca_to_use <- if (bg == "occ") pca_all else pca_all_hesse
   
    #Project onto 2PCA-Axis density  grid
    grid.hlnug <-ecospat.grid.clim.dyn(glob=glob, glob1=glob, sp=sp_hlnug, R=150, th.sp=0)
    grid.gbif <-ecospat.grid.clim.dyn(glob=glob, glob1=glob, sp=sp_gbif, R=150, th.sp=0)

  #Schoeners D
  ##########################################################################
  D.overlap<- ecospat.niche.overlap (grid.hlnug, grid.gbif, cor = TRUE)$D; D.overlap
  D_list[[species]] <- D.overlap
  
  #Warren Niche Equicalency test
  ##########################################################################
 rep = 1000
 
 eq_test <- ecospat.niche.equivalency.test(
   grid.hlnug,
   grid.gbif,
   rep = rep,
   overlap.alternative = "lower",    # test if overlap is LESS than expected
 )
eq_test_list[[species]] <- eq_test

 #save results
 niche_summary <- bind_rows(
   data.frame(
     Test = "Equivalency",
     Metric = c("D", "I", "Expansion", "Stability", "Unfilling"),
     Observed = c(eq_test$obs$D, eq_test$obs$I, eq_test$obs$expansion, eq_test$obs$stability, eq_test$obs$unfilling),
     p_value = c(eq_test$p.D, eq_test$p.I, eq_test$p.expansion, eq_test$p.stability, eq_test$p.unfilling),
     Direction = "NA"
   )
 )

 # Save as CSV
 write.csv(niche_summary,paste0("PCA_results/niche_overlap_summary_table",species,".csv"), row.names = FALSE)
}

#One DF for loadings
##################################
combined_loadings <- do.call(cbind, lapply(names(loadings), function(sp) {
  df <- loadings[[sp]]
  
  # If the df already has a 'Variable' column, remove it
  df <- df[ , !colnames(df) %in% c("Variable"), drop = FALSE]
  
  # Rename columns with species name
  colnames(df) <- paste(colnames(df), sp, sep = "_")
  df
}))

# Add a single variable column from the first species
combined_loadings <- data.frame(Variable = rownames(loadings[[1]]), combined_loadings, row.names = NULL) %>%
  select(-Variable, -variable_3190653, -variable_3034825) %>%
  rename(
    Variable = variable_2891770)

combined_loadings[order(-abs(combined_loadings$PC1_2891770)), ]

write.csv(combined_loadings,paste0("PCA_results/PCA_loadings/PCA_loadings_",species,"_all.csv"), row.names = FALSE)

#combined PCA Plot
###################################
D_value <- list()
p_value <- list()

for (species in species_list){
values <- read.csv(paste0("PCA_results/niche_overlap_summary_table",species,".csv"))

D_value[[species]] <- values$Observed[values$Metric == "D"]
p_value[[species]] <- values$p_value[values$Metric == "D"]

}

global_xlim <- range(unlist(xlims))
global_ylim <- range(unlist(ylims))

padding_x <- 0.2 * diff(global_xlim)
padding_y <- 0.2 * diff(global_ylim)

global_xlim_pad <- c(global_xlim[1] - padding_x, global_xlim[2] + padding_x)
global_ylim_pad <- c(global_ylim[1] - padding_y, global_ylim[2] + padding_y)

p_fixed <- lapply(p, function(pl) {
  pl + coord_cartesian(xlim = global_xlim_pad, ylim = global_ylim_pad, expand = TRUE)
})

p_fixed <- lapply(p_fixed, function(pl) {
  pl +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1), # thinner black frame
      panel.background = element_rect(fill = "white") # optional: white background
    )
})

#add Schoeners D values and p values to plot
add_schoeners_box <- function(plot, species_id) {
  D_obs <- D_value[[species_id]]
  p_D <- p_value[[species_id]]
  
  var_text <- paste0(
    "Schoener's D: ", round(D_obs, 3), "\n",
    "Equivalency p-value: ", round(p_D, 3)
  )

  plot + annotate(
    "label",
    x = -Inf, y = -Inf,
    label = var_text,
    hjust = 0, vjust = 0,
    size = 7,                  # bigger box text
    fill = "white",
    color = "black",
    label.size = 0.2,
    label.r = grid::unit(0, "pt"),
    lineheight = 1.2
  )
}

# Apply to all plots
  for (sp in names(p_fixed)) {
    p_fixed[[sp]] <- add_schoeners_box(p_fixed[[sp]], sp)
  }


plots_list <- list(
  A = p_fixed[["3190653"]],  # A. altissima (a)
  B = p_fixed[["3034825"]],  # Heracleum (b)
  C = p_fixed[["2891770"]]   # Impatiens (c)
)

# Design als Text: erste Reihe "AB", zweite Reihe "C." (Punkt = leere Zelle)
design <- "
AB
C.
"

# wrap_plots mit dem Text-Design — danach plot_layout für widths/heights/Abstände
combined_plot <- wrap_plots(plots_list, design = design) +
  plot_layout(widths = c(1, 1), heights = c(1, 1), guides = "collect") &
  theme(plot.margin = margin(5,5,5,5))

# anzeigen
combined_plot

ggsave(
  filename = "PCA_results/PCA_plots/PCA_all_three_panel.png",  # file path/name
  plot = combined_plot,
  width = 20,                      # very wide (in inches)
  height = 20,                      # adjust height as needed
  units = "in",                     # dimensions in inches
  dpi = 600,                        # high resolution for print quality
  bg = "white"                      # white background
)

