# Full script: read local CSVs, compute per-observation means, then
# produce combined continuous facets + categorical group barplots side-by-side
# with a single shared legend (one PNG per species).

# ---------------------------
# Libraries
# ---------------------------
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Model_results/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  Model_results/response_data/res__local_<origin>_<species>_<run>.csv (10_SDM_Code_Model_final.R)
# Output: Model_results/response_data/Plots/response_combined_side_ALL_SPECIES.png (Suppl. App. 5.1),
#         response_combined_side-<species>.png (Suppl. App. 5.2-5.4)

# ---------------------------
# Paths & species
# ---------------------------
path <- "Model_results/response_data"
species_list <- c("2891770","3190653","3034825")
save_path <- file.path(path, "Plots")
dir.create(save_path, showWarnings = FALSE)

# ---------------------------
# Continuous predictor renaming
# ---------------------------
rename_map <- c(
  aspect = "Aspect",
  DGM = "DEM",
  flowing_water_distance = "Distance_Flowing_Water",
  railroads_density = "Railroad_Density",
  railroads_distance = "Distance_Railroads",
  roads_distance = "Distance_Roads",
  standing_water_distance = "Distance_Standing_Water",
  travel_time = "Travel_Time_to_Cities",
  twi = "TWI",
  urban_builup_density = "Urban_Density"
)

# ---------------------------
# Read only local CSVs and assemble raw_data
# ---------------------------
all_data <- list()

for (sp in species_list) {
  files <- list.files(path, pattern = paste0("local.*", sp), full.names = TRUE)
  if (length(files) == 0) next
  sp_data <- lapply(files, function(f) {
    # read_csv2: files use comma as decimal mark
    df <- read_csv2(f, locale = locale(decimal_mark = ","))
    origin <- str_extract(basename(f), "(local_gbif|local_hlnug|local_combined)")
    iter   <- as.integer(str_match(basename(f), paste0(sp, "_([0-9]+)\\.csv"))[,2])
    df$origin <- origin
    df$species <- sp
    df$iteration <- iter
    df$rowid <- seq_len(nrow(df))
    df
  })
  all_data[[sp]] <- bind_rows(sp_data)
}

raw_data <- bind_rows(all_data)

# ---------------------------
# Compute per-observation means across iterations
# ---------------------------
first_non_na <- function(x) {
  x_non_na <- x[!is.na(x)]
  if (length(x_non_na) == 0) return(NA)
  x_non_na[1]
}

mean_data <- raw_data %>%
  group_by(origin, species, rowid) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    across(where(~ is.character(.x) || is.factor(.x)), ~ first_non_na(.x)),
    .groups = "drop"
  )

# drop rowid
mean_data <- dplyr::select(mean_data, -rowid)

# ---------------------------
# Categorical & pretty naming helpers
# ---------------------------
corine_map <- c(
  CORINE_X12 = "Industrial units",
  CORINE_X14 = "Artificial vegetated areas",
  CORINE_X21 = "Arable land",
  CORINE_X22 = "Permanent crops",
  CORINE_X23 = "Pastures",
  CORINE_X24 = "Heterogeneous agricultural land",
  CORINE_X31 = "Forests",
  CORINE_X32 = "Shrub & herbaceous vegetations",
  CORINE_X33 = "Open spaces with little vegetation",
  CORINE_X51 = "Water bodies"
)

flood_map <- c(
  floodplains_X2 = "Low",
  floodplains_X3 = "Medium",
  floodplains_X4 = "High"
)

pretty_whc  <- c(
  whc_Siedlung.und.Verkehr = "Settlement & Transportation",
  whc_extrem.gering..0...80.mm. = "Extremely low",
  whc_gering...110...150.mm. = "Low",
  whc_hoch...200...260.mm. = "High",    
  whc_mittel...150...200.mm. = "Medium",
  whc_sehr.gering...80...110.mm. = "Very low",
  whc_sehr.hoch...260.mm. = "Very high"
)

pretty_heat <-  c(
  heat_index_heiß = "hot",
  heat_index_kühl = "cool",
  heat_index_mäßig = "moderate",
  heat_index_sehr.warm = "very warm",
  heat_index_warm = "warm")

# group patterns for categorical sets
group_patterns <- list(
  CORINE = "^CORINE_",
  floodplains = "^floodplains_",
  whc = "^whc_",
  heat_index = "^heat_index_"
)

# ---------------------------
# Build continuous long dataframe
# ---------------------------
continuous_preds <- names(rename_map)
cont_long <- mean_data %>%
  select(any_of(c("origin", "species", "prob", continuous_preds))) %>%
  pivot_longer(
    cols = all_of(intersect(continuous_preds, names(.))),
    names_to = "predictor",
    values_to = "value"
  ) %>%
  mutate(type = "continuous")

# ---------------------------
# Build categorical assigned-by-max summaries
# ---------------------------
cat_rows <- list()
for (grp in names(group_patterns)) {
  pat <- group_patterns[[grp]]
  cols_grp <- grep(pat, names(mean_data), value = TRUE)
  if (length(cols_grp) == 0) next
  
  mat_all <- as.matrix(mean_data[, cols_grp, drop = FALSE])
  mat_safe <- ifelse(is.na(mat_all), -Inf, mat_all)
  max_idx <- apply(mat_safe, 1, which.max)
  all_na_rows <- apply(mat_all, 1, function(r) all(is.na(r)))
  assigned <- cols_grp[max_idx]
  assigned[all_na_rows] <- NA_character_
  
  tmp <- mean_data %>% mutate(assigned = assigned)
  summary_assigned <- tmp %>%
    filter(!is.na(assigned)) %>%
    group_by(species, assigned, origin) %>%
    summarise(mean_prob = mean(prob, na.rm = TRUE), n = n(), .groups = "drop") %>%
    mutate(group = grp)
  cat_rows[[grp]] <- summary_assigned
}

cat_summary_all <- if (length(cat_rows) > 0) bind_rows(cat_rows) else tibble()

# prepare pretty class labels for categorical plotting
if (nrow(cat_summary_all) > 0) {
  cat_plot_df <- cat_summary_all %>%
    mutate(
      predictor = group,
      class_pretty = case_when(
        grepl("^CORINE_", assigned) ~ 
          ifelse(assigned %in% names(corine_map), corine_map[assigned], sub("^CORINE_", "", assigned)),
        grepl("^floodplains_", assigned) ~ 
          ifelse(assigned %in% names(flood_map), flood_map[assigned], sub("^floodplains_", "", assigned)),
        grepl("^whc_", assigned) ~ 
          ifelse(assigned %in% names(pretty_whc), pretty_whc[assigned], sub("^whc_", "", assigned)),
        grepl("^heat_index_", assigned) ~ 
          ifelse(assigned %in% names(pretty_heat), pretty_heat[assigned], sub("^heat_index_", "", assigned)),
        TRUE ~ assigned
      )
    )
} else {
  cat_plot_df <- tibble(species = character(), assigned = character(), origin = character(), mean_prob = numeric(), n = integer(), group = character(), predictor = character(), class_pretty = character())
}

# facet labels map (continuous + categorical groups)
facet_label_map <- c(
  rename_map,
  CORINE = "CORINE land cover",
  floodplains = "Floodplains",
  whc = "WHC",
  heat_index = "Heat Index"
)

# ---------------------------
# Build combined side-by-side figure per species (continuous left, categorical right),
# no empty facets, only categorical legend
# ---------------------------
origin_labels <- c(
  local_combined = "combined",
  local_gbif    = "Citizen Science",
  local_hlnug   = "State agency"
)

ncol_facets <- 3

for (sp in species_list) {
  # ---------------------------
  # Continuous data (left)
  # ---------------------------
  cont_df <- cont_long %>%
    filter(species == sp) %>%
    group_by(predictor) %>%
    filter(!all(is.na(prob)) & !all(is.na(value))) %>%
    ungroup()
  
  # if there is really nothing to plot, skip species
  if (nrow(cont_df) == 0 && (!exists("cat_plot_df") || nrow(cat_plot_df %>% filter(species == sp)) == 0)) {
    message("No data to plot for species ", sp, "; skipping.")
    next
  }
  
  p_cont <- ggplot(cont_df, aes(x = value, y = prob, color = origin)) +
    geom_smooth(method = "loess", se = TRUE) +
    facet_wrap(
      ~ predictor,
      scales = "free_x",
      ncol = ncol_facets,
      labeller = as_labeller(facet_label_map)
    ) +
    labs(
      x = NULL,
      y = "Mean suitability",   # same as categorical
      color = "Origin"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    # hide continuous legend, but keep consistent labels
    scale_color_discrete(
      breaks = names(origin_labels),
      labels = origin_labels,
      guide = "none"
    )
  
  # ---------------------------
  # Categorical data (right)
  # ---------------------------
  cat_df  <- cat_plot_df %>%
    filter(species == sp)
  
  cat_plots <- list()
  
  for (grp in names(group_patterns)) {
    grp_df <- cat_df %>%
      filter(group == grp) %>%
      filter(!is.na(mean_prob))   # drop empty/NA rows
    
    if (nrow(grp_df) == 0) next
    
    # order classes by overall mean_prob
    order_classes <- grp_df %>%
      group_by(class_pretty) %>%
      summarise(overall = mean(mean_prob, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(overall)) %>%
      pull(class_pretty)
    
    grp_df <- grp_df %>%
      mutate(class_pretty = factor(class_pretty, levels = unique(order_classes)))
    
    # use facet_label_map for the title so capitalization matches
    grp_title <- if (!is.null(facet_label_map[[grp]])) facet_label_map[[grp]] else grp
    
    p_cat <- ggplot(grp_df, aes(x = class_pretty, y = mean_prob, fill = origin)) +
      geom_col(
        position = position_dodge(width = 0.8),
        width = 0.7,
        na.rm = TRUE
      ) +
      geom_text(
        aes(label = sprintf("%.2f", mean_prob), group = origin),
        position = position_dodge(width = 0.8),
        vjust = -0.4,
        size = 3,
        na.rm = TRUE
      ) +
      labs(
        title = grp_title,
        x = NULL,
        y = "Mean suitability",   # Mean suitability always on y-axis
        fill = "Origin"
      ) +
      scale_fill_discrete(
        breaks = names(origin_labels),
        labels = origin_labels
      ) +
      theme_minimal(base_size = 10) +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      )
    
    # IMPORTANT: no coord_flip() -> Mean suitability stays on y-axis
    cat_plots[[grp]] <- p_cat
  }
  
  # ---------------------------
  # Combine and save
  # ---------------------------
  if (length(cat_plots) > 0) {
    cats_wrap <- wrap_plots(cat_plots, ncol = 1)
    final <- (p_cont | cats_wrap) +
      plot_layout(widths = c(3, 1), guides = "collect")
    
    final <- final & theme(legend.position = "bottom")
  } else {
    # only continuous side, still without its own legend
    final <- p_cont
  }
  
  n_continuous <- length(unique(cont_df$predictor))
  nrow_facets <- ceiling((n_continuous) / ncol_facets)
  cont_h <- max(4, nrow_facets * 4)
  w <- 18
  h <- cont_h
  
  out_name <- file.path(save_path, paste0("response_combined_side-", sp, ".png"))
  ggsave(out_name, final, width = w, height = h, dpi = 300, bg = "white")
  message("Saved combined side-by-side panel for species ", sp, " -> ", out_name)
}


# =========================================================
# OVERALL SUMMARY PLOT: species-averaged responses per origin
# (same style as per-species panel, but one figure for all species)
# =========================================================

# ---------------------------
# Continuous: average over species
# ---------------------------
cont_long_all <- cont_long %>%
  group_by(origin, predictor, value) %>%
  summarise(
    prob = mean(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(predictor) %>%
  filter(!all(is.na(prob)) & !all(is.na(value))) %>%
  ungroup()

p_cont_all <- ggplot(cont_long_all, aes(x = value, y = prob, color = origin)) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(
    ~ predictor,
    scales = "free_x",
    ncol = ncol_facets,
    labeller = as_labeller(facet_label_map)
  ) +
  labs(
    x = NULL,
    y = "Mean suitability",
    color = "Origin"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  scale_color_discrete(
    breaks = names(origin_labels),
    labels = origin_labels
  )

# ---------------------------
# Categorical: average over species
# ---------------------------
cat_df_all <- cat_plot_df %>%
  group_by(origin, group, predictor, class_pretty) %>%
  summarise(
    mean_prob = mean(mean_prob, na.rm = TRUE),
    n = sum(n, na.rm = TRUE),
    .groups = "drop"
  )

cat_plots_all <- list()

for (grp in names(group_patterns)) {
  grp_df <- cat_df_all %>%
    filter(group == grp) %>%
    filter(!is.na(mean_prob))
  
  if (nrow(grp_df) == 0) next
  
  # order classes by overall mean_prob
  order_classes <- grp_df %>%
    group_by(class_pretty) %>%
    summarise(overall = mean(mean_prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(overall)) %>%
    pull(class_pretty)
  
  grp_df <- grp_df %>%
    mutate(class_pretty = factor(class_pretty, levels = unique(order_classes)))
  
  grp_title <- if (!is.null(facet_label_map[[grp]])) facet_label_map[[grp]] else grp
  
  p_cat_all <- ggplot(grp_df, aes(x = class_pretty, y = mean_prob, fill = origin)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      na.rm = TRUE
    ) +
    geom_text(
      aes(label = sprintf("%.2f", mean_prob), group = origin),
      position = position_dodge(width = 0.8),
      vjust = -0.4,
      size = 3,
      na.rm = TRUE
    ) +
    labs(
      title = grp_title,
      x = NULL,
      y = "Mean suitability",
      fill = "Origin"
    ) +
    scale_fill_discrete(
      breaks = names(origin_labels),
      labels = origin_labels
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  
  cat_plots_all[[grp]] <- p_cat_all
}

# ---------------------------
# Combine and save overall figure
# ---------------------------
if (length(cat_plots_all) > 0) {
  cats_wrap_all <- wrap_plots(cat_plots_all, ncol = 1)
  final_all <- (p_cont_all | cats_wrap_all) +
    plot_layout(widths = c(3, 1), guides = "collect")
  
  final_all <- final_all & theme(legend.position = "bottom")
} else {
  final_all <- p_cont_all
}

# Size: scale with number of continuous predictors
n_continuous_all <- length(unique(cont_long_all$predictor))
nrow_facets_all <- ceiling(n_continuous_all / ncol_facets)
cont_h_all <- max(4, nrow_facets_all * 4)
w_all <- 18
h_all <- cont_h_all

out_name_all <- file.path(save_path, "response_combined_side_ALL_SPECIES.png")
ggsave(out_name_all, final_all, width = w_all, height = h_all, dpi = 300, bg = "white")
message("Saved overall combined side-by-side panel for ALL species -> ", out_name_all)
