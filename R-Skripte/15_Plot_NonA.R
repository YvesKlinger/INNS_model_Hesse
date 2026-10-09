# ============================
# Libraries
# ============================
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Model_results/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  Model_results/Model_analyses/Non_applicable_analyses/LM_summary__local_<origin>_<species>_mean.csv
#         (10_SDM_Code_Model_final.R)
# Output: .../Non_applicable_analyses/plots/Dumbbell_all_species.png (Figure 3),
#         .../Non_applicable_analyses/plots/Dumbbell_meaned_over_origins.png (Suppl. App. 6)

# ============================
# Path to the *_mean.csv files
# ============================
path_nonapp <- "Model_results/Model_analyses/Non_applicable_analyses"
dir.create(file.path(path_nonapp, "plots"), showWarnings = FALSE)

# Only mean files, e.g.:
# LM_summary__local_hlnug_3190653_mean.csv
na_files <- list.files(
  path_nonapp,
  pattern = "^LM_summary__local_.*_mean\\.csv$",
  full.names = TRUE
)

stopifnot("No LM_summary__local_*_mean.csv files found!" = length(na_files) > 0)

# ============================
# Helper: read one mean file & extract origin/species
# ============================
read_mean_file <- function(fp) {
  bn <- basename(fp)
  # Expected pattern: LM_summary__local_<origin>_<species>_mean.csv
  m <- stringr::str_match(bn, "^LM_summary__local_([^_]+)_([0-9]+)_mean\\.csv$")
  origin  <- m[, 2]
  species <- m[, 3]
  
  if (is.na(origin) || is.na(species)) {
    stop("Could not parse origin/species from filename: ", bn)
  }
  
  dat <- readr::read_csv(fp, show_col_types = FALSE)
  cols <- tolower(names(dat))
  
  # predictor column (should be "predictor"; matched case-insensitively)
  pred_col <- names(dat)[cols == "predictor"][1]
  
  # coefficient: try mean_estimate first, then estimate
  coef_col <- names(dat)[cols %in% c("mean_estimate", "estimate_mean", "estimate")][1]
  
  # importance: try mean_weight first, then weight
  w_col <- names(dat)[cols %in% c("mean_weight", "weight")][1]
  
  if (is.na(pred_col) || is.na(coef_col) || is.na(w_col)) {
    stop(
      "Missing required columns in file: ", bn,
      " (need predictor + coefficient + importance). Got columns: ",
      paste(names(dat), collapse = ", ")
    )
  }
  
  dat %>%
    dplyr::transmute(
      species   = species,
      origin    = origin,
      predictor = .data[[pred_col]],
      coef      = as.numeric(.data[[coef_col]]),
      weight    = as.numeric(.data[[w_col]])
    )
}

# ============================
# Read & combine all mean files
# ============================
raw_tbl <- purrr::map_dfr(na_files, read_mean_file)

# (optional) map origin codes to nicer labels
origin_map <- c(
  "combined" = "COM",
  "gbif"     = "CS",
  "hlnug"    = "StAg"
)

raw_tbl <- raw_tbl %>%
  dplyr::mutate(
    origin = origin_map[origin]
  )

# ============================
# Mean over files per species × origin × predictor
# ============================
big_table <- raw_tbl %>%
  dplyr::group_by(species, origin, predictor) %>%
  dplyr::summarise(
    mean_coef   = mean(coef,   na.rm = TRUE),
    sd_coef     = sd(coef,     na.rm = TRUE),
    mean_weight = mean(weight, na.rm = TRUE),
    sd_weight   = sd(weight,   na.rm = TRUE),
    n_files     = dplyr::n(),
    .groups     = "drop"
  ) %>%
  dplyr::arrange(species, origin, dplyr::desc(mean_weight)) %>%
  select(-sd_coef, -sd_weight, -n_files)

print(big_table, n = 50)

# ============================
# Mean over species for each origin × predictor
# ============================

origin_summary <- big_table %>%
  group_by(origin, predictor) %>%
  summarise(
    mean_coef_origin   = mean(mean_coef,   na.rm = TRUE),
    mean_weight_origin = mean(mean_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(origin, desc(mean_weight_origin))

print(origin_summary, n = 50)

# ============================
# Keep only predictors above the median importance
# ============================

# global median across origins and predictors
global_median <- median(origin_summary$mean_weight_origin, na.rm = TRUE)

message("Global median of mean_weight_origin = ", round(global_median, 4))

origin_summary_filtered <- origin_summary %>%
  filter(mean_weight_origin > global_median) %>%
  arrange(origin, desc(mean_weight_origin))

print(origin_summary_filtered, n = 50)

# ============================
# Pretty labels
# ============================
pretty_var <- function(v) {
  
  # CORINE
  if (grepl("^CORINE_X", v)) {
    code <- sub("^CORINE_X", "", v)
    name <- dplyr::case_when(
      code == "12" ~ "Ind/Comm/Transp",
      code == "14" ~ "Artificial vegetated",
      code == "21" ~ "Arable",
      code == "22" ~ "Perm crops",
      code == "23" ~ "Pastures",
      code == "24" ~ "Het. agri",
      code == "31" ~ "Forests",
      code == "32" ~ "Scrub/herb",
      code == "51" ~ "Inland waters",
      TRUE ~ paste0("CORINE ", code)
    )
    return(paste0("COR - ", name))
  }
  
  # Floodplains
  if (grepl("^floodplains_X", v)) {
    k <- sub("^floodplains_X", "", v)
    level <- dplyr::case_when(
      k == "4" ~ "High",
      k == "3" ~ "Medium",
      k == "2" ~ "Low",
      k == "1" ~ "None",
      TRUE ~ k
    )
    return(paste0("FP - ", level))
  }
  
  # Heat Index
  if (grepl("^heat_index", v)) {
    cls <- sub("^heat_index_?", "", v)
    cls <- gsub("\\.", " ", cls)
    eng <- dplyr::case_when(
      cls %in% c("kühl","kuehl") ~ "Cool",
      cls %in% c("mäßig","maessig") ~ "Moderate",
      cls == "warm" ~ "Warm",
      cls == "sehr warm" ~ "Very warm",
      cls %in% c("heiß","heiss") ~ "Hot",
      TRUE ~ stringr::str_to_title(cls)
    )
    return(paste0("HI - ", eng))
  }
  
  # WHC
  if (grepl("^whc", v)) {
    cls <- v
    cls <- sub("^whc_", "", cls)
    cls <- gsub("\\.+", " ", cls)
    cls <- trimws(sub("\\(.*$", "", cls))
    
    eng <- dplyr::case_when(
      grepl("^extrem gering|sehr gering", cls) ~ "Very low",
      grepl("^gering", cls)                    ~ "Low",
      grepl("^mittel", cls)                    ~ "Medium",
      grepl("^sehr hoch", cls)                 ~ "Very high",
      grepl("^hoch", cls)                      ~ "High",
      grepl("Siedlung und Verkehr|Siedlung.und.Verkehr", v, ignore.case = TRUE) ~ "Settlement & Transport",
      TRUE ~ stringr::str_to_title(cls)
    )
    return(paste0("WHC - ", eng))
  }
  
  # Hydro / Infra
  if (v %in% c("flowing_water_distance","standing_water_distance",
               "roads_distance","railroads_distance","railroads_density",
               "urban_builup_density","travel_time")) {
    s <- gsub("_", " ", v)
    return(stringr::str_to_sentence(s))
  }
  
  # Topography
  if (v %in% c("DGM","twi","aspect")) {
    return(dplyr::case_when(
      v == "DGM"    ~ "DEM",
      v == "twi"    ~ "TWI",
      v == "aspect" ~ "Aspect",
      TRUE ~ v
    ))
  }
  
  # default fallback
  gsub("_", " ", v)
}

origin_summary_filtered <- origin_summary_filtered %>%
  mutate(
    var_pretty = vapply(predictor, pretty_var, character(1))
  )

# ============================
# Prepare data for dumbbell plot (origins)
# ============================
coef_range <- origin_summary_filtered %>%
  mutate(
    group = dplyr::case_when(
      grepl("^COR -", var_pretty) ~ "COR",
      grepl("^WHC -", var_pretty) ~ "WHC",
      grepl("^HI -",  var_pretty) ~ "HI",
      grepl("^FP -",  var_pretty) ~ "FP",
      TRUE                        ~ "Other"
    )
  ) %>%
  group_by(group, var_pretty) %>%
  summarise(
    min_coef = min(mean_coef_origin, na.rm = TRUE),
    max_coef = max(mean_coef_origin, na.rm = TRUE),
    # backup ordering by largest absolute coefficient
    order    = max(abs(mean_coef_origin), na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  # group order: OTHER (continuous) first, then COR, FP, WHC, HI
  mutate(
    group = factor(group, levels = c("Other", "COR", "FP", "WHC", "HI")),
    height_rank = dplyr::case_when(
      # FP levels (within FP)
      group == "FP"  & var_pretty == "FP - High"   ~ 3,
      group == "FP"  & var_pretty == "FP - Medium" ~ 2,
      group == "FP"  & var_pretty == "FP - Low"    ~ 1,
      group == "FP"  & var_pretty == "FP - None"   ~ 0,
      # WHC levels (Settlement & Transport as lowest)
      group == "WHC" & var_pretty == "WHC - Settlement & Transport" ~ 0,
      group == "WHC" & var_pretty == "WHC - Very low"               ~ 1,
      group == "WHC" & var_pretty == "WHC - Low"                    ~ 2,
      group == "WHC" & var_pretty == "WHC - Medium"                 ~ 3,
      group == "WHC" & var_pretty == "WHC - High"                   ~ 4,
      group == "WHC" & var_pretty == "WHC - Very high"              ~ 5,
      # HI levels
      group == "HI"  & var_pretty == "HI - Cool"      ~ 0,
      group == "HI"  & var_pretty == "HI - Warm"      ~ 1,
      group == "HI"  & var_pretty == "HI - Very warm" ~ 2,
      TRUE ~ NA_real_
    ),
    # fallback to coefficient-based order where not specified
    height_rank = ifelse(is.na(height_rank), order, height_rank)
  ) %>%
  arrange(group, desc(height_rank)) %>%
  mutate(
    var_fac = factor(var_pretty, levels = rev(unique(var_pretty)))
  )

df_plot <- origin_summary_filtered %>%
  inner_join(coef_range %>% select(var_pretty, var_fac, min_coef, max_coef, group),
             by = "var_pretty")

# ============================
# Dumbbell plot (all origins)
# ============================
gg_coef_origins <- ggplot(df_plot, aes(y = var_fac)) +
  geom_segment(aes(x = min_coef, xend = max_coef, yend = var_fac),
               colour = "grey80", linewidth = 0.6) +
  # main reference line at 0
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(x = mean_coef_origin, colour = origin),
             size = 2.3) +
  scale_colour_manual(
    values = c(COM = "#d95f02", CS = "#1b9e77" , StAg = "#7570b3"),
    name   = "Origin",
    labels = c(COM = "Combined", CS = "Citizen Science", StAg = "State agency")
  ) +
  labs(
    x = "Coefficient (DI)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background    = element_rect(fill = "white", colour = NA),
    plot.background     = element_rect(fill = "white", colour = NA),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    axis.text.y         = element_text(size = 7),
    axis.title.x        = element_text(size = 10),
    legend.position     = "bottom"
  )

ggsave(
  file.path(path_nonapp, "plots/Dumbbell_all_species.png"),
  gg_coef_origins, width = 9, height = 6, dpi = 300, bg = "white"
)

gg_coef_origins

# ============================
# Species names (used later)
# ============================
species_names <- c(
  "2891770" = "I. glandulifera",
  "3190653" = "A. altissima",
  "3034825" = "H. mantegazzianum"
)

# ============================
# SPECIES-BASED DUMBBELL (origins averaged, points = species)
# ============================

species_summary <- big_table %>%
  group_by(species, predictor) %>%
  summarise(
    mean_coef_species   = mean(mean_coef,   na.rm = TRUE),
    mean_weight_species = mean(mean_weight, na.rm = TRUE),
    .groups = "drop"
  )

species_summary <- species_summary %>%
  mutate(
    species_label = species_names[species]
  )

species_summary_filtered <- species_summary %>%
  filter(mean_weight_species > global_median) %>%
  mutate(
    var_pretty = vapply(predictor, pretty_var, character(1))
  )

coef_range_species <- species_summary_filtered %>%
  mutate(
    group = case_when(
      grepl("^COR -", var_pretty) ~ "COR",
      grepl("^WHC -", var_pretty) ~ "WHC",
      grepl("^HI -",  var_pretty) ~ "HI",
      grepl("^FP -",  var_pretty) ~ "FP",
      TRUE                        ~ "Other"
    )
  ) %>%
  group_by(group, var_pretty) %>%
  summarise(
    min_coef = min(mean_coef_species, na.rm = TRUE),
    max_coef = max(mean_coef_species, na.rm = TRUE),
    order    = max(abs(mean_coef_species), na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    group = factor(group, levels = c("Other", "COR", "FP", "WHC", "HI")),
    height_rank = case_when(
      group == "FP"  & var_pretty == "FP - High"   ~ 3,
      group == "FP"  & var_pretty == "FP - Medium" ~ 2,
      group == "FP"  & var_pretty == "FP - Low"    ~ 1,
      group == "FP"  & var_pretty == "FP - None"   ~ 0,
      group == "WHC" & var_pretty == "WHC - Settlement & Transport" ~ 0,
      group == "WHC" & var_pretty == "WHC - Very low"               ~ 1,
      group == "WHC" & var_pretty == "WHC - Low"                    ~ 2,
      group == "WHC" & var_pretty == "WHC - Medium"                 ~ 3,
      group == "WHC" & var_pretty == "WHC - High"                   ~ 4,
      group == "WHC" & var_pretty == "WHC - Very high"              ~ 5,
      group == "HI"  & var_pretty == "HI - Cool"      ~ 0,
      group == "HI"  & var_pretty == "HI - Warm"      ~ 1,
      group == "HI"  & var_pretty == "HI - Very warm" ~ 2,
      TRUE ~ NA_real_
    ),
    height_rank = ifelse(is.na(height_rank), order, height_rank)
  ) %>%
  arrange(group, desc(height_rank)) %>%
  mutate(
    var_fac = factor(var_pretty, levels = rev(unique(var_pretty)))
  )

df_species_plot <- species_summary_filtered %>%
  inner_join(
    coef_range_species %>% select(var_pretty, var_fac, min_coef, max_coef, group),
    by = "var_pretty"
  )

species_cols <- c(
  "I. glandulifera"   = "#7E6AAE",
  "A. altissima"      = "#D4C95F",
  "H. mantegazzianum" = "#4F8F6B"
)

gg_coef_species <- ggplot(df_species_plot, aes(y = var_fac)) +
  geom_segment(aes(x = min_coef, xend = max_coef, yend = var_fac),
               colour = "grey80", linewidth = 0.6) +
  # main reference line at 0
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(x = mean_coef_species, colour = species_label),
             size = 2.3) +
  scale_colour_manual(
    values = species_cols,
    name   = "Species"
  ) +
  labs(
    x = "Coefficient (DI)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background    = element_rect(fill = "white", colour = NA),
    plot.background     = element_rect(fill = "white", colour = NA),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    axis.text.y         = element_text(size = 7),
    axis.title.x        = element_text(size = 10),
    legend.position     = "bottom"
  )

gg_coef_species

ggsave(
  file.path(path_nonapp, "plots/Dumbbell_meaned_over_origins.png"),
  gg_coef_species,
  width = 9, height = 6, dpi = 300, bg = "white"
)
