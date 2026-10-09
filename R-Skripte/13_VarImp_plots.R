# ============================
# Libraries
# ============================
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)
library(forcats)
library(patchwork)  # for panel layouts

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Model_results/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  Model_results/Feature_importance/fi__local_<origin>_<species>_mean.csv (10_SDM_Code_Model_final.R)
# Output: Model_results/Feature_importance/plots/AllSpecies_Averaged_FI_Panel.png (Figure 2);
#         printed category shares of importance (Results 3.3, "> 85 %")

theme_set(theme_minimal(base_size = 12) +
            theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white", colour = NA),
              plot.background  = element_rect(fill = "white", colour = NA)
            ))

# ============================
# Paths & output folders
# ============================
# Base path WITHOUT trailing slash
path_fi <- "Model_results/Feature_importance"

# Output folders
plots_dir   <- file.path(path_fi, "plots")
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)

# ============================
# Species & origin maps
# ============================
species_names <- c(
  "2891770" = "I. glandulifera",
  "3190653" = "A. altissima",
  "3034825" = "H. mantegazzianum"
)
origin_map  <- c("combined"="COM","gbif"="CS","hlnug"="StAg")
orig_levels <- c("COM","CS","StAg")

# ============================
# Category mapping
# ============================
map_category <- function(v) {
  vl <- tolower(v)
  if (vl %in% c("dgm","twi","aspect")) return("Topography")
  if (grepl("^floodplains", vl) || grepl("flowing_water_distance|standing_water_distance", vl)) return("Hydrology")
  if (grepl("^whc", vl)) return("Soil (WHC)")
  if (grepl("^corine", vl)) return("Land Use (CORINE)")
  if (grepl("roads_distance|railroads_distance|railroads_density|urban_builup_density|travel_time", vl)) return("Infrastructure")
  if (grepl("^heat_index", vl)) return("Climate (Heat Index)")
  return(NA_character_)
}

# ============================
# Pretty variable labels
# ============================
pretty_var <- function(v) {
  # CORINE abbreviations
  if (grepl("^CORINE_X", v)) {
    code <- sub("^CORINE_X", "", v)
    return(dplyr::case_when(
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
    ))
  }
  # Floodplains: hyphen, no parentheses
  if (grepl("^floodplains_X", v)) {
    k <- sub("^floodplains_X", "", v)
    level <- dplyr::case_when(
      k == "4" ~ "high",
      k == "3" ~ "medium",
      k == "2" ~ "low",
      k == "1" ~ "none",
      TRUE ~ k
    )
    return(paste0("Floodplains - ", level))
  }
  # Climate: class only
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
    return(eng)
  }
  # WHC: drop 'whc' and parentheses; keep class label only
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
      grepl("Siedlung und Verkehr|Siedlung.und.Verkehr", v, ignore.case = TRUE) ~ "Settlement and transport",
      TRUE ~ stringr::str_to_title(cls)
    )
    return(eng)
  }
  # Infra/Hydro: underscores → spaces; sentence case
  if (v %in% c("flowing_water_distance","standing_water_distance",
               "roads_distance","railroads_distance","railroads_density",
               "urban_builup_density","travel_time")) {
    s <- gsub("_", " ", v)
    return(stringr::str_to_sentence(s))
  }
  # Topography: nicer case
  if (v %in% c("DGM","twi","aspect")) {
    return(dplyr::case_when(
      v == "DGM" ~ "DGM",
      v == "twi" ~ "TWI",
      v == "aspect" ~ "Aspect",
      TRUE ~ v
    ))
  }
  gsub("_"," ", v)
}

# ============================
# Read FI CSVs (semicolon + comma decimals)
# ============================
fi_files <- list.files(path_fi, pattern="^fi__local_.*_\\d+_mean\\.csv$", full.names=TRUE)
stopifnot(length(fi_files) > 0)

read_fi <- function(fp){
  bn  <- basename(fp)
  origin  <- stringr::str_match(bn, "^fi__local_([^_]+)_\\d+_mean\\.csv$")[,2]
  species <- stringr::str_match(bn, "_(\\d+)_mean\\.csv$")[,2]
  dat <- readr::read_delim(fp, delim=";", locale=readr::locale(decimal_mark=","), show_col_types=FALSE, trim_ws=TRUE)
  cols <- tolower(names(dat))
  var_col  <- names(dat)[which(cols=="variable")[1]]
  mean_col <- names(dat)[which(cols %in% c("mean_overall","mean","fi","importance","meanimportance"))[1]]
  if (is.na(var_col) || is.na(mean_col)) stop("Bad columns in: ", bn)
  dat %>%
    dplyr::transmute(variable=.data[[var_col]], fi=as.numeric(.data[[mean_col]]),
                     species=species, origin=origin)
}

fi_tbl <- purrr::map_dfr(fi_files, read_fi) %>%
  dplyr::mutate(
    origin        = origin_map[origin],
    species_label = species_names[species],
    category      = vapply(variable, map_category, character(1)),
    var_pretty    = vapply(variable, pretty_var,   character(1))
  ) %>%
  dplyr::filter(!is.na(category))

# ============================
# Palettes (3 tones: COM, CS, StAg)
# ============================
color_palettes <- list(
  "Topography"            = c("#e9d8a6", "#ca9c68", "#8c510a"),
  "Hydrology"             = c("#67a9cf", "#2166ac", "#053061"),
  "Climate (Heat index)"  = c("#fcae91", "#fb6a4a", "#cb181d"),
  "Soil (WHC)"            = c("#fff7bc", "#fec44f", "#d95f0e"),
  "Land Use (CORINE)"     = c("#d8b9f2", "#9e77c2", "#6a3d9a"),
  "Infrastructure"        = c("#b2dfdb", "#1abc9c", "#00695c")
)
orig_levels <- c("COM","CS","StAg")

get_palette <- function(cat, palettes) {
  pal <- palettes[[cat]]
  if (!is.null(pal)) return(setNames(pal, orig_levels))
  keys <- names(palettes)
  hit  <- match(tolower(cat), tolower(keys))
  if (!is.na(hit)) return(setNames(palettes[[keys[hit]]], orig_levels))
  setNames(c("#bbbbbb","#888888","#555555"), orig_levels) # fallback
}

# ============================
# Global y-axis across ALL plots
# ============================
fi_min <- min(fi_tbl$fi, na.rm = TRUE)
fi_max <- max(fi_tbl$fi, na.rm = TRUE)

# ============================
# Category order
# ============================
cat_order <- c("Topography","Hydrology","Infrastructure",
               "Land Use (CORINE)","Soil (WHC)","Climate (Heat Index)")

# ============================
# Panel helper (legend top-right, y-title conditional)
# ============================
build_cat_plot <- function(df_cat, cat, show_y_title = FALSE) {
  pal_named <- get_palette(cat, color_palettes)
  xcol <- if ("var_pretty" %in% names(df_cat)) "var_pretty" else "variable"
  
  if (nrow(df_cat) == 0 || all(is.na(df_cat$fi))) return(NULL)
  
  df_cat <- df_cat %>%
    dplyr::mutate(origin = factor(origin, levels = orig_levels)) %>%
    tidyr::complete(!!rlang::sym(xcol), origin, fill = list(fi = NA_real_)) %>%
    dplyr::group_by(!!rlang::sym(xcol)) %>%
    dplyr::mutate(order_stat = mean(fi, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x_fac = forcats::fct_reorder(!!rlang::sym(xcol), order_stat, .desc = TRUE))
  
  ggplot(df_cat, aes(x = x_fac, y = fi, fill = origin)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, na.rm = TRUE) +
    scale_fill_manual(
      values = pal_named,
      breaks = orig_levels,
      labels = c(COM="Combined", CS="Citizen Science", StAg="State Agency"),
      name   = "Origin"
    ) +
    scale_y_continuous(limits = c(fi_min, fi_max)) +
    labs(
      title = cat,
      x = NULL,
      y = if (show_y_title) "Mean Predictor importance" else NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position      = c(1, 1),
      legend.justification = c(1, 1),
      legend.background    = element_rect(fill = "white", colour = NA),
      legend.direction     = "vertical",
      legend.box           = "vertical",
      axis.text.x          = element_text(angle = 45, hjust = 1),
      panel.grid.major     = element_blank(),
      panel.grid.minor     = element_blank(),
      plot.title           = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

# ============================
# NEW: Panel with ALL species averaged per variable × origin
# → 6 category panels in one figure, legend stays top-right per panel
# ============================

# 1) Average over species for each variable × origin × category
fi_all_mean <- fi_tbl %>%
  dplyr::group_by(category, var_pretty, origin) %>%
  dplyr::summarise(fi = mean(fi, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(!is.na(category), !is.na(var_pretty))

#mean per category
mean_per_category <- fi_all_mean %>%
  group_by(category) %>%
  summarise(mean_fi = mean(fi, na.rm = TRUE)) %>%
  arrange(desc(mean_fi))

mean_per_category

top5_per_origin <- fi_all_mean %>%
  group_by(origin) %>%
  arrange(desc(fi)) %>%
  slice_head(n = 5)

top5_per_origin


#mean per category per origin
mean_per_category_pct_all <- fi_all_mean %>%
  group_by(origin, category) %>%
  summarise(mean_fi = mean(fi, na.rm = TRUE), .groups = "drop") %>%
  group_by(origin) %>%
  mutate(
    total_fi = sum(mean_fi),
    pct = mean_fi / total_fi * 100
  )

main_categories <- c("Hydrology", "Infrastructure", "Topography")

mean_per_category_origin_main_pct <- mean_per_category_pct_all %>%
  filter(category %in% main_categories) %>%
  arrange(origin, desc(pct))

mean_per_category_origin_main_pct

# sum per category per origin
sum_per_category_pct_all <- fi_all_mean %>%
  group_by(origin, category) %>%
  summarise(sum_fi = sum(fi, na.rm = TRUE), .groups = "drop") %>%
  group_by(origin) %>%
  mutate(
    total_fi = sum(sum_fi),
    pct = sum_fi / total_fi * 100
  )

# extract only the main categories
main_categories <- c("Hydrology", "Infrastructure", "Topography")

sum_per_category_origin_main_pct <- sum_per_category_pct_all %>%
  filter(category %in% main_categories) %>%
  arrange(origin, desc(pct))

sum_per_category_origin_main_pct

# Global y-range for the averaged data
fi_min_all <- min(fi_all_mean$fi, na.rm = TRUE)
fi_max_all <- max(fi_all_mean$fi, na.rm = TRUE)

# Keep consistent category order
fi_all_mean <- fi_all_mean %>%
  dplyr::mutate(category = factor(category, levels = cat_order))

# 2) Build per-category plots (reuse helper with legend in top-right)
plots_raw_all <- vector("list", length(cat_order))
for (i in seq_along(cat_order)) {
  cat_i <- cat_order[i]
  df_cat <- fi_all_mean %>% dplyr::filter(tolower(category) == tolower(cat_i))
  
  # Show y-axis title only on left column (like before)
  show_y <- ((i - 1) %% 3 == 0)
  
  p_i <- build_cat_plot(df_cat, as.character(cat_i), show_y_title = show_y)
  
  if (!is.null(p_i)) {
    p_i <- p_i + scale_y_continuous(limits = c(fi_min_all, fi_max_all))
  }
  plots_raw_all[[i]] <- p_i
}

# Drop NULLs (categories with no data)
plots_all <- plots_raw_all[!vapply(plots_raw_all, is.null, logical(1))]

if (length(plots_all) > 0) {
  layout_all <- patchwork::wrap_plots(plots_all, ncol = 3, guides = "keep") +
    patchwork::plot_annotation(
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
    )
  
  out_all <- file.path(plots_dir, "AllSpecies_Averaged_FI_Panel.png")
  ggsave(
    filename = out_all,
    plot = layout_all,
    width = 14,
    height = 9,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  message("✅ Saved ALL-SPECIES averaged panel → ", out_all)
} else {
  message("⚠️ No data to plot for ALL-SPECIES averaged panel.")
}
