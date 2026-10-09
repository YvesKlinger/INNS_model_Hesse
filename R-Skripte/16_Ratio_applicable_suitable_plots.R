library(terra)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Model_results/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  Model_results/Maps/binary/predicted_raster_binary__local_<origin>_<species>_<run>.tif and
#         Model_results/AOA/aoa_binary_local_<origin>_<species>_mean.tif (10_SDM_Code_Model_final.R)
# Output: Model_results/Suitability_ratio_plot.{png,pdf} (Figure 4: a) suitable share of each
#         model's AOA, b) extent of the AOA), Model_results/Suitability_AOA_data.csv,
#         Model_results/Maps/binary/binary__local_<origin>_<species>_vote8of10.tif

#Suitable - non suitable cell count
#####################################################################################
#####################################################################################

# --- 1) Files (single-iteration binaries only; excludes *_mean/_sd) ---
pred_files <- list.files(
  path = "Model_results/Maps/binary/",
  pattern = "binary__local.*_[0-9]+\\.tif$",
  full.names = TRUE
)

aoa_files <- list.files(
  path = "Model_results/AOA/",
  pattern = "aoa_binary.*_mean\\.tif$",
  full.names = TRUE
)

# --- 2) Parse species + origin from filenames ---
pred_df <- tibble(file = pred_files) %>%
  mutate(
    species = str_extract(file, "\\d{7}"),
    origin  = str_extract(file, "(?<=binary__).*?(?=_[0-9]{7})")
  )

# --- 3) Majority-vote consensus per species × origin (>= 8/10) ---
vote_threshold <- 0.8  # 8 of 10

results <- pred_df %>%
  group_by(species, origin) %>%
  group_split() %>%
  map_dfr(function(g) {
    # stack the binary maps for this species-origin
    rlist  <- lapply(g$file, rast)
    rstack <- rast(rlist)
    agree  <- mean(rstack, na.rm = TRUE)  # per-pixel agreement proportion (0..1)

    # AOA mask for this species-origin (if present)
    aoa_file <- aoa_files[grepl(g$species[1], aoa_files) & grepl(g$origin[1], aoa_files)]
    if (length(aoa_file) == 1) {
      aoa_r <- rast(aoa_file)
      aoa_r <- resample(aoa_r, agree, method = "near")
      agree <- mask(agree, aoa_r, maskvalues = 0)  # keep only inside AOA
    }

    # consensus binary by majority vote
    consensus_bin <- ifel(agree >= vote_threshold, 1, 0)

    # (optional) save consensus raster
    out_dir <- "Model_results/Maps/binary/"
    out_fn  <- file.path(out_dir, paste0("binary__", g$origin[1], "_", g$species[1], "_vote8of10.tif"))
    writeRaster(consensus_bin, out_fn, overwrite = TRUE)

    # counts for 0/1 inside AOA
    counts <- freq(consensus_bin, digits = 0)
    tibble(
      species = g$species[1],
      origin  = g$origin[1],
      value   = counts$value,
      count   = counts$count
    )
  })

# --- 4) Build table with percentages (suitable = value==1) ---
pred_table <- results %>%
  group_by(species, origin) %>%
  mutate(percent = round(count / sum(count) * 100, 1)) %>%
  ungroup()

pred_suitable <- pred_table %>% filter(value == 1)

# --- 5) Extent of each model's AOA (denominator of the percentages) ---
aoa_table <- tibble(file = aoa_files) %>%
  mutate(species = str_extract(file, "\\d{7}"),
         origin  = str_extract(file, "local_[a-z]+")) %>%
  mutate(aoa_pct = map_dbl(file, function(f) {
    v <- values(rast(f), mat = FALSE)
    round(100 * sum(v == 1, na.rm = TRUE) / sum(!is.na(v)), 1)
  })) %>%
  select(-file)

species_names <- c(
  "3190653" = "A. altissima",
  "3034825" = "H. mantegazzianum",
  "2891770" = "I. glandulifera"
)
source_labels <- c(
  "local_hlnug"    = "State Agency",
  "local_gbif"     = "Citizen Science",
  "local_combined" = "Combined"
)
source_colours <- c("State Agency" = "#9467BD", "Citizen Science" = "#59A14F", "Combined" = "#F28E2B")

fig_data <- pred_suitable %>%
  select(species, origin, suitable_pct = percent) %>%
  left_join(aoa_table, by = c("species", "origin")) %>%
  mutate(species = factor(species_names[species], levels = species_names),
         source  = factor(source_labels[origin], levels = rev(source_labels))) %>%
  arrange(species, desc(source))
write.csv(select(fig_data, species, source, suitable_pct, aoa_pct),
          "Model_results/Suitability_AOA_data.csv", row.names = FALSE)
print(fig_data)

# --- 6) Plot: grouped bars, percentages relative to each model's own AOA ---
panel_plot <- function(var, x_label, title, x_max, x_break, show_y = TRUE) {
  p <- ggplot(fig_data, aes(x = .data[[var]], y = source, fill = source)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = sprintf("%.1f%%", .data[[var]])), hjust = -0.15, size = 3.2) +
    facet_grid(species ~ ., switch = "y") +
    scale_fill_manual(values = source_colours, guide = "none") +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, 100, by = x_break),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = x_label, y = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1, size = 11),
      panel.spacing.y = unit(0.8, "lines"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = 11),
      plot.title.position = "plot"
    )
  if (!show_y) {
    p <- p + theme(strip.text.y.left = element_blank(), axis.text.y = element_blank())
  }
  p
}

p_suitable <- panel_plot("suitable_pct", "Cells predicted as suitable (% of the model's AOA)",
                         "a) Suitable proportion of the AOA", x_max = 50, x_break = 10)
p_aoa <- panel_plot("aoa_pct", "Applicable area (% of Hesse)",
                    "b) Extent of the AOA", x_max = 110, x_break = 25, show_y = FALSE)

p_fig <- p_suitable + p_aoa + plot_layout(widths = c(1, 1))

plot(p_fig)

ggsave("Model_results/Suitability_ratio_plot.png", p_fig, width = 22, height = 11, units = "cm", dpi = 600, bg = "white")
ggsave("Model_results/Suitability_ratio_plot.pdf", p_fig, width = 22, height = 11, units = "cm")
