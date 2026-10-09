library(tools)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Model_results/, etc.). Called at the end of 10_SDM_Code_Model_final.R.
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Paths
auc_path <- "Model_results/Metrics/AUC/"
tss_path <- "Model_results/Metrics/TSS/"

# --- AUC ---
# Get all mean CSV files
auc_files <- list.files(auc_path, pattern = "mean.*\\.csv$", full.names = TRUE)
# Extract filenames without extension
auc_names <- file_path_sans_ext(basename(auc_files))
# Split into parts
auc_parts <- strsplit(auc_names, "_")

# Build metadata for AUC
auc_meta <- data.frame(
  file    = auc_files,
  extent  = sapply(auc_parts, `[`, 3),
  origin  = sapply(auc_parts, `[`, 4),
  species = sapply(auc_parts, `[`, 5),
  stringsAsFactors = FALSE
)

# Read values and round
auc_values <- sapply(auc_meta$file, function(f) {
  df <- read.csv2(f, dec = ",")
  round(as.numeric(df$x[1]), 3)
})

auc_meta$AUC_mean <- auc_values

# --- TSS ---
# Get all mean CSV files
tss_files <- list.files(tss_path, pattern = "mean.*\\.csv$", full.names = TRUE)
tss_names <- file_path_sans_ext(basename(tss_files))
tss_parts <- strsplit(tss_names, "_")

# Build metadata for TSS
tss_meta <- data.frame(
  file    = tss_files,
  extent  = sapply(tss_parts, `[`, 3),
  origin  = sapply(tss_parts, `[`, 4),
  species = sapply(tss_parts, `[`, 5),
  stringsAsFactors = FALSE
)

# Read values and round
tss_values <- sapply(tss_meta$file, function(f) {
  df <- read.csv2(f, dec = ",")
  round(as.numeric(df$x[1]), 3)
})

tss_meta$TSS_mean <- tss_values

# --- Combine AUC and TSS ---
# Merge by extent, origin, species
result <- merge(
  auc_meta[, c("extent", "origin", "species", "AUC_mean")],
  tss_meta[, c("extent", "origin", "species", "TSS_mean")],
  by = c("extent", "origin", "species")
)

print(result)

# --- Export ---
output_file <- "Model_results/Metrics/metrics_summary.txt"
write.table(result, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

