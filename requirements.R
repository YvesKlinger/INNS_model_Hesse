# requirements.R
#
# Installs the R packages used across all scripts in R-Skripte/.

cran_packages <- c(
  "rgbif",           # GBIF occurrence download (01_SDM_Code_GBIF_download.R)
  "dismo",           # SDM utilities, gridSample (02_SDM_Code_Occurrence_filter.R)
  "raster",          # legacy raster handling (02_SDM_Code_Occurrence_filter.R, 03_CHELSA.R)
  "terra",           # raster handling (all scripts)
  "sf",              # vector/spatial data
  "sp",              # legacy spatial classes
  "readr",           # fast text file reading
  "data.table",      # fast large-table reading (fread)
  "dplyr",           # data wrangling
  "tidyr",           # data wrangling (pivot_longer etc.)
  "tibble",          # rownames_to_column etc.
  "CoordinateCleaner", # occurrence coordinate cleaning
  "ggplot2",         # plotting
  "ggrepel",         # non-overlapping plot labels
  "patchwork",       # combining ggplot panels
  "ecospat",         # niche overlap (Schoener's D, equivalency test)
  "caret",           # model training framework
  "gbm",             # gradient boosting machine (SDM algorithm used)
  "pROC",            # ROC/AUC calculation
  "CAST",            # spatial CV, forward feature selection, AOA
  "blockCV",         # spatial block cross-validation
  "gstat",           # variogram fitting (blockCV block size)
  "recipes",         # preprocessing (scaling, dummy coding)
  "broom",           # tidy model summaries (tidy(lm(...)))
  "glmmTMB",         # beta GLMMs (12_Model_performance_comparison.R)
  "emmeans",         # marginal means / pairwise contrasts (12_Model_performance_comparison.R)
  "DHARMa",          # simulated residual diagnostics (12_Model_performance_comparison.R)
  "flextable",       # supplementary Word table (12_Model_performance_comparison.R)
  "officer",         # supplementary Word table (12_Model_performance_comparison.R)
  "whitebox",        # aspect, slope, TWI (04_DGM_und_abhaengige_Variablen.R); run whitebox::install_whitebox() once
  "ggcorrplot",      # correlation plots (08_correlation_matrix_complete_model_area.R)
  "writexl",         # correlation tables (08_correlation_matrix_complete_model_area.R)
  "stringr",         # string handling (plot scripts 13-17)
  "purrr",           # iteration (plot scripts 13-16)
  "forcats",         # factor ordering (13_VarImp_plots.R)
  "stars",           # raster plotting (17_Suitability_maps.R)
  "cowplot"          # figure layout (17_Suitability_maps.R)
)

installed <- rownames(installed.packages())
to_install <- setdiff(cran_packages, installed)
if (length(to_install) > 0) {
  install.packages(to_install)
}
