library(ggcorrplot)
library(writexl)
library(terra)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Predictors/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

# Input:  Predictors/CHELSA/bio<1-19>.tif (03_CHELSA.R), Predictors/*_hesse.tif (04-07)
# Output: Plots/Correlationplots/ (correlation plots and tables of |r| >= 0.7)

dir.create("Plots/Correlationplots/Corr_table", recursive = TRUE, showWarnings = FALSE)

spearman_fun <- function(x, y) cor(x, y, method = "spearman", use = "pairwise.complete.obs")

#global
##########################################################################
results_list <- list()
bioclim_list <- 1:19

for (bioclim in bioclim_list) {

  bioclim_path <- paste0("Predictors/CHELSA/bio",bioclim,".tif")
  raster_layer <- rast(bioclim_path)
  names(raster_layer) <- paste0("bio", bioclim)
  results_list[[bioclim]] <- raster_layer
}

predictors_global <- rast(results_list)

valid_cells <- which(!is.na(values(predictors_global[[1]])))

# number of repetitions
n_reps <- 10

# initialize list to store correlation matrices
corr_matrices <- vector("list", n_reps)

# loop over repetitions
for (i in 1:n_reps) {
  message("Running i: ", i)
  set.seed(i)
  sampled_cells <- sample(valid_cells, 1000000)
  sampled_values <- terra::extract(predictors_global, sampled_cells)

  corr_res <- cor(sampled_values, method = "spearman", use = "pairwise.complete.obs")
  corr_matrices[[i]] <- as.matrix(corr_res)
}

mean_corr_matrix <- Reduce("+", corr_matrices) / n_reps

cor_long <- as.data.frame(as.table(mean_corr_matrix))
cor_long <- cor_long[cor_long$Var1 != cor_long$Var2, ]
cor_filtered <- cor_long[abs(cor_long$Freq) >= 0.7, ]
cor_filtered_unique <- cor_filtered[!duplicated(t(apply(cor_filtered[,1:2], 1, sort))), ]
write_xlsx(cor_filtered_unique, paste0("Plots/Correlationplots/Corr_table/corr_table_global.xlsx"))

korr_plot_global <- ggcorrplot(
  mean_corr_matrix,
  hc.order   = TRUE,
  outline.col= "black",
  ggtheme     = theme_gray(),
  colors      = c("#6D9EC1", "white", "#E46726"),
  lab         = TRUE,
  legend.title= "",
  type        = "lower"
)

korr_plot_global <- korr_plot_global +
  aes(alpha = abs(value) >= 0.7) +
  scale_alpha_manual(
    values = c("TRUE"  = 1,
               "FALSE" = 0.5),
    guide = "none"
  ); plot(korr_plot_global)

ggsave(paste0("Plots/Correlationplots/correlation_global.jpeg"), plot = korr_plot_global, width = 15, height = 15)



#local
##########################################################################
main_directory <- file.path(getwd(),"Predictors/")
raster_files <- list.files(path = main_directory, pattern = "_hesse.*\\.tif$", recursive = TRUE, full.names = TRUE)
predictors_local <- rast(raster_files)
raster_names <- tools::file_path_sans_ext(basename(raster_files))
names(predictors_local) <- raster_names
values(predictors_local$whc_hesse) <- as.numeric(as.factor(values(predictors_local$whc_hesse)))
values(predictors_local$heat_index_hesse) <- as.numeric(as.factor(values(predictors_local$heat_index_hesse)))

corr_res_local <- layerCor(
  predictors_local,
  fun = spearman_fun
)

corr_mat <- as.matrix(corr_res_local)

cor_long <- as.data.frame(as.table(corr_mat))
cor_long <- cor_long[cor_long$Var1 != cor_long$Var2, ]
cor_filtered <- cor_long[abs(cor_long$Freq) >= 0.7, ]
cor_filtered_unique <- cor_filtered[!duplicated(t(apply(cor_filtered[,1:2], 1, sort))), ]
write_xlsx(cor_filtered_unique, paste0("Plots/Correlationplots/Corr_table/corr_table_local.xlsx"))

layer_names <- names(predictors_local)
colnames(corr_mat) <- layer_names
rownames(corr_mat) <- layer_names

korr_plot <- ggcorrplot(
  corr_mat,
  hc.order   = TRUE,
  outline.col= "black",
  ggtheme     = theme_gray(),
  colors      = c("#6D9EC1", "white", "#E46726"),
  lab         = TRUE,
  legend.title= "",
  type        = "lower"
)

korr_plot <- korr_plot +
  aes(alpha = abs(value) >= 0.7) +
  scale_alpha_manual(
    values = c("TRUE"  = 1,
               "FALSE" = 0.5),
    guide = "none"
  ); plot(korr_plot)

ggsave(paste0("Plots/Correlationplots/correlation_local.jpeg"), plot = korr_plot, width = 15, height = 15)
