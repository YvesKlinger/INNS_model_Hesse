#Setup
##########################################################################
##########################################################################

#Packages
library(caret)
library(terra) 
library(sp)
library(dplyr)
library(gbm)
library(pROC)
library(CAST)
library(gstat)
library(sf)
library(ggplot2)
library(recipes)
library(tibble)
library(blockCV)
library(tidyr)
library(broom)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Occurrences/, Predictors/, Model_results/, etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")


#load predictors
##########################################################################
##########################################################################

# load global predictors
##########################################################################
results_list <- list()
bioclim_list <- c(1, 4, 15, 16, 17, 18)

for (bioclim in bioclim_list) {
  bioclim_path <- paste0("Predictors/CHELSA/bio", bioclim, ".tif")
  raster_layer <- rast(bioclim_path)
  name <- paste0("bio", bioclim)
  names(raster_layer) <- name
  results_list[[name]] <- raster_layer
}
predictors_global <- rast(results_list)

# load local predictors
##########################################################################
main_directory <- file.path(getwd(),"Predictors/")
raster_files <- list.files(path = main_directory, pattern = "_hesse.*\\.tif$", recursive = TRUE, full.names = TRUE)
raster_stack <- rast(raster_files)
raster_names <- tools::file_path_sans_ext(basename(raster_files))
raster_names_clean <- gsub("_hesse$", "", raster_names)
names(raster_stack) <- raster_names_clean
predictors_local <- raster_stack[[!duplicated(names(raster_stack))]] %>% na.omit()

#stack species-specific
predictors_local_2891770 <- predictors_local[[names(predictors_local) != "heat_index"]]
predictors_local_3034825 <- predictors_local[[names(predictors_local) != "heat_index"]]
predictors_local_3190653 <- predictors_local[[names(predictors_local) != "whc"]]


#load occurrences
##########################################################################
##########################################################################

# load global occurrences
##########################################################################
occ_global_gbif <- read.csv2("Occurrences/Occurrence_data_GBIF/global/filtered/all_species_gbif_global.csv") %>% 
  mutate(presence = "present") %>%
  mutate(decimalLatitude = as.numeric(decimalLatitude),
    decimalLongitude = as.numeric(decimalLongitude))

# load local occurrences
##########################################################################
occ_local_gbif <- read.csv2("Occurrences/Occurrence_data_GBIF/local/filtered/all_species_gbif_local.csv") %>% 
  mutate(presence = "present") %>% 
  mutate(weight = 1) %>%
  mutate(decimalLatitude = as.numeric(decimalLatitude),
         decimalLongitude = as.numeric(decimalLongitude))

occ_local_hlnug <- read.csv2("Occurrences/Occurrence_data_HLNUG/filtered/all_species_hlnug.csv") %>% 
  mutate(presence = "present") %>% 
  mutate(weight = 1) %>%
  mutate(decimalLatitude = as.numeric(decimalLatitude),
         decimalLongitude = as.numeric(decimalLongitude))

occ_local_combined <- read.csv2("Occurrences/Occurrence_data_combined/filtered/all_species_combined.csv") %>% 
  mutate(presence = "present") %>% 
  mutate(weight = 1) %>%
  mutate(decimalLatitude = as.numeric(decimalLatitude),
         decimalLongitude = as.numeric(decimalLongitude))

#load model extent
##########################################################################
##########################################################################

hesse_raster <- rast("Predictors/DGM_hesse.tif")

model_extent <- list("global","local")


# SDM
##########################################################################
##########################################################################

# set loop for extent (global, local), origin (CS, StAg, COM), 
#species (I. glandulifera, H. mantegazzianum, A. altissima)
##############################################################
for (extent in model_extent){
  message("Running extent: ", extent)
  if (extent == "global") {
    origin_list <- "gbif" 
  } else {
    origin_list <- c("gbif", "hlnug", "combined")
  }
  
  for (origin in origin_list) {
    message("Running data origin: ", origin)
    occ <- get(paste0("occ_",extent,"_",origin))
    occ[c("decimalLongitude", "decimalLatitude")] <- lapply(occ[c("decimalLongitude", "decimalLatitude")], as.numeric)
    species_list <- unique(occ$speciesKey)
    
    for (species in species_list) {
      message("Running species: ", species)
      
      #set up result lists for current extent_origin_species
      ##########################################################################
      pa_result_list <- list()
      final_occ_list <- list()
      best_model_list <- list()
      model_train_roc_list <- list()
      feature_importance_list <- list()
      response_table_list <- list()
      predicted_raster_list <- list()
      predicted_raster_mean <- list()
      predicted_binary_list <- list()
      confusion_matrix_list <- list()
      auc_list <- list()
      best_threshold_list <- list()
      tss_list <- list()
      auc_mean_list <- list()
      tss_mean_list <- list()
      aoa_result_list <- list()
      NA_LM_list <- list()
      aoa_pred_weight_list <- list()

      #load species data for current extent_origin_species
      species_data <- occ %>%
        filter(speciesKey == species) %>%
        dplyr::select(speciesKey, decimalLongitude, decimalLatitude, presence, any_of("weight")) %>%
        dplyr::select(-speciesKey)

      #get global or local predictors depending on extent
      if (extent == "global") {
        predictors <- get(paste0("predictors_",extent))
      } else {
        predictors <- get(paste0("predictors_",extent,"_",species))
      }
      
      #compute space for PA selection
      ##########################################################################
      #for extent = global pre-run PA raster < 20 km around occurrences
      if (extent == "global") {
        base_raster <- rast("Predictors/CHELSA/bio1.tif")
        pa_area <- rast(paste0("Pseudo_absences/global/pa_area/",species,"_ras.tif"))
        
      } else {
        #for extent = local compute PCA and calculate euclidian distance in PCA space
        #select threshold to only use "known" space for PA selection
        factor_cols <- c("CORINE", "floodplains", "whc", "heat_index")
        
        occ_pca <- species_data %>%
          dplyr::select(decimalLongitude, decimalLatitude)
        
        occ_env <- terra::extract(predictors, occ_pca) %>%
          dplyr::select(-ID)
        
        all_env <- terra::as.data.frame(predictors, na.rm = FALSE)
        
        occ_env <- occ_env %>% mutate(across(any_of(factor_cols), as.factor))
        all_env <- all_env %>% mutate(across(any_of(factor_cols), as.factor))
        
        rec_pca <- recipe(~ ., data = all_env) %>%
          step_unknown(all_nominal_predictors(), new_level = "NA") %>%
          step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
          prep()
        
        all_env_dummy <- bake(rec_pca, new_data = all_env)
        occ_env_dummy <- bake(rec_pca, new_data = occ_env)
        
        common_cols <- intersect(names(all_env_dummy), names(occ_env_dummy))
        all_env_dummy  <- all_env_dummy[, common_cols, drop = FALSE]
        occ_env_dummy  <- occ_env_dummy[, common_cols, drop = FALSE]
        
        var_cols <- sapply(all_env_dummy, function(x) var(x, na.rm = TRUE))
        keep_cols <- names(var_cols)[!is.na(var_cols) & var_cols > 0]
        all_env_dummy <- all_env_dummy[, keep_cols, drop = FALSE]
        occ_env_dummy <- occ_env_dummy[, keep_cols, drop = FALSE]
        
        valid_rows <- complete.cases(all_env_dummy)
        
        dat_for_pca <- all_env_dummy[valid_rows, , drop = FALSE]
        
        dat_for_pca <- data.frame(lapply(dat_for_pca, function(x) as.numeric(x)))
        
        col_var <- apply(dat_for_pca, 2, function(x) var(x, na.rm = TRUE))
        keep_cols <- names(col_var)[!is.na(col_var) & col_var > 0]
        
        if (length(keep_cols) < 2) stop("Too few variable columns for PCA")
        
        dat_for_pca <- dat_for_pca[, keep_cols, drop = FALSE]
        
        missing_in_occ <- setdiff(keep_cols, names(occ_env_dummy))
        if (length(missing_in_occ) > 0) occ_env_dummy[, missing_in_occ] <- 0
        occ_for_pca <- occ_env_dummy[, keep_cols, drop = FALSE]
        
        pca_all <- prcomp(dat_for_pca, center = TRUE, scale. = TRUE)

        occ_scores <- predict(pca_all, newdata = occ_env_dummy)

        occ_centroid <- colMeans(occ_scores, na.rm = TRUE)

        grid_scores <- pca_all$x

        env_dist_to_occ <- sqrt(rowSums((grid_scores - matrix(occ_centroid, nrow = nrow(grid_scores), ncol = ncol(grid_scores), byrow = TRUE))^2))

        occ_distances_to_centroid <- sqrt(rowSums((occ_scores - matrix(occ_centroid, nrow = nrow(occ_scores), ncol = ncol(occ_scores), byrow = TRUE))^2))
        # distance threshold
        if (species == "2891770") {
          threshold <- quantile(occ_distances_to_centroid, 0.95, na.rm = TRUE)
        } else if (species == "2891770") {
          threshold <- quantile(occ_distances_to_centroid, 0.75, na.rm = TRUE)
        } else {
          threshold <- quantile(occ_distances_to_centroid, 0.55, na.rm = TRUE)
        }
        
        similar_env_mask <- rep(NA_real_, nrow(all_env_dummy))
        similar_env_mask[valid_rows] <- ifelse(env_dist_to_occ <= threshold, 1, 0)
        
        mask_raster <- predictors_local_2891770[[1]]
        values(mask_raster) <- similar_env_mask
 
        #exclude distances from occs < 1 km and greater > 20 km
        species_data_sf <- st_as_sf(species_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
        buffer_inner <- st_transform(st_buffer(st_union(st_convex_hull(species_data_sf)), dist = 1000), crs = st_crs(mask_raster))
        buffer_outer <- st_transform(st_buffer(st_union(st_convex_hull(species_data_sf)), dist = 20000),crs = st_crs(mask_raster))
        
        r_df <- terra::as.data.frame(mask_raster, xy = TRUE, cells = TRUE)   # columns: x,y,cells,layer
        names(r_df)[which(names(r_df) == names(r_df)[4])] <- "mask_value"     # ensure clear name
  
        cand <- r_df %>% filter(!is.na(mask_value) & mask_value == 1)
      
        cand_sf <- st_as_sf(cand, coords = c("x", "y"), crs = terra::crs(mask_raster))
        cand_sf <- st_transform(cand_sf, st_crs(buffer_inner))
     
        inside_outer <- st_within(cand_sf, buffer_outer, sparse = FALSE)[,1]
        inside_inner <- st_within(cand_sf, buffer_inner, sparse = FALSE)[,1]
        
        cand_keep <- cand %>% filter(inside_outer & !inside_inner)
        
        # create final pa_area raster
        pa_area <- mask_raster
        values(pa_area) <- NA_real_
        pa_area[ cand_keep$cell] <- 1
        sum_cells <- sum(values(pa_area) == 1, na.rm = TRUE)

        write.table(sum_cells, paste0("model_results/pa_area_local/num_cells_pa_area_",origin,"_",species,".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        terra::writeRaster(pa_area, paste0("Model_results/pa_area_local/",origin,"_",species,".tif"), overwrite = TRUE)
      }
      
      #Number of runs (global = 5 and local = 10)
      n_runs <- if (extent == "global") 5 else 10
      
      #Start Loop for number of runs (different random PA absence selection)
      ##########################################################################
      for (i in 1:n_runs) {
        message("Running i: ", i)
        
        tryCatch({
        #sample PA absences
        set.seed(526 + i)
        NoP <- nrow(species_data)
        name_it <- paste0("_",extent,"_",origin,"_",species, "_", i)
        
        available_cells <- which(!is.na(values(pa_area)))
        sampled_cells <- sample(available_cells, size = NoP, replace = FALSE)

        sampled_coords <- xyFromCell(pa_area, sampled_cells)
        
        pseudoabs_pts <- vect(sampled_coords, crs = crs(pa_area), type = "points")
        
        pseudo_abs <- terra::as.data.frame(pseudoabs_pts, geom = "XY") %>% 
          rename(decimalLongitude = x, decimalLatitude = y) %>% 
          mutate(presence = "absent")
        
        pa_result_list[[name_it]] <- pseudo_abs
        write.csv2(pseudo_abs, file = paste0("Pseudo_absences/",extent,"/pas_",name_it,".csv"), row.names = FALSE)
        
        if (extent == "local"){
          # for extent = local weigh pseudo-absences depending on global model results
          ############################################################################
          weighing_raster <- rast(paste0("Model_results/Maps/_global_gbif_",species,"_mean.tif"))
          extracted_weight <- terra::extract(weighing_raster, pa_result_list[[name_it]][, c("decimalLongitude", "decimalLatitude")])
          
          pa_with_weight <- cbind(pa_result_list[[name_it]][, c("decimalLongitude", "decimalLatitude", "presence")], extracted_weight) %>%
            dplyr::select(-ID) %>%
            rename(weight = mean)
          pa_with_weight$weight <- 1 - pa_with_weight$weight
          
          pseudo_abs <- pa_with_weight
        }
        
        #Combine Occurences with PAs
        #########################################################################
        common_cols <- intersect(c("decimalLongitude", "decimalLatitude", "presence", "weight"),
                                 intersect(names(species_data), names(pseudo_abs)))
        
        occurrence_with_pa <- full_join(species_data, pseudo_abs, by = common_cols) %>% 
          dplyr::select(decimalLongitude, decimalLatitude, presence, any_of("weight"))
        
        final_occ_list[[name_it]] <- occurrence_with_pa %>%
          mutate(extracted_values = terra::extract(predictors, dplyr::select(., decimalLongitude, decimalLatitude))) %>%
          bind_cols(as.data.frame(.$extracted_values)) %>%
          dplyr::select(-any_of(c("extracted_values", "ID"))) 
        
        final_occ <- final_occ_list[[name_it]] %>%
          na.omit()
        if(any(c("CORINE", "floodplains", "whc", "heat_index") %in% colnames(final_occ))) {
          final_occ <- final_occ %>%
            mutate(across(any_of(c("CORINE", "floodplains", "whc", "heat_index")), as.factor)) 
        }
        #blockCV - test and train data
        ####################################################################
        if (extent == "local"){

          occ_data_sf <- st_as_sf(final_occ,
                                  coords    = c("decimalLongitude", "decimalLatitude"),
                                  crs       = 4326)
          
          predictors_utm <- terra::project(predictors, "EPSG:25832")
          
          occ_data_sf_utm <- sf::st_transform(occ_data_sf, crs = sf::st_crs(predictors_utm)) %>%
            mutate(presence = recode(presence, "present" = 1,
                                     "absent" = 0))
          
          vgram <- variogram(presence ~ 1, data = occ_data_sf_utm)
          model_vgram <- fit.variogram(vgram, model = vgm(psill = 0.5, model = "Exp", range = 3000, nugget = 0.1))
          
          size <- model_vgram[model_vgram$model != "Nug", "range"]
          
          cv_result <- cv_spatial(
            x         = occ_data_sf_utm,
            column    = "presence",
            r         = predictors_utm,
            size      = size,
            k         = 6,
            selection = "systematic",
            hexagon   = TRUE,
            plot      = FALSE
          )
          
          occ_data_sf_utm$folds <- cv_result$folds_ids
          
          test_fold <- sample(1:6, 1)
          
          test_indices  <- which(occ_data_sf_utm$folds == test_fold)
          train_indices <- which(occ_data_sf_utm$folds != test_fold)
          
          test_data_pre <- occ_data_sf_utm[test_indices, ]
          coords_test_data <- st_coordinates(test_data_pre) %>% as.data.frame() %>%
            rename(decimalLongitude = X, decimalLatitude = Y)
          
          test_data_pre <- test_data_pre %>% 
            dplyr::select(-folds) %>%
            st_transform(crs = 4326) %>%
            mutate(presence = recode(presence, "1" = "present", "0" = "absent")) %>%            
            st_drop_geometry()
          
          train_data_pre <- occ_data_sf_utm[train_indices, ] %>%
            st_transform(crs = 4326) %>%
            mutate(presence = recode(presence, "1" = "present", "0" = "absent")) %>%            
            st_drop_geometry()
          
          predictors_train_pre <- train_data_pre[, !names(train_data_pre) %in% c("decimalLongitude", "decimalLatitude", "presence")]
          predictors_test_pre <- test_data_pre[, !names(test_data_pre) %in% c("decimalLongitude", "decimalLatitude", "presence")]
          
        } else {
          #split in test and train data for global
          ##########################################################################
          train_index <- createDataPartition(final_occ$presence, p = 0.8, list = FALSE)
          train_data_pre <- final_occ[train_index, ] %>%
            mutate(presence = factor(presence, levels = c("absent", "present")))
          test_data_pre <- final_occ[-train_index, ] %>%
            mutate(presence = factor(presence, levels = c("absent", "present")))
          coords_test_data <- test_data_pre %>%
            dplyr::select(decimalLongitude, decimalLatitude)
          
          predictors_train_pre <- train_data_pre[, !names(train_data_pre) %in% c("decimalLongitude", "decimalLatitude", "presence")]
          predictors_test_pre <- test_data_pre[, !names(test_data_pre) %in% c("decimalLongitude", "decimalLatitude", "presence")]
        }
        
        #stardise train-set (z-score)
        ##########################################################################
        if (extent == "global"){
          rec_pre <- recipe(presence ~ ., data = train_data_pre) %>%
            step_normalize(all_numeric(), -all_outcomes(), -has_role("ID"),-any_of(c("decimalLongitude", "decimalLatitude")))
        } else {
          rec_pre <- recipe(presence ~ ., data = train_data_pre) %>%
            update_role(any_of(c("weight", "folds")), new_role = "ID") %>%
            step_dummy(all_nominal_predictors(), -all_outcomes()) %>%
            step_normalize(all_numeric(), -all_outcomes(), -has_role("ID"))
        }
        
        rec_prep_pre <- prep(rec_pre, training = train_data_pre)
        train_data_scaled_pre <- bake(rec_prep_pre, new_data = train_data_pre) %>%
          as.data.frame()
        predictors_train_scaled_pre <- train_data_scaled_pre[, !names(train_data_scaled_pre) %in% c("weight","folds", "presence")]
        
        #Set Control Parameters and CV
        ##########################################################################
        if (extent == "global"){
          ctrl <- trainControl(method = "cv", 
                               number = 5,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               savePredictions = "final")
          
        } else {
          index <- split(
            seq_len(nrow(train_data_scaled_pre)),
            train_data_scaled_pre$folds)
          
          index_df <- do.call(rbind, lapply(names(index), function(fold) {
            data.frame(Fold = fold, RowNumber = index[[fold]])
          }))
          fold_counts <- as.data.frame(table(index_df$Fold))
          write.csv(fold_counts, paste0("Model_results/folds_count/fold_counts_",name_it,".csv"), row.names = FALSE)
          
          ctrl <- trainControl(
            method      = "cv",
            index       = index,
            summaryFunction = twoClassSummary,
            classProbs     = TRUE,
            savePredictions = "final"
          )
        }
        
        #grid for tuning
        ##########################################################################
        tuneGrid <- expand.grid(
          n.trees = c(1000, 1500),
          interaction.depth = c(2,3),
          shrinkage = c(0.01),
          n.minobsinnode = c(5, 10))
        
        #Predictor selection
        ##########################################################################
        predictors_train_ffs <- predictors_train_pre %>%
          dplyr::select(-any_of(c("folds", "weight")))
        
        ffs_result <- ffs(predictors_train_ffs,
                          train_data_scaled_pre$presence,
                          method = "gbm",
                          trControl = ctrl,
                          metric = "ROC",
                          maximize = TRUE,
                          tuneGrid = tuneGrid)

        chosen_predictors <- ffs_result[["finalModel"]][["xNames"]]
        predictors_ffs <- predictors[[chosen_predictors]]
        
        #Filter only selected predictors from test and train data
        ##########################################################################
        train_data <- train_data_pre %>%
          dplyr::select(presence, any_of(c("weight","folds")), all_of(chosen_predictors)) 
        
        predictors_train <- predictors_train_pre %>%
          dplyr::select(all_of(chosen_predictors))
        
        test_data <- test_data_pre %>%
          dplyr::select(presence, all_of(chosen_predictors))
        
        predictors_test <- predictors_test_pre %>%
          dplyr::select(all_of(chosen_predictors))
        
        if (extent == "local"){
          #weight as vector for model
          ##########################################################################
          weight_vector <- train_data$weight
          train_data <- train_data %>% dplyr::select(-any_of(c("weight")))
        }
        
        #stardise again for the new predictor choice of the train-set (z-score)
        ##########################################################################
        #########################################################################
        if (extent == "global"){
          rec <- recipe(presence ~ ., data = train_data) %>%
            step_dummy(all_nominal_predictors(), -all_outcomes()) %>%
            step_normalize(all_numeric(), -all_outcomes())
          
          rec_prep <- prep(rec, training = train_data)
          train_data_scaled <- bake(rec_prep, new_data = train_data) %>% as.data.frame()
        } else {
          rec <- recipe(presence ~ ., data = train_data) %>%
            update_role(folds, new_role = "id") %>% 
            step_dummy(all_nominal_predictors(), -all_outcomes(), -all_of("folds")) %>%
            step_normalize(all_numeric(), -all_outcomes(), -all_of("folds")) 
          rec_prep <- prep(rec, training = train_data)
          train_data_scaled <- bake(rec_prep, new_data = train_data) %>% as.data.frame()
        
        #remove columns with zero variance in one of the blocks
        numeric_cols <- names(train_data_scaled)[sapply(train_data_scaled, is.numeric)]
        exclude_cols <- c("presence", "folds")
        num_predictors <- setdiff(numeric_cols, exclude_cols)
        
        fold_var <- train_data_scaled %>%
          group_by(folds) %>%
          summarise(across(all_of(num_predictors), ~ var(.x, na.rm = TRUE)), .groups = "drop")
        
        zero_in_any_fold <- sapply(fold_var[num_predictors], function(col) {
          any(is.na(col) | col == 0)
        })
        
        drop_cols <- names(zero_in_any_fold)[zero_in_any_fold]
            message("Removing predictors that have zero variance in at least one fold:\n",
            paste(drop_cols, collapse = ", "))
        train_data_scaled <- train_data_scaled %>% dplyr::select(-all_of(c(drop_cols,"folds")))
        }
        
        #training model
        ##########################################################################
        args <- list(
          presence ~ .,
          data = train_data_scaled, 
          method = "gbm",
          trControl = ctrl,
          tuneGrid = tuneGrid,
          metric = "ROC"
        )
        
        if(extent == "local") {
          args$weights <- weight_vector
        }
        
        gbm_model <- do.call(train, args)
        
        model_best_params <- gbm_model$bestTune
        
        #train best model
        ##########################################################################
        args2 <- list(
          presence ~ .,
          data = train_data_scaled,
          method = "gbm",
          trControl = ctrl,
          tuneGrid = model_best_params
        )
        
        if(extent == "local") {
          args2$weights <- weight_vector
        }
        
        best_model_list[[name_it]] <- do.call(train, args2)
        best_model <- best_model_list[[name_it]]
        saveRDS(best_model, file = paste0("Model_results/best_models/best_model_", name_it, ".rds"))
        
        model_train_roc_list[[name_it]] <- best_model_list[[name_it]][["results"]][["ROC"]]
        
        #feature importance
        ##########################################################################
        importance <- varImp(best_model)
        importance_df <- as.data.frame(importance$importance)
        importance_df <- rownames_to_column(importance_df , var = "variable")
        
        feature_importance_list[[name_it]] <- importance_df
        write.csv2(importance_df, file = paste0("Model_results/Feature_importance/fi_",name_it,".csv"), row.names = FALSE)
        
        #data for response curves
        ########################################################################
        norm_step <- rec_prep$steps[[2]]
        means <- norm_step$means
        sds   <- norm_step$sds
        
        cols_to_unscale <- intersect(names(means), names(train_data_scaled))
        
        train_data_unscaled <- train_data_scaled %>%
          dplyr::mutate(
            dplyr::across(
              .cols = all_of(cols_to_unscale),
              .fns  = ~ .x * sds[[cur_column()]] + means[[cur_column()]]
            ))
    
        train_data_unscaled$prob[best_model$pred$rowIndex] <- best_model$pred$present
        response_table_list[[name_it]] <- train_data_unscaled
        write.csv2(train_data_unscaled, file = paste0("Model_results/response_data/res_",name_it,".csv"), row.names = FALSE)
        
        #stardise test-data with train fit (z-score)
        ##########################################################################
        if (extent == "local") {
        test_data$folds <- 0 #needed for baking because in recipe 
        }
        test_data_scaled <- bake(rec_prep, new_data = test_data) %>% as.data.frame() %>%
          dplyr::select(-any_of(c("weight","folds")))
        
        #Test model (AUC & TSS)
        ##########################################################################
        probabilities <- predict(best_model, newdata = test_data_scaled, type = "prob")[, "present"]
        roc_curve <- roc(test_data_scaled$presence, probabilities)
        auc_value <- auc(roc_curve)
        auc_list[[name_it]] <- auc_value
        write.csv2(auc_value, file = paste0("Model_results/Metrics/AUC/auc_",name_it,".csv"), row.names = FALSE)
        
        best_threshold <- roc_curve$thresholds[which.max(roc_curve$sensitivities + roc_curve$specificities)]
        best_threshold_list[[name_it]] <- best_threshold
        write.csv2(best_threshold, file = paste0("Model_results/Metrics/threshold_binary/threshold_binary_",name_it,".csv"), row.names = FALSE)
        predicted_class <- ifelse(probabilities >= best_threshold, "present", "absent")
        predicted_class <- factor(predicted_class, levels = c("absent", "present"))
        observed_class <- factor(test_data$presence, levels = c("absent", "present"))
        conf_matrix <- confusionMatrix(predicted_class, observed_class)
        confusion_matrix_list[[name_it]] <- conf_matrix
        write.csv2(conf_matrix$table, file = paste0("Model_results/Metrics/confusion_matrix/cm_",name_it,".csv"), row.names = FALSE)
        
        tss <- conf_matrix$byClass["Sensitivity"] + conf_matrix$byClass["Specificity"] - 1
        tss_list[[name_it]] <- tss
        write.csv2(tss, file = paste0("Model_results/Metrics/TSS/tss_",name_it,".csv"), row.names = FALSE)
        
        #Predict model to whole raster
        ##########################################################################
        
        gc() #free storage for raster calc
        
        final_predictors <- predictors[[chosen_predictors]] %>% 
          crop(hesse_raster)
        final_predictors_df <- as.data.frame(final_predictors, xy = TRUE, na.rm = TRUE)
        
        if(any(c("CORINE", "floodplains", "whc", "heat_index") %in% colnames(final_occ))) {
          final_predictors_df <- final_predictors_df %>%
            mutate(across(any_of(c("CORINE", "floodplains", "whc", "heat_index")), as.factor))
        }
        
        coords <- final_predictors_df[, c("x", "y")]
        
        final_predictors_df$weight <- 1 # only for bake
        final_predictors_df$folds <- 1
        
        final_predictors_df <- final_predictors_df %>%
          dplyr::select(-x,-y)
        
        final_pred_scaled <- bake(rec_prep, new_data = final_predictors_df) %>%
          as.data.frame() %>%
          dplyr::select(-any_of(c("weight","folds")))
        
        common_cols <- intersect(names(final_pred_scaled), names(train_data_scaled))
        
        ordered_common_cols <- names(train_data_scaled)[names(train_data_scaled) %in% common_cols]
        final_pred_scaled <- final_pred_scaled[, ordered_common_cols, drop = FALSE]
        
        final_pred_scaled$x <- coords$x
        final_pred_scaled$y <- coords$y
        
        coordinates(final_pred_scaled) <- ~x + y
        gridded(final_pred_scaled) <- TRUE 
        scaled_raster <- rast(final_pred_scaled) 
        scaled_raster_clean <- scaled_raster[[best_model$finalModel$xNames]] %>% na.omit()
        if (extent == "local") {
          terra::writeRaster(scaled_raster_clean, paste0("Model_results/Model_analyses/Predictor_raster_scaled/scaled_preds_",name_it,".tif"), overwrite = TRUE)
        }
        rm(scaled_raster)
        
        predicted_ras <- scaled_raster_clean %>% terra::predict(best_model, type = "prob", index = 2, na.rm = TRUE)
        predicted_raster_list[[name_it]] <-  predicted_ras
        name_raster <- paste0("predicted_raster",name_it)
        terra::writeRaster(predicted_ras, paste0("Model_results/Maps/",name_raster,".tif"), overwrite = TRUE)
          
        if (extent == "local") {
          #binär
          predicted_bin <- ifel(predicted_ras >= best_threshold, 1, 0)
          names(predicted_bin) <- paste0(name_raster, "_bin")
          predicted_binary_list[[name_it]] <- predicted_bin
          name_raster_bin <- paste0("predicted_raster_binary_",name_it)
          terra::writeRaster(predicted_bin,
                             paste0("Model_results/Maps/binary/", name_raster_bin, ".tif"),
                             overwrite = TRUE)
        }

        if(extent == "local") {
          #Calculate AoA for local models
          #####################################################################################

            aoa_results <- aoa(
              newdata   = scaled_raster_clean,
              model     = best_model,
              variables = "all",
              train     = train_data_scaled,
              method    = "L2",
              useWeight = TRUE,
              useCV     = FALSE,
              LPD       = TRUE,
              verbose   = FALSE,
              algorithm = "kd_tree"
            )

            plot(aoa_results$AOA)

            aoa <- aoa_results$AOA
            aoa_di <- aoparametersaoa_di <- aoa_results$DI
            aoa_pred_weight <- aoa_results$parameters$weight
            aoa_pred_weight <- aoa_pred_weight %>%
              pivot_longer(
                cols = everything(),
                names_to = "predictor",
                values_to = "weight"
              )
            write.csv(aoa_pred_weight, paste0("Model_results/Model_analyses/Non_applicable_analyses/AOA_pred_weight_", name_it, ".csv"), row.names = FALSE)
            aoa_pred_weight_list[[name_it]] <- aoa_pred_weight
            
            aoa_result_list[[name_it]] <- aoa
            writeRaster(aoa, paste0("Model_results/AOA/aoa_", name_it, ".tif"), overwrite = TRUE)
            writeRaster(aoa_di, paste0("Model_results/AOA/DI/aoa_DI_", name_it, ".tif"), overwrite = TRUE)
            
            #Analyses non-applicable space
            ##########################################################################
            predictor_names   <- aoa_pred_weight$predictor        
            pred_rast_scaled  <- scaled_raster_clean[[predictor_names]]  
            train_df_scaled   <- train_data_scaled              
            DI_rast           <- aoa_results$DI                 
            thr               <- aoa_results$parameters$threshold
            top_n             <- 5                               
            
            # build a dataframe of cellwise predictor values + DI using raster values
            pred_vals_mat <- values(pred_rast_scaled)    
            DI_vec <- values(DI_rast)                  
            df_all <- as.data.frame(pred_vals_mat)
            df_all$cell <- seq_len(ncell(DI_rast))
            df_all$DI <- DI_vec
            
            # remove rows where DI is NA or all predictors NA
            df_all <- df_all %>% filter(!is.na(as.vector(DI)))
            
            # keep only non-applicable pixels (DI > threshold)
            vals_nonapp <- df_all %>% filter(DI > thr)
            n_nonapp <- nrow(vals_nonapp)
            
            # 4) Fit LM (DI ~ all scaled predictors)
            lm_mod <- lm(DI ~ ., data = vals_nonapp %>% select(all_of(predictor_names), DI))

            coef_tbl <- broom::tidy(lm_mod) %>%
              filter(term != "(Intercept)") %>%
              rename(predictor = term, estimate = estimate, std_error = std.error, pvalue = p.value) %>%
              mutate(abs_est = abs(estimate)) %>%
              left_join(aoa_pred_weight %>% select(predictor, weight), by = "predictor") %>%
              arrange(desc(abs_est))
            write.csv(coef_tbl, paste0("Model_results/Model_analyses/Non_applicable_analyses/LM_summary_", name_it, ".csv"), row.names = FALSE)
            NA_LM_list[[name_it]] <- coef_tbl
            
            #Select predictors by VarImp
            topN <- coef_tbl %>%
              slice_max(weight, n = top_n) %>%
              mutate(
                predictor_label = paste0(predictor, " (VarImp=", round(weight, 1), ")"),
                predictor_label = factor(predictor_label, levels = rev(predictor_label)),
                direction = ifelse(estimate > 0, "higher → higher DI", "lower → higher DI")
              )
            
            #Plot
            plot_NoNA <- ggplot(topN, aes(x = predictor_label, y = estimate, fill = estimate > 0)) +
              geom_col() +
              coord_flip() +
              scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                                labels = c("higher → higher DI", "lower → higher DI"),
                                name = "Direction") +
              geom_text(aes(
                label = paste0(round(estimate,3), "\n(p=", signif(pvalue,3), ")"),
                y = ifelse(estimate > 0, estimate + 0.02*max(abs(estimate), na.rm=TRUE),
                           estimate - 0.02*max(abs(estimate), na.rm=TRUE))
              ),
              hjust = ifelse(topN$estimate > 0, 0, 1), size = 3  # hier hjust außerhalb aes(), okay so
              ) +
              labs(
                x = "Predictor (VarImp)", 
                y = "Conditional coefficient (scaled predictors)",
                title = paste0("Top ", top_n, " predictors by VarImp: conditional association with DI (non-app pixels)"),
                subtitle = paste0("N non-app pixels = ", n_nonapp)
              ) +
              theme_minimal(base_size = 13)
            ggsave(paste0("Model_results/Model_analyses/Non_applicable_analyses/LM_plot_with_weights_", name_it, ".png"), plot = plot_NoNA, bg = "white", width = 22, height = 12)

        }
        #ERROR messages
        }, error = function(e) {
          message("Error in run ", i, ": ", e$message)

          cat(
            paste0("Run: ", i,
                   " | Species: ", species,
                   " | Extent: ", extent,
                   " | Origin: ", origin,
                   " | Error: ", e$message, "\n"),
            file = paste0("Model_results/Error_messages/error_log_",name_it,".txt"),
            append = TRUE)
        })

      }
    
      #save meaned results to files
      ###############################################################################
      name_species_mean <- paste0("_",extent,"_",origin,"_",species, "_mean")
      
      #ROC
      model_train_roc_df <- data.frame(
        Species = names(model_train_roc_list),
        ROC = unlist(model_train_roc_list))
      write.csv(model_train_roc_df, paste0("Model_results/Metrics/Train_AUC/model_train_roc_values_",name_species_mean,".csv"), row.names = FALSE)
      
      #VarImp
      importance_combined <- bind_rows(feature_importance_list)
      mean_importance <- importance_combined %>%
        group_by(variable) %>%
        summarise(mean_overall = mean(Overall, na.rm = TRUE)) %>%
        arrange(desc(mean_overall))
      write.csv2(mean_importance, file = paste0("Model_results/Feature_importance/fi_",name_species_mean,".csv"), row.names = FALSE)
      
      #AUC
      auc_mean <- mean(unlist(auc_list))
      auc_mean_list[[name_species_mean]] <- auc_mean
      write.csv2(auc_mean, file = paste0("Model_results/Metrics/AUC/auc_",name_species_mean,".csv"), row.names = FALSE)
      tss_mean <- mean(unlist(tss_list))
      tss_mean_list[[name_species_mean]] <- tss_mean
      write.csv2(tss_mean, file = paste0("Model_results/Metrics/TSS/tss_",name_species_mean,".csv"), row.names = FALSE)
      
      #Predicted Raster
      pred_ras_stack <- rast(predicted_raster_list)
      pred_ras_sd <- app(pred_ras_stack, fun = sd, na.rm = TRUE)
      file_path_sd_ras <- paste0("Model_results/Maps/", name_species_mean, "_sd.tif")
      writeRaster(pred_ras_sd, file_path_sd_ras, overwrite = TRUE)
      pred_ras_mean <- app(pred_ras_stack, fun = mean, na.rm = TRUE)
      file_path_mean_ras <- paste0("Model_results/Maps/", name_species_mean, ".tif")
      writeRaster(pred_ras_mean, file_path_mean_ras, overwrite = TRUE)
      predicted_raster_mean[[name_species_mean]] <- pred_ras_mean

      if(extent == "local") {
        #binary maps
        best_threshold_mean <- mean(unlist(best_threshold_list))
        pred_ras_bin_stack <- rast(predicted_binary_list)
        pred_ras_bin_sd <- app(pred_ras_bin_stack, fun = sd, na.rm = TRUE)
        pred_ras_bin_mean <- mean(pred_ras_bin_stack, na.rm = TRUE)
        pred_ras_mean_bin <- ifel(pred_ras_bin_mean >= best_threshold_mean, 1, 0)
        file_path_mean_ras_bin <- paste0("Model_results/Maps/binary/binary_", name_species_mean, ".tif")
        writeRaster(pred_ras_mean_bin, file_path_mean_ras_bin, overwrite = TRUE)
        file_path_sd_ras_bin <- paste0("Model_results/Maps/binary/binary_sd_", name_species_mean, ".tif")
        writeRaster(pred_ras_bin_sd, file_path_sd_ras_bin, overwrite = TRUE)
        
        #AOA
        aoa_stack <- rast(aoa_result_list)
        aoa_sd <- app(aoa_stack, fun = sd, na.rm = TRUE)
        aoa_rast_mean <- mean(aoa_stack, na.rm = TRUE)
        file_aoa_mean <- paste0("Model_results/AOA/aoa_", name_species_mean, ".tif")
        writeRaster(aoa_rast_mean, file_aoa_mean, overwrite = TRUE)
        aoa_rast_bin <- ifel(aoa_rast_mean >= 0.8, 1, 0)
        file_aoa_bin <- paste0("Model_results/AOA/aoa_binary", name_species_mean, ".tif")
        writeRaster(aoa_rast_bin, file_aoa_bin, overwrite = TRUE)
        writeRaster(aoa_sd, paste0("Model_results/AOA/aoa_sd_", name_species_mean, ".tif"), overwrite = TRUE)
        
        #mask applicable area of predicted raster 
        aoa_rast_bin[aoa_rast_bin == 0] <- NA
        pred_mean_masked_aoa <- mask(pred_ras_mean, aoa_rast_bin)
        writeRaster(pred_mean_masked_aoa, paste0("Model_results/Maps/inside_aoa/pred_ras_mean_aoa", name_species_mean, ".tif"), overwrite = TRUE)
        
        #mean aoa predictors weight
        aoa_pred_all_weights <- bind_rows(aoa_pred_weight_list)
        mean_weights_aoa <- aoa_pred_all_weights %>%
          group_by(predictor) %>%                # group by predictor
          summarise(mean_weight = mean(weight, na.rm = TRUE)) %>%
          arrange(desc(mean_weight))
        write.csv(mean_weights_aoa, paste0("Model_results/Model_analyses/Non_applicable_analyses/aoa_predictor_weight",name_species_mean,".csv"), row.names = FALSE)
        
        #LM NoNA mean
        LM_all <- dplyr::bind_rows(NA_LM_list, .id = "iteration")
        
        LM_summ_mean <- LM_all %>%
          group_by(predictor) %>%
          summarise(
            mean_estimate = mean(estimate, na.rm = TRUE),
            sd_estimate   = sd(estimate, na.rm = TRUE),
            mean_pvalue   = mean(pvalue, na.rm = TRUE),
            mean_abs_est  = mean(abs_est, na.rm = TRUE),
            n_iter        = n()   
          ) %>%
          arrange(desc(abs(mean_estimate)))
        
        LM_summ_mean  <- LM_summ_mean %>%
          left_join(mean_weights_aoa, by = "predictor")   
        
        write.csv(LM_summ_mean, paste0("Model_results/Model_analyses/Non_applicable_analyses/LM_summary_",name_species_mean,".csv"), row.names = FALSE)
        
        # top predictors by mean_weight
        topN <- LM_summ_mean %>%
          slice_max(mean_weight, n = top_n) %>%
          mutate(
            predictor_label = paste0(predictor, " (VarImp=", round(mean_weight, 1), ")"),
            predictor_label = factor(predictor_label, levels = rev(predictor_label)),
            direction = ifelse(mean_estimate > 0, "higher → higher DI", "lower → higher DI")
          )
        
        # plot
        p_NonA_mean <- ggplot(topN, aes(x = predictor_label, y = mean_estimate, fill = mean_estimate > 0)) +
          geom_col() +
          coord_flip() +
          scale_fill_manual(
            values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
            labels = c("higher → higher DI", "lower → higher DI"),
            name = "Direction"
          ) +
          geom_text(aes(
            label = paste0(round(mean_estimate, 3), "\n(p≈", signif(mean_pvalue, 3), ")"),
            y = ifelse(mean_estimate > 0,
                       mean_estimate + 0.02*max(abs(mean_estimate), na.rm=TRUE),
                       mean_estimate - 0.02*max(abs(mean_estimate), na.rm=TRUE))
          ),
          hjust = ifelse(topN$mean_estimate > 0, 0, 1), size = 3
          ) +
          labs(
            x = "Predictor (VarImp)",
            y = "Conditional coefficient (scaled predictors)",
            title = paste0("Top ", top_n, " predictors by VarImp (mean over iterations)"),
            subtitle = paste0("N iterations = ", unique(LM_summ_mean$n_iter))
          ) +
          theme_minimal(base_size = 13)
        plot(p_NonA_mean)
        
        ggsave(paste0("Model_results/Model_analyses/Non_applicable_analyses/LM_plot_with_weights_",name_species_mean,".png"), plot = p_NonA_mean, bg = "white", width = 22, height = 12)

      }
    }
  }
  if (extent == "global") {
    if (!file.exists("Workspace_results/workspace_final_models_global.RData")) {
      save.image(file = "Workspace_results/workspace_final_models_global.RData")
    } else {
      message("File already exists")
    }
  } else {
    save.image(file = "Workspace_results/workspace_final_models_local.RData")
  }
}

#save metrics as one table
message(paste("saving mean metrics tables"))
system2(command = "Rscript", args = c("R-Skripte/11_Table_metrics.R"))
