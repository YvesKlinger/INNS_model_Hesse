# INNS_model_Hesse — Invasive plant species distribution models

Species distribution models (SDMs) for three invasive non-native plant
species in the German state of Hesse, comparing predictions built from
citizen-science occurrence data (GBIF) against official state-agency
monitoring data (HLNUG), at both a local (Hesse) and global extent.

Species modeled:

| Species                         | GBIF `speciesKey` |
|----------------------------------|--------------------|
| *Ailanthus altissima*            | 3190653            |
| *Heracleum mantegazzianum*       | 3034825            |
| *Impatiens glandulifera*         | 2891770            |

## Repository structure

```
R-Skripte/                         All analysis scripts, numbered in execution order (see below)
Finale_occ_Tabellen/               Final, filtered GBIF occurrence tables (CSV), global and Hesse
Pseudo_absences/global/pa_area/    Areas for global pseudo-absence selection (input data, made in ArcGIS Pro)
Model result metrics/Metrics.zip   Per-run model performance metrics (AUC, TSS, thresholds, confusion matrices)
requirements.R                     Installs the R packages used across all scripts
INNS_model_Hesse.Rproj             RStudio project file — open this to set the working directory
```

## Setup

1. Open `INNS_model_Hesse.Rproj` in RStudio (this sets the working directory
   to the project root automatically — the scripts assume they run with
   this folder, not `R-Skripte/`, as the working directory). If you're not
   using RStudio, `setwd()` into this folder manually before running
   anything; each script has a comment marking where to do this.
2. Run `source("requirements.R")` to install the required packages.
3. Obtain the input data this pipeline expects (not all of it ships in
   this repository — see "Data availability" below) and recreate the
   directory structure the scripts read from and write to (see "Required
   directory structure" below).

## Pipeline / execution order

| Script | Purpose | Manuscript / Supplementary material |
|---|---|---|
| `01_SDM_Code_GBIF_download.R` | GBIF download queries (global; Hesse) | 2.2 |
| `02_SDM_Code_Occurrence_filter.R` | Cleaning and thinning of the occurrence data | 2.2, App. 2 |
| `03_CHELSA.R` | Scaling of the 19 CHELSA bioclim variables | 2.3, App. 3.2 |
| `04_DGM_und_abhaengige_Variablen.R` | DEM, aspect, TWI, heat index | 2.3, App. 3.2 |
| `05_distance_density_variables.R` | Distance and density variables (roads, railroads, water bodies, urban build-up) | 2.3, App. 3.2 |
| `06_land_use.R` | CORINE land cover | 2.3, App. 3.2 |
| `07_whc_traveltime.R` | Water holding capacity, travel time to cities | 2.3, App. 3.2 |
| `08_correlation_matrix_complete_model_area.R` | Predictor correlation (\|r\| ≥ 0.7), global and regional | 2.3, App. 3.1 |
| `09_PCA_occ_analysis.R` | Sampled environmental space: PCA, Schoener's D, equivalency test | 2.4, Figure 1 |
| `10_SDM_Code_Model_final.R` | Global and regional SDMs (BRT), evaluation, prediction, AOA, DI regression | 2.5, 2.6 |
| `11_Table_metrics.R` | Mean AUC/TSS table (called at the end of script 10) | 3.2 |
| `12_Model_performance_comparison.R` | Beta regressions of AUC/TSS between data sources; AUC/TSS boxplots | 2.5, 3.2, App. 4.1–4.3 |
| `13_VarImp_plots.R` | Predictor importance by category | 3.3, Figure 2 |
| `14_Response_plots.R` | Response plots | App. 5.1–5.4 |
| `15_Plot_NonA.R` | Regression coefficients of DI in non-applicable areas | 3.4, Figure 3, App. 6 |
| `16_Ratio_applicable_suitable_plots.R` | Proportion of each model's AOA classified as suitable (≥ 8 of 10 runs) and extent of the AOA, as grouped bars | 3.4, 3.5, Figure 4 |
| `17_Suitability_maps.R` | Suitability maps (COM) and differences CS − COM, StAg − COM | 3.5, Figure 5 |
| `18_Occurrence_maps.R` | Maps of the number of cleaned, unthinned CS and StAg records per 1-km cell and species | App. 1.2 |

Notes on individual steps:

- **01** — `occ_download()` prepares the download server-side; fetch the
  archive yourself (e.g. via `occ_download_get()`). The downloads used for
  the published results are `0021796-250717081556266`
  (global, https://doi.org/10.15468/dl.ns4htj) and `0021801-250717081556266`
  (Hesse, https://doi.org/10.15468/dl.865r6h), both from 21.07.2025.
  Unpack the downloaded archives to
  `Occurrences/Occurrence_data_GBIF/global/0021796-250717081556266/` and
  `Occurrences/Occurrence_data_GBIF/local/0021801-250717081556266/`
  (scripts 02 and 18 read `occurrence.txt` from there).
- **02** — Writes the filtered tables to `Occurrences/.../filtered/`
  (copies of the GBIF tables in `Finale_occ_Tabellen/`); the printed
  `table(...)` counts are the numbers in App. 2. To run script 10 without
  re-running script 02, copy `Finale_occ_Tabellen/all_species_gbif_global.csv`
  to `Occurrences/Occurrence_data_GBIF/global/filtered/` and
  `all_species_gbif_local.csv` to `Occurrences/Occurrence_data_GBIF/local/filtered/`
  (the HLNUG and combined tables have to be created with script 02).
- **03–07** — Prepare the predictors from the raw data (see "Data
  availability"). The 14 regional predictors are written to
  `Predictors/<name>_hesse.tif`; script 10 reads every `*_hesse.tif` in
  `Predictors/`, so no other `*_hesse.tif` may be placed there. The CHELSA
  layers are written to `Predictors/CHELSA/bio<1-19>.tif`.
- **Floodplains** (`Predictors/floodplains_hesse.tif`) were prepared in
  ArcGIS from the flooding areas HQ10, HQ100 and HQextrem of the flood risk
  management plans (HLNUG, 2025c: Überflutungsflächen HQ10, HQ100, HQextrem
  (HWRMP), https://www.hlnug.de/themen/geografische-informationssysteme/geodienste/wasser)
  as a categorical raster on the 1-km model grid (values 1–4, labelled none,
  low, medium and high flooding in the figures). The **CORINE crop**
  (`Predictors_raw/land_use/CORINE/corine_crop.tif`) was also prepared in
  ArcGIS. No script exists for these steps.
- **Global pseudo-absence area** (`Pseudo_absences/global/pa_area/<species>_ras.tif`,
  read by script 10): made in ArcGIS Pro (13.03.2025) and included here as
  input data. Steps recorded in the raster metadata: XYTableToPoint (global
  filtered occurrences) → Buffer 20 000 m (dissolve all) → Erase of a
  1 000 m buffer → PolygonToRaster (cell centre, snapped to the CHELSA
  bio1 grid). The area is thus a ring 1–20 km around the occurrences.
- **10** — For every extent (global / local) × data origin (gbif / hlnug /
  combined — local only) × species: pseudo-absences, spatially-blocked
  (local) or random (global) cross-validated GBM models with forward
  feature selection, 10 (local) / 5 (global) runs, evaluation (AUC, TSS),
  suitability rasters, AOA (`CAST::aoa`) and regression of the DI on the
  predictors in non-applicable areas. Writes all outputs to `Model_results/`.
- **12** — Reads the per-run metrics from `Model_results/Metrics/` or, if
  that folder does not exist, from `Model result metrics/Metrics.zip`, so
  it can be run without re-running the models. Writes all tables,
  diagnostic plots, model summaries and the App. 4.1 boxplots to
  `Model_results/Model_comparison/` (created by the script).
- **13–17** — Read the outputs of script 10. The difference maps of
  script 17 are computed from suitability maps masked by each model's own
  AOA, so differences are shown only for cells inside the AOAs of both
  compared models (intersection).
- **18** — Reads the raw GBIF download (Hesse) and the HLNUG shapefiles,
  applies the same cleaning as script 02 but without the thinning to one
  record per 1-km cell, and maps the number of records per 1-km cell (model
  grid, `Predictors/DGM_hesse.tif`) per species and data source. Writes to `Model_results/Figures/`; the
  number of occupied cells per panel equals the number of thinned records
  in `Finale_occ_Tabellen/`.

## Data availability and restrictions

- **GBIF occurrence data**: publicly available; see the download DOIs
  above or re-run `01_SDM_Code_GBIF_download.R`. The filtered GBIF tables
  are included in `Finale_occ_Tabellen/`. Each record keeps its original
  licence as published on GBIF (CC0 1.0, CC BY 4.0 or CC BY-NC 4.0; column
  `license`); records under CC BY-NC 4.0 may not be used commercially. When
  using the tables, cite the GBIF downloads
  (https://doi.org/10.15468/dl.ns4htj, https://doi.org/10.15468/dl.865r6h).
- **HLNUG (state agency) occurrence data**: **not included** — it may not
  be redistributed. The data are available from the Hessian Biodiversity
  Database (HEBID) of the Hessian State Agency for Nature Conservation,
  Environment and Geology (HLNUG) on request. This also applies to the
  combined data set (GBIF + HLNUG), which is created by
  `02_SDM_Code_Occurrence_filter.R`. To run the pipeline with the HLNUG
  data, place the three shapefiles as
  `Occurrences/Occurrence_data_HLNUG/unfiltered/<name>/<name>.shp` with
  `<name>` = `Chinesischer_Götterbaum`, `Riesen_Bärenklau`, `Drüsiges_Springkraut`
  (as read by script 02).
- **Environmental predictors**: **not included** — they may not be
  redistributed; the sources are listed below. Scripts 03–07 expect the raw
  data in `Predictors_raw/` and the Hesse outline in `Shape_Hesse/`:

  | Raw data (path below `Predictors_raw/`) | Source (App. 3.2 of the Supplementary material) |
  |---|---|
  | `CHELSA/CHELSA_bio<1-19>_1981-2010_V.2.1.tif` | CHELSA v2.1 (Brun et al., 2022) |
  | `DGM/dgm_original_data/dgm1000/dgm1000_utm32s.asc` | DEM 1000 m (BKG, 2021) |
  | `MITHBEL20J/MITHBEL20J.shp` | Heat index (Wan et al., 2015; HLNUG) |
  | `land_use/ATKIS/ATKIS Basis DLM Hessen (shape)/` (roads, railroads, standing water) | ATKIS Basis-DLM (HVBG, 2025) |
  | `Distance_density_variables/{roads,railroads,urban_buildup,flowing_water,standing_water}/` | BKG (2024), BfG (2024), HLNUG (2025a) |
  | Flooding areas HQ10, HQ100, HQextrem (processed in ArcGIS to `Predictors/floodplains_hesse.tif`) | HLNUG (2025c), https://www.hlnug.de/themen/geografische-informationssysteme/geodienste/wasser |
  | `land_use/CORINE/corine_crop.tif` | CORINE Land Cover 2018 (EEA, 2019), cropped in ArcGIS |
  | `water_holding_capacity/Nutzbare_Feldkapazitaet100cm_50000__epsg25832/` | Usable field capacity (HLNUG, 2025b) |
  | `Travel_time/201501_Global_Travel_Time_to_Cities_2015.tif` | Weiss et al. (2018) |
  | `Shape_Hesse/Hessen.shp` (project root) | Outline of Hesse |

## Required directory structure

The scripts read from and write to a number of subfolders that are **not**
all auto-created and are **not** all present in this repository (most are
regenerated by running the pipeline, or depend on data not included here).
Before running the full pipeline, create:

```
Predictors_raw/                    (raw data, see table above)
Predictors/                        (written by scripts 03-07; floodplains_hesse.tif added by hand)
Shape_Hesse/
Occurrences/Occurrence_data_GBIF/{global,local}/filtered/
Occurrences/Occurrence_data_HLNUG/{unfiltered,filtered}/
Occurrences/Occurrence_data_combined/filtered/
Pseudo_absences/global/pa_area/    (included)
Pseudo_absences/local/
Model_results/best_models/
Model_results/Feature_importance/
Model_results/response_data/
Model_results/Metrics/{AUC,threshold_binary,confusion_matrix,TSS,Train_AUC}/
Model_results/Maps/{binary,inside_aoa}/
Model_results/AOA/DI/
Model_results/Model_analyses/{Predictor_raster_scaled,Non_applicable_analyses}/
Model_results/{pa_area_local,folds_count,Error_messages}/
Workspace_results/
PCA_results/{PCA_loadings,PCA_plots}/
```

## Large files / data availability note

The per-run performance metrics (`Model result metrics/Metrics.zip`,
130 KB) are included; `12_Model_performance_comparison.R` uses them when
`Model_results/Metrics/` does not exist. 

## Not included

No script is available for the following results of the manuscript:

- Table 1 (Schoener's D between the suitability maps).
- The mean absolute DI coefficients given in Results 3.4.
- The selection of candidate species (> 50 records in both data sources, 2.1).
- The preparation of the floodplains predictor and of the CORINE crop (ArcGIS, see above).


## License

The code in this repository is released under the MIT License (see
`LICENSE`). The MIT License applies to the code only. The GBIF occurrence
records in `Finale_occ_Tabellen/` keep their original licences (CC0 1.0,
CC BY 4.0 or CC BY-NC 4.0, column `license`) — see "Data availability and
restrictions" above.

## Citation

See `CITATION.cff`.
