# Regional predictor layers: sources and licences

This folder contains the 14 environmental predictor layers of the regional
(Hesse) models, as used in `10_SDM_Code_Model_final.R`. They were produced
from the original data by scripts `04`–`07` (floodplains and the CORINE crop
in ArcGIS, see README). All layers share one grid: 30 arc-seconds (≈ 1 km),
WGS 84 (EPSG:4326), aligned to CHELSA v2.1 and masked to Hesse.

The standardised (z-scored) predictors are not provided as separate files;
they are computed from these layers within each model run of script 10.

## Licence of the processed layers

The processed layers in this folder are provided under the
[Creative Commons Attribution 4.0 International licence (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
When using them, cite this repository and the original data providers listed
below. The layers are derived products: the original data were rasterised,
resampled, reclassified and/or used to calculate distances and densities
(column "Processing"), and the providers are not responsible for these
changes.

## Sources

| Layer file | Predictor (unit) | Original data | Provider (reference) | Original licence | Processing |
|---|---|---|---|---|---|
| `DGM_hesse.tif` | Elevation (m a.s.l.) | Digital elevation model 1000 m (DGM1000) | © GeoBasis-DE / BKG (2021) | [Datenlizenz Deutschland – Namensnennung 2.0](https://www.govdata.de/dl-de/by-2-0) | Reprojected, bilinearly resampled to the CHELSA grid, masked to Hesse |
| `aspect_hesse.tif` | Aspect (°) | DGM1000 | © GeoBasis-DE / BKG (2021) | Datenlizenz Deutschland – Namensnennung 2.0 | Calculated from DGM1000 (WhiteboxTools), then as elevation |
| `twi_hesse.tif` | Topographic wetness index (dimensionless) | DGM1000 | © GeoBasis-DE / BKG (2021) | Datenlizenz Deutschland – Namensnennung 2.0 | Calculated from DGM1000 (D-infinity flow accumulation and slope, WhiteboxTools), then as elevation |
| `roads_distance_hesse.tif` | Distance to roads (m) | Hesse: ATKIS Basis-DLM, roads (AX_Strasse); outside Hesse: DLM250, roads (ver01_l) | HVBG (2025); © GeoBasis-DE / BKG (2024) | Hesse: [Datenlizenz Deutschland – Zero 2.0](https://www.govdata.de/dl-de/zero-2-0); DLM250: Datenlizenz Deutschland – Namensnennung 2.0 | Layers merged; distance from each cell centre to the nearest road |
| `railroads_distance_hesse.tif` | Distance to railroads (m) | Hesse: ATKIS Basis-DLM, railway lines (AX_Bahnstrecke); outside Hesse: DLM250, railways (ver03_l) | HVBG (2025); © GeoBasis-DE / BKG (2024) | Datenlizenz Deutschland – Zero 2.0; Datenlizenz Deutschland – Namensnennung 2.0 | Layers merged; distance to the nearest railway line |
| `railroads_density_hesse.tif` | Railroad density (m railway line per cell) | ATKIS Basis-DLM, railway lines (AX_Bahnstrecke) | HVBG (2025) | Datenlizenz Deutschland – Zero 2.0 | Length of railway lines per cell |
| `urban_builup_density_hesse.tif` | Urban built-up density (share of cell, 0–1) | ATKIS Basis-DLM, feature group "Siedlung" (settlement areas) | HVBG (2025) | Datenlizenz Deutschland – Zero 2.0 | Polygons merged and dissolved; share of each cell covered |
| `flowing_water_distance_hesse.tif` | Distance to flowing water (m) | Hesse: river network based on the digital landscape model 1:25,000 ([Gewaessernetz_DLM25](https://www.geoportal.hessen.de/mapbender/php/mod_dataISOMetadata.php?outputFormat=iso19139&id=03d99fcb-b11e-5bc8-1981-afe12abe87fb)); outside Hesse: WFD river water bodies ([Wasserkörper-DE, WFD 3rd cycle 2022–2027](https://registry.gdi-de.org/id/de.bund.bfg.csw/268d52a5-8409-4fb4-94ed-98f0c7b71d8d)) | HLNUG (2025a), source note: „Darstellung auf der Grundlage von Daten des Hessischen Landesamtes für Naturschutz, Umwelt und Geologie (HLNUG)"; BfG (2024), source note: „© WasserBLIcK/BfG und Zuständige Behörden der Länder, 2024" | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (both) | Layers merged; distance to the nearest watercourse |
| `standing_water_distance_hesse.tif` | Distance to standing water (m) | Hesse: ATKIS Basis-DLM, standing water bodies (AX_StehendesGewaesser); outside Hesse: WFD lake water bodies ([Wasserkörper-DE, WFD 3rd cycle 2022–2027](https://registry.gdi-de.org/id/de.bund.bfg.csw/268d52a5-8409-4fb4-94ed-98f0c7b71d8d)) | HVBG (2025); BfG (2024), source note: „© WasserBLIcK/BfG und Zuständige Behörden der Länder, 2024" | Datenlizenz Deutschland – Zero 2.0; BfG: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Layers merged; distance to the nearest standing water body |
| `floodplains_hesse.tif` | Flood zone (4 classes) | Flood areas HQ10, HQ100, HQextrem | HLNUG (2025c); source note: „Darstellung auf der Grundlage von Daten des Hessischen Landesamtes für Naturschutz, Umwelt und Geologie (HLNUG)" | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Combined and rasterised to the model grid in ArcGIS |
| `heat_index_hesse.tif` | Heat index: mean summer heat load 2001–2020 (7 classes, "sehr kühl" to "dauerheiß") | Mean heat load MITHBEL20J ([HLNUG climate geodata](https://umweltdaten.hessen.de/klima/geodaten/)), derived from MODIS/Terra land surface temperature (MOD11A2 V006; June–August 2001–2020, 1 km) | HLNUG, Kompetenzzentrum Fernerkundung; source data: Wan, Hook & Hulley (2015), NASA EOSDIS Land Processes DAAC | No access or use restrictions (HLNUG metadata); MODIS data: NASA open data policy, no restrictions on reuse | Rasterised to the model grid (class of each cell) |
| `whc_hesse.tif` | Usable field capacity, 0–100 cm (8 classes: "extrem gering" 0–80 mm to "sehr hoch" > 260 mm, plus settlement/traffic and water) | Usable field capacity in the first metre, soil map 1:50,000 (BFD50; [Nutzbare_Feldkapazitaet100cm_50000](https://www.hlnug.de/registry/spatial/dataset/b76214f0-1d61-9d6a-2c66-0297b69b7ee5)) | HLNUG (2025b); source note: „Darstellung auf der Grundlage von Daten des Hessischen Landesamtes für Naturschutz, Umwelt und Geologie (HLNUG)" | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Settlement and excavation classes merged into one class; rasterised to the model grid |
| `CORINE_hesse.tif` | Land cover (CORINE level-2 classes) | CORINE Land Cover 2018 | © European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA, 2019) | [Copernicus data policy](https://land.copernicus.eu/en/data-policy) (free, full and open access) | Cropped in ArcGIS; reclassified to level-2 classes; nearest-neighbour resampling to the model grid |
| `travel_time_hesse.tif` | Travel time to the nearest city (min) | Global travel time to cities 2015 | Weiss et al. (2018), Malaria Atlas Project | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Cropped and masked to Hesse |

Full references: see Appendix 3.2 of the Supplementary material.

## Not included

- **Predictors of the global model** (CHELSA v2.1 bioclimatic variables
  bio1, bio4, bio15–bio18; Brun et al., 2022; licence CC0 1.0): not
  included because of their size (≈ 3.5 GB). `03_CHELSA.R` recreates them
  from the original CHELSA download.
- **Original (raw) data** in `Predictors_raw/`: available from the
  providers listed above.
