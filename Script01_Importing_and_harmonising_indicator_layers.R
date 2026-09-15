####################################################################################################
# ARCHETYPE ANALYSIS SCRIPT 1 — IMPORTING AND HARMONISING INDICATOR LAYERS
#
# This script prepares the spatial indicator data used in:
# Timberlake et al., "Archetypes of desertification vulnerability in Mediterranean drylands
# and their implications for socio-ecological resilience".
#
# The workflow mirrors the indicator structure used in the manuscript:
#   1. Define the Mediterranean study grid (0.5°, EPSG:4326)
#   2. Prepare candidate exposure indicators
#   3. Prepare candidate sensitivity indicators
#   4. Prepare candidate adaptive-capacity indicators
#   5. Combine all candidate indicators into a single grid-cell table for Script 2
#   6. Produce the indicator maps presented in the Supporting Information
#
# Script 1 deliberately prepares the full candidate indicator set. The final indicator selection
# used for the main clustering is made later, after the correlation checks in Script 2 and the
# explicit selection settings in Script 3. Aridity index is used only to define the dryland analysis
# domain and is not itself a clustering indicator.
#
# All analytical layers are aligned to the same 0.5° longitude-latitude grid. Country polygons
# clipped to the study extent are used for masking, while full country outlines are retained for
# plotting so that the study boundary is not drawn as an artificial national border.
####################################################################################################


# DATA AVAILABILITY AND REPRODUCIBILITY
# --------------------------------------------------------------------------------------------------
# This script documents the processing of the publicly available spatial datasets used to construct
# the Mediterranean dryland indicator dataset. The original source datasets are not included in the
# GitHub repository because several are large geospatial files and are available directly from their
# respective public data providers.
#
# To allow the analytical workflow to be reproduced without downloading and reprocessing all raw
# spatial layers, the repository includes the combined 0.5° grid-cell dataset generated at the end
# of this script:
#
#   output_data/Script1/med_05deg_all_indicators_combined.csv
#
# This file provides the starting point for Script 2 and contains all candidate indicators together
# with the aridity index used to define the dryland analysis domain. Users wishing to reproduce only
# the preprocessing, clustering, sensitivity analyses and manuscript outputs can therefore begin
# directly with Script 2.
#
# Users wishing to reproduce the complete data-construction workflow from the original spatial
# sources should obtain the datasets listed below from the cited public repositories and preserve
# the folder/file structure expected by this script.


rm(list = ls())

#===================================================================================================
# Packages
#===================================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(terra)
  library(sf)
  library(tidyverse)
  library(exactextractr)
  library(rnaturalearth)
  library(geodata)
  library(patchwork)
})

#===================================================================================================
# Folder structure
#===================================================================================================
# I define folders once here and reuse them throughout.
dir_input  <- "input_data"
dir_output <- "output_data/Script1"
dir_plots  <- "plots_script1"

dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots,  showWarnings = FALSE, recursive = TRUE)

#===================================================================================================
# Mediterranean extent + template grid
#===================================================================================================
# I define a single canonical Mediterranean extent for all indicators.
extent_med <- ext(-10, 40, 28, 45)  # xmin, xmax, ymin, ymax

# I also keep an sf bbox version because I use it for clipping/filtering vector layers.
med_bbox <- st_as_sfc(
  st_bbox(
    c(xmin = xmin(extent_med), xmax = xmax(extent_med),
      ymin = ymin(extent_med), ymax = ymax(extent_med)),
    crs = 4326
  )
)

# This is my 0.5° TEMPLATE GRID that everything gets aligned to.
raster_med <- rast(extent_med, resolution = 0.5, crs = "EPSG:4326")

#===================================================================================================
# Country outlines and land mask
#===================================================================================================
# I create two country layers:
# - world_mask is clipped to the study extent and is used only to mask raster values.
# - world_outline is not clipped and is used only to draw national borders.
#
# This prevents the artificial lower edge of the study extent, 28°N, from being drawn as if it were
# a real national border.

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(4326)

# I correct the Cyprus country code so it joins cleanly to World Bank datasets.
world$adm0_a3[world$adm0_a3 == "CYN"] <- "CYP"

# For raster masking only.
world_mask <- st_intersection(world, med_bbox)

# For plotting borders only.
world_outline <- world %>%
  filter(st_intersects(geometry, med_bbox, sparse = FALSE))


#===================================================================================================
# Small plotting/export helper
#===================================================================================================
# I use this helper to save the same plot as both PNG and SVG with consistent settings.
save_png_svg <- function(plot_obj, png_path, svg_path, width = 8, height = 6, dpi = 600) {
  ggsave(plot = plot_obj, filename = png_path, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(plot = plot_obj, filename = svg_path, width = width, height = height, dpi = dpi, bg = "white")
}

#===================================================================================================
# CANDIDATE EXPOSURE INDICATORS
#===================================================================================================

#---------------------------------------------------------------------------------------------------
# Interannual variability in available water supply (Aqueduct 4.0)
#---------------------------------------------------------------------------------------------------
# Here I join Aqueduct interannual variability scores to HydroBASINS level 6 polygons,
# then assign each 0.5° grid cell the Aqueduct score of the HydroBASINS polygon
# containing that grid-cell centre.

aqueduct_data <- fread(
  file.path(dir_input, "Water_stress_index/CVS/Aqueduct40_baseline_annual_y2023m07d05.csv")
)

# The Aqueduct table contains multiple rows per pfaf_id, so I summarise to one value per basin.
aqueduct_iav <- aqueduct_data %>%
  mutate(
    pfaf_id = as.character(pfaf_id),
    iav_score = suppressWarnings(as.numeric(iav_score)),
    iav_score = ifelse(iav_score <= -9999, NA_real_, iav_score)
  ) %>%
  group_by(pfaf_id) %>%
  summarise(iav_score = mean(iav_score, na.rm = TRUE), .groups = "drop")

# Read HydroBASINS level 6 polygons for Africa and Europe/Middle East.
hydro_af <- st_read(
  file.path(dir_input, "Water_stress_index/Hydrobasins_geometry/hybas_af_lev06_v1c.shp"),
  quiet = TRUE
)

hydro_eu <- st_read(
  file.path(dir_input, "Water_stress_index/Hydrobasins_geometry/hybas_eu_lev06_v1c.shp"),
  quiet = TRUE
)

# Combine HydroBASINS regions and make the ID compatible with the Aqueduct pfaf_id column.
hydro6 <- bind_rows(hydro_af, hydro_eu) %>%
  st_transform(4326) %>%
  mutate(pfaf_id = as.character(PFAF_ID)) %>%
  left_join(aqueduct_iav, by = "pfaf_id") %>%
  st_make_valid()

# Make a dataframe of the 0.5° grid-cell centres.
# The warning "[readValues] raster has no values" is harmless because raster_med is only a template.
raster_df_iav <- as.data.frame(raster_med, xy = TRUE, na.rm = FALSE)[, c("x", "y")]

grid_points <- st_as_sf(
  raster_df_iav,
  coords = c("x", "y"),
  crs = 4326,
  remove = FALSE
)

# Assign each grid-cell centre the iav_score of the HydroBASINS polygon it falls inside.
# This avoids clipping/rasterising the HydroBASINS polygons directly.
sf_use_s2(FALSE)

grid_points_iav <- st_join(
  grid_points,
  hydro6 %>% dplyr::select(iav_score),
  join = st_within,
  left = TRUE
)

sf_use_s2(TRUE)

raster_df_iav <- grid_points_iav %>%
  st_drop_geometry() %>%
  dplyr::select(x, y, interannual_variability = iav_score) %>%
  mutate(interannual_variability = as.numeric(interannual_variability))

# Convert back to raster on the same 0.5° template.
iav_raster <- rast(
  raster_df_iav,
  type = "xyz",
  crs = "EPSG:4326"
)

iav_raster <- resample(iav_raster, raster_med, method = "near")
iav_raster <- mask(iav_raster, vect(world_mask))
iav_raster <- clamp(iav_raster, lower = 0, upper = 5, values = TRUE)

raster_df_iav <- as.data.frame(iav_raster, xy = TRUE, na.rm = FALSE)
colnames(raster_df_iav)[3] <- "interannual_variability"

write.csv(
  raster_df_iav,
  file.path(dir_output, "interannual_variability_raster_data.csv"),
  row.names = FALSE
)

iav_plot <- ggplot(raster_df_iav) +
  geom_raster(aes(x = x, y = y, fill = interannual_variability)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D", na.value = "gray90", limits = c(0, 5)) +
  coord_sf(
    xlim = c(xmin(extent_med), xmax(extent_med)),
    ylim = c(ymin(extent_med), ymax(extent_med)),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Interannual variability in water supply",
    fill = "Score"
  )

save_png_svg(
  iav_plot,
  file.path(dir_plots, "Mediterranean_interannual_variability.png"),
  file.path(dir_plots, "Mediterranean_interannual_variability.svg")
)

#---------------------------------------------------------------------------------------------------
# Mean annual precipitation (WorldClim v2.1 BIO12)
#---------------------------------------------------------------------------------------------------
# Here I use WorldClim v2.1 BIO12 because I only need a simple, general representation of
# long-term annual rainfall patterns across the Mediterranean.
#
# BIO12 = annual precipitation, based on the WorldClim 1970–2000 climatology.
#
# Workflow:
# - Download WorldClim bioclimatic variables
# - Extract BIO12, annual precipitation
# - Crop to the Mediterranean extent
# - Resample to the shared 0.5° template grid
# - Mask to land
# - Export raster/CSV/map
#
# Units:
# - BIO12 is annual precipitation in mm/year


dir_worldclim <- file.path(dir_input, "WorldClim")
dir.create(dir_worldclim, showWarnings = FALSE, recursive = TRUE)

# Download WorldClim v2.1 bioclimatic variables.
# res = 10 means 10 arc-minutes, approximately 0.167 degrees.
# This is finer than your 0.5° template but still light enough for easy processing.
worldclim_bio <- geodata::worldclim_global(
  var  = "bio",
  res  = 10,
  path = dir_worldclim
)

# BIO12 = annual precipitation.
mean_annual_ppt <- worldclim_bio[[12]]
names(mean_annual_ppt) <- "mean_annual_precipitation"

# Crop to Mediterranean extent.
mean_annual_ppt_med <- crop(mean_annual_ppt, extent_med)

# Align to the shared 0.5° template grid.
# Since WorldClim is finer than 0.5°, aggregate using mean before resampling.
agg_fact_x <- round(xres(raster_med) / xres(mean_annual_ppt_med))
agg_fact_y <- round(yres(raster_med) / yres(mean_annual_ppt_med))

mean_annual_ppt_agg <- aggregate(
  mean_annual_ppt_med,
  fact = c(agg_fact_x, agg_fact_y),
  fun = mean,
  na.rm = TRUE
)

mean_annual_ppt_raster <- resample(
  mean_annual_ppt_agg,
  raster_med,
  method = "bilinear"
)

# Mask to land only.
mean_annual_ppt_raster <- mask(mean_annual_ppt_raster, vect(world_mask))
names(mean_annual_ppt_raster) <- "mean_annual_precipitation"

# Export raster as GeoTIFF.
writeRaster(
  mean_annual_ppt_raster,
  file.path(dir_output, "mean_annual_precipitation_WorldClim_BIO12_0p5deg.tif"),
  overwrite = TRUE
)

# Export raster values as CSV.
raster_df_ppt <- as.data.frame(mean_annual_ppt_raster, xy = TRUE, na.rm = FALSE)
colnames(raster_df_ppt)[3] <- "mean_annual_precipitation"

write.csv(
  raster_df_ppt,
  file.path(dir_output, "mean_annual_precipitation_WorldClim_BIO12_raster_data.csv"),
  row.names = FALSE
)

# Plot.
ppt_plot <- ggplot(raster_df_ppt) +
  geom_raster(aes(x = x, y = y, fill = mean_annual_precipitation)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(
    option = "D",
    direction = -1, 
    na.value = "gray90",
    name = "mm/year"
  ) +
  coord_sf(
    xlim = c(xmin(extent_med), xmax(extent_med)),
    ylim = c(ymin(extent_med), ymax(extent_med)),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Mean annual precipitation",
    fill = "mm/year"
  )

save_png_svg(
  ppt_plot,
  file.path(dir_plots, "Mediterranean_mean_annual_precipitation_WorldClim_BIO12.png"),
  file.path(dir_plots, "Mediterranean_mean_annual_precipitation_WorldClim_BIO12.svg")
)


#---------------------------------------------------------------------------------------------------
# Baseline water stress (Aqueduct 4.0; candidate indicator retained for sensitivity analysis)
#---------------------------------------------------------------------------------------------------
# Here I join Aqueduct baseline water stress scores to HydroBASINS level 6 polygons,
# then assign each 0.5° grid cell the Aqueduct score of the HydroBASINS polygon
# containing that grid-cell centre.
#
# IMPORTANT:
# Aqueduct uses -9999 for no data. Values of 0 are real values and should be kept.

aqueduct_bws <- aqueduct_data %>%
  mutate(
    pfaf_id = as.character(pfaf_id),
    bws_score = suppressWarnings(as.numeric(bws_score)),
    bws_score = ifelse(is.na(bws_score) | bws_score <= -9999,
                       NA_real_,
                       bws_score)
  ) %>%
  group_by(pfaf_id) %>%
  summarise(
    bws_score = ifelse(
      all(is.na(bws_score)),
      NA_real_,
      mean(bws_score, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Join water stress scores to the same HydroBASINS polygons used above.
hydro6_bws <- hydro6 %>%
  dplyr::select(-any_of("bws_score")) %>%
  left_join(aqueduct_bws, by = "pfaf_id") %>%
  st_make_valid()

# Make a dataframe of the 0.5° grid-cell centres.
# The warning "[readValues] raster has no values" is harmless because raster_med is only a template.
raster_df_bws <- as.data.frame(raster_med, xy = TRUE, na.rm = FALSE)[, c("x", "y")]

grid_points_bws <- st_as_sf(
  raster_df_bws,
  coords = c("x", "y"),
  crs = 4326,
  remove = FALSE
)

# Assign each grid-cell centre the bws_score of the HydroBASINS polygon it falls inside.
sf_use_s2(FALSE)

grid_points_bws <- st_join(
  grid_points_bws,
  hydro6_bws %>% dplyr::select(bws_score),
  join = st_within,
  left = TRUE
)

sf_use_s2(TRUE)

raster_df_bws <- grid_points_bws %>%
  st_drop_geometry() %>%
  dplyr::select(x, y, water_stress = bws_score) %>%
  mutate(
    water_stress = as.numeric(water_stress),
    water_stress = ifelse(is.na(water_stress) | water_stress <= -9999,
                          NA_real_,
                          water_stress)
  )

# Convert back to raster on the same 0.5° template.
bws_raster <- rast(
  raster_df_bws,
  type = "xyz",
  crs = "EPSG:4326"
)

bws_raster <- resample(bws_raster, raster_med, method = "near")
bws_raster <- mask(bws_raster, vect(world_mask))
bws_raster <- clamp(bws_raster, lower = 0, upper = 5, values = TRUE)

raster_df_bws <- as.data.frame(bws_raster, xy = TRUE, na.rm = FALSE)
colnames(raster_df_bws)[3] <- "water_stress"

# Final safety clean: keep 0 as a real value, but remove -9999 no-data values.
raster_df_bws <- raster_df_bws %>%
  mutate(
    water_stress = ifelse(is.na(water_stress) | water_stress <= -9999,
                          NA_real_,
                          water_stress)
  )

write.csv(
  raster_df_bws,
  file.path(dir_output, "water_stress_raster_data.csv"),
  row.names = FALSE
)

water_stress_plot <- ggplot(raster_df_bws) +
  geom_raster(aes(x = x, y = y, fill = water_stress)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D", na.value = "gray90", limits = c(0, 5)) +
  coord_sf(
    xlim = c(xmin(extent_med), xmax(extent_med)),
    ylim = c(ymin(extent_med), ymax(extent_med)),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Water stress index",
    fill = "Score"
  )

save_png_svg(
  water_stress_plot,
  file.path(dir_plots, "Mediterranean_water_stress.png"),
  file.path(dir_plots, "Mediterranean_water_stress.svg")
)


#---------------------------------------------------------------------------------------------------
# Population density (log transformed)
#---------------------------------------------------------------------------------------------------
# Here I calculate the area-weighted mean population density within each 0.5° grid cell,
# then log-transform the resulting population density for clustering.
#
# The original GPW raster is at 30 arc-second resolution and reports population density
# in people/km2. Using an area-weighted mean ensures that all fine-resolution pixels within
# each 0.5° grid cell contribute to the final value, rather than interpolating values using
# bilinear resampling.

population <- rast(file.path(
  dir_input,
  "gpw_v4_population_density_rev11_2020_30_sec.tif"
))

# Crop first, then only project if needed.
population_clipped <- crop(population, extent_med)

if (!grepl("4326", crs(population_clipped))) {
  population_clipped <- project(
    population_clipped,
    "EPSG:4326",
    method = "bilinear"
  )
  population_clipped <- crop(population_clipped, extent_med)
}

### Create 0.5° grid polygons

# Reuse the existing template-grid polygons if they have already been created elsewhere in the script.
if (!exists("grid_polys") || !exists("grid_sf")) {
  grid_polys <- as.polygons(raster_med, dissolve = FALSE)
  grid_sf <- st_as_sf(grid_polys)
}

### Calculate source-pixel area

# Longitude/latitude pixels vary slightly in area with latitude, so calculate the area
# of each source pixel in km2 and use this as a weighting raster.
population_pixel_area_km2 <- cellSize(
  population_clipped,
  unit = "km"
)

### Aggregate population density to the 0.5° grid

# Calculate the area-weighted mean population density within each 0.5° grid cell.
population_density_05 <- exactextractr::exact_extract(
  population_clipped,
  grid_sf,
  "weighted_mean",
  weights = population_pixel_area_km2,
  progress = TRUE
)

# Put the aggregated values back onto the shared 0.5° template grid.
population_density_raster <- raster_med
values(population_density_raster) <- population_density_05
names(population_density_raster) <- "population_density_km2"

# Mask to land within the Mediterranean study extent.
population_density_raster <- mask(
  population_density_raster,
  vect(world_mask)
)

# Remove any negative artefacts if present.
population_density_raster <- clamp(
  population_density_raster,
  lower = 0,
  values = TRUE
)

### Log-transform for clustering

population_resampled_log <- log(population_density_raster + 1)
names(population_resampled_log) <- "pop_density_log"


### Export population-density data

# Retain both the original population density and the log-transformed version.
raster_df_population <- as.data.frame(
  c(population_resampled_log, population_density_raster),
  xy = TRUE,
  na.rm = FALSE
)

write.csv(
  raster_df_population,
  file.path(dir_output, "population_raster_data.csv"),
  row.names = FALSE
)

### Plot

pop_density_med <- ggplot(raster_df_population) +
  geom_raster(aes(x = x, y = y, fill = pop_density_log)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D", na.value = "gray90") +
  coord_sf(xlim = c(xmin(extent_med), xmax(extent_med)),
           ylim = c(ymin(extent_med), ymax(extent_med)),
           expand = FALSE) +
  theme_minimal() +
  labs(title = "Population density (log)", fill = "Log Pop. Density")

save_png_svg(
  pop_density_med,
  file.path(dir_plots, "Mediterranean_population_density.png"),
  file.path(dir_plots, "Mediterranean_population.svg")
)

#---------------------------------------------------------------------------------------------------
# Aridity index (0.5° for analysis; 0.1° for plotting)
#---------------------------------------------------------------------------------------------------
# Here I prepare the aridity index for analysis on the 0.5° grid,
# and I also generate the 0.1° maps for visual interpretation.
# This layer is exported for downstream clipping/filtering, but is not used as an archetype indicator.

aridity_index <- rast(file.path(dir_input, "Zomer_et_al_Aridty_Index_data.tif"))

# Crop first, then only project if needed.
aridity_index_clipped <- crop(aridity_index, extent_med)

if (!grepl("4326", crs(aridity_index_clipped))) {
  aridity_index_clipped <- project(aridity_index_clipped, "EPSG:4326", method = "bilinear")
  aridity_index_clipped <- crop(aridity_index_clipped, extent_med)
}

# I apply the scale factor and remove non-terrestrial zeros.
aridity_index_clipped <- aridity_index_clipped * 0.0001
aridity_index_clipped[aridity_index_clipped == 0] <- NA

# 0.5° resample for analysis
aridity_index_resampled <- resample(aridity_index_clipped, raster_med, method = "bilinear")

# Mask to land within the Mediterranean study extent.
aridity_index_resampled <- mask(aridity_index_resampled, vect(world_mask))

raster_df_aridity_index <- as.data.frame(aridity_index_resampled, xy = TRUE, na.rm = FALSE)
colnames(raster_df_aridity_index)[3] <- "aridity_index"
write.csv(raster_df_aridity_index, file.path(dir_output, "Aridity_Index_raster_data.csv"), row.names = FALSE)

# 0.1° resample only for plotting detail
raster_med_0p1 <- rast(ext(raster_med), crs = crs(raster_med), resolution = 0.1)
aridity_index_resampled_0p1 <- resample(aridity_index_clipped, raster_med_0p1, method = "bilinear")

# Mask plotting raster to land as well.
aridity_index_resampled_0p1 <- mask(aridity_index_resampled_0p1, vect(world_mask))

raster_df_aridity_index_0p1 <- as.data.frame(aridity_index_resampled_0p1, xy = TRUE, na.rm = FALSE)
colnames(raster_df_aridity_index_0p1)[3] <- "aridity_index"

# Categorical version for interpretability (drylands bands)
raster_df_aridity_index_0p1_cats <- raster_df_aridity_index_0p1 %>%
  mutate(category = case_when(
    is.na(aridity_index) ~ NA_character_,
    aridity_index < 0.03 ~ "Hyper-Arid",
    aridity_index >= 0.03 & aridity_index < 0.2 ~ "Arid",
    aridity_index >= 0.2 & aridity_index < 0.5 ~ "Semi-Arid",
    aridity_index >= 0.5 & aridity_index < 0.65 ~ "Dry Sub-Humid"
  )) %>%
  mutate(category = factor(category,
                           levels = c("Hyper-Arid", "Arid", "Semi-Arid", "Dry Sub-Humid")))

color_palette <- c("Hyper-Arid" = "#FF7300",
                   "Arid" = "#FFD700",
                   "Semi-Arid" = "#FFFACD",
                   "Dry Sub-Humid" = "#ADD8E6")

aridity_index_med <- ggplot(raster_df_aridity_index_0p1_cats) +
  geom_raster(aes(x = x, y = y, fill = category)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_manual(values = color_palette, na.value = "white") +
  coord_sf(xlim = c(xmin(extent_med), xmax(extent_med)),
           ylim = c(ymin(extent_med), ymax(extent_med)),
           expand = FALSE) +
  theme_minimal() +
  labs(title = "Aridity Index - 0.1° resolution", fill = "Aridity class")

save_png_svg(
  aridity_index_med,
  file.path(dir_plots, "Mediterranean_aridity_index.png"),
  file.path(dir_plots, "Mediterranean_aridity_index.svg")
)

# Continuous version
aridity_index_med_cont <- ggplot(raster_df_aridity_index_0p1) +
  geom_raster(aes(x = x, y = y, fill = aridity_index)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D", direction = -1, na.value = "white", name = "Aridity Index") +
  coord_sf(xlim = c(xmin(extent_med), xmax(extent_med)),
           ylim = c(ymin(extent_med), ymax(extent_med)),
           expand = FALSE) +
  theme_minimal() +
  labs(title = "Aridity index")

save_png_svg(
  aridity_index_med_cont,
  file.path(dir_plots, "Mediterranean_aridity_index_continuous.png"),
  file.path(dir_plots, "Mediterranean_aridity_index_continuous.svg")
)

#===================================================================================================
# CANDIDATE SENSITIVITY INDICATORS
#===================================================================================================

#---------------------------------------------------------------------------------------------------
# Land-cover processing: cropland percentage + irrigated percentage + land-cover diversity
# Land-cover diversity is prepared here as a candidate indicator and retained for sensitivity analysis.
#---------------------------------------------------------------------------------------------------
# Here I crop the Copernicus land-cover raster once, group detailed land-cover classes,
# then calculate cropland percentage, irrigated percentage and land-cover diversity per 0.5° cell
# in one extraction.

landcover <- rast(file.path(dir_input, "Copernicus_Landcover_data.nc"))
lc <- landcover[["lccs_class"]]

# Crop first, then only project if needed.
lc_med <- crop(lc, extent_med)

if (!grepl("4326", crs(lc_med))) {
  lc_med <- project(lc_med, "EPSG:4326", method = "near")
  lc_med <- crop(lc_med, extent_med)
}

# If processed_flag exists, use it to mask out unprocessed pixels.
if ("processed_flag" %in% names(landcover)) {
  proc <- landcover[["processed_flag"]]
  
  proc_m <- crop(proc, extent_med)
  
  if (!grepl("4326", crs(proc_m))) {
    proc_m <- project(proc_m, "EPSG:4326", method = "near")
    proc_m <- crop(proc_m, extent_med)
  }
  
  lc_med <- mask(lc_med, proc_m, maskvalues = 0)
}

# Check which Copernicus land-cover codes are present after cropping.
fv <- freq(lc_med)
if (is.matrix(fv)) fv <- as.data.frame(fv)
fv <- fv %>% arrange(value) %>% dplyr::select(value, count)

lut_c3s <- tribble(
  ~value, ~class,
  10, "Cropland, rainfed",
  11, "Cropland, herbaceous cover",
  12, "Cropland, tree or shrub cover",
  20, "Cropland, irrigated or post-flooding",
  30, "Mosaic cropland (>50%) / natural vegetation (<50%)",
  40, "Mosaic natural vegetation (>50%) / cropland (<50%)",
  50, "Tree cover, broadleaved, evergreen, closed to open (>15%)",
  60, "Tree cover, broadleaved, deciduous, closed to open (>15%)",
  61, "Tree cover, broadleaved, deciduous, closed (>40%)",
  62, "Tree cover, broadleaved, deciduous, open (15–40%)",
  70, "Tree cover, needleleaved, evergreen, closed to open (>15%)",
  71, "Tree cover, needleleaved, evergreen, closed (>40%)",
  80, "Tree cover, needleleaved, evergreen, open (15–40%)",
  90, "Tree cover, needleleaved, deciduous",
  100, "Tree cover, mixed leaf type",
  110, "Mosaic tree & shrub (>50%) / herbaceous (<50%)",
  120, "Mosaic herbaceous (>50%) / tree & shrub (<50%)",
  121, "Shrubland, evergreen",
  122, "Shrubland, deciduous",
  130, "Grassland",
  140, "Lichens and mosses",
  150, "Sparse vegetation (tree/shrub/herbaceous)",
  153, "Sparse herbaceous vegetation",
  160, "Tree cover, flooded, fresh or brackish water",
  170, "Tree cover, flooded, saline water (mangroves)",
  180, "Bare areas",
  190, "Consolidated bare areas",
  200, "Unconsolidated bare areas",
  201, "Sandy areas",
  202, "Stony areas",
  210, "Water bodies",
  220, "Permanent snow and ice"
)

codes_table <- fv %>%
  left_join(lut_c3s, by = "value") %>%
  mutate(
    class = ifelse(is.na(class), paste0("UNMAPPED_", value), class),
    pct = 100 * count / sum(count, na.rm = TRUE)
  ) %>%
  arrange(desc(count)) %>%
  as_tibble()

print(codes_table, n = nrow(codes_table))
setdiff(fv$value, lut_c3s$value)

# Reclassify detailed Copernicus classes into broader functional groups.
# Cropland is split into subgroups because mixed agricultural landscapes
# may be more resilient than single-crop landscapes.
#
# Groups:
# 1 = Rainfed/herbaceous cropland
# 2 = Tree/shrub cropland
# 3 = Irrigated cropland
# 4 = Mosaic agriculture/natural vegetation
# 5 = Forest/tree cover
# 6 = Shrubland
# 7 = Grassland/herbaceous
# 8 = Sparse/bare
# Water/snow/ice classes are set to NA and excluded from terrestrial diversity.

landcover_group_lut <- matrix(c(
  10,  1,
  11,  1,
  12,  2,
  20,  3,
  30,  4,
  40,  4,
  50,  5,
  60,  5,
  61,  5,
  62,  5,
  70,  5,
  71,  5,
  80,  5,
  90,  5,
  100, 5,
  160, 5,
  170, 5,
  110, 6,
  120, 6,
  121, 6,
  122, 6,
  130, 7,
  140, 7,
  150, 8,
  153, 8,
  180, 8,
  190, 8,
  200, 8,
  201, 8,
  202, 8,
  210, NA,
  220, NA
), ncol = 2, byrow = TRUE)

lc_grouped <- classify(lc_med, landcover_group_lut, others = NA)

# I turn the template raster into polygons for exact_extract, reusing if already created.
if (!exists("grid_polys") || !exists("grid_sf")) {
  grid_polys <- as.polygons(raster_med, dissolve = FALSE)
  grid_sf <- st_as_sf(grid_polys)
}

n_landcover_groups <- 8

# Calculate all land-cover metrics in one exact_extract pass:
# - prop_cropland: proportion of valid terrestrial land-cover area in groups 1, 2, or 3
# - prop_irrigated: proportion of valid terrestrial land-cover area in group 3
# - landcover_diversity: normalised Shannon diversity across groups 1–8
lc_metrics_list <- exact_extract(
  lc_grouped,
  grid_sf,
  function(values, coverage_fractions) {
    
    ok <- !is.na(values) & !is.na(coverage_fractions) & coverage_fractions > 0
    
    if (!any(ok)) {
      return(data.frame(
        prop_cropland = NA_real_,
        prop_irrigated = NA_real_,
        landcover_diversity = NA_real_
      ))
    }
    
    values <- values[ok]
    coverage_fractions <- coverage_fractions[ok]
    
    total_area <- sum(coverage_fractions, na.rm = TRUE)
    
    if (is.na(total_area) || total_area <= 0) {
      return(data.frame(
        prop_cropland = NA_real_,
        prop_irrigated = NA_real_,
        landcover_diversity = NA_real_
      ))
    }
    
    prop_cropland <- sum(
      coverage_fractions[values %in% c(1, 2, 3)],
      na.rm = TRUE
    ) / total_area
    
    prop_irrigated <- sum(
      coverage_fractions[values == 3],
      na.rm = TRUE
    ) / total_area
    
    group_area <- tapply(
      coverage_fractions,
      factor(values, levels = 1:n_landcover_groups),
      sum,
      na.rm = TRUE
    )
    
    if (all(is.na(group_area)) || sum(group_area, na.rm = TRUE) <= 0) {
      return(data.frame(
        prop_cropland = prop_cropland,
        prop_irrigated = prop_irrigated,
        landcover_diversity = NA_real_
      ))
    }
    
    p <- group_area[group_area > 0] / sum(group_area, na.rm = TRUE)
    
    shannon <- -sum(p * log(p), na.rm = TRUE)
    shannon_norm <- shannon / log(n_landcover_groups)
    
    data.frame(
      prop_cropland = prop_cropland,
      prop_irrigated = prop_irrigated,
      landcover_diversity = shannon_norm
    )
  },
  progress = TRUE
)

lc_metrics <- dplyr::bind_rows(lc_metrics_list)

grid_sf$prop_cropland <- lc_metrics$prop_cropland
grid_sf$prop_irrigated <- lc_metrics$prop_irrigated
grid_sf$landcover_diversity <- lc_metrics$landcover_diversity

#---------------------------------------------------------------------------------------------------
# Export and plot % cropland
#---------------------------------------------------------------------------------------------------

prop_crop_raster <- rasterize(vect(grid_sf), raster_med, field = "prop_cropland")
prop_crop_percent <- prop_crop_raster * 100

prop_crop_percent <- mask(prop_crop_percent, vect(world_mask))
prop_crop_percent <- clamp(prop_crop_percent, lower = 0, upper = 100, values = TRUE)

raster_df_cropland <- as.data.frame(prop_crop_percent, xy = TRUE, na.rm = FALSE)
colnames(raster_df_cropland)[3] <- "cropland_percent"

write.csv(
  raster_df_cropland,
  file.path(dir_output, "crop_perc_raster_data.csv"),
  row.names = FALSE
)

crop_med_plot <- ggplot(raster_df_cropland) +
  geom_raster(aes(x = x, y = y, fill = cropland_percent)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(
    option = "D",
    na.value = "gray90",
    name = "% Cropland",
    limits = c(0, 100)
  ) +
  coord_sf(
    xlim = c(xmin(extent_med), xmax(extent_med)),
    ylim = c(ymin(extent_med), ymax(extent_med)),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(title = "Percentage cropland")

save_png_svg(
  crop_med_plot,
  file.path(dir_plots, "Mediterranean_crop_perc_cover.png"),
  file.path(dir_plots, "Mediterranean_crop_perc_cover.svg")
)

#---------------------------------------------------------------------------------------------------
# Export and plot % irrigated land
#---------------------------------------------------------------------------------------------------

prop_irrigated_raster <- rasterize(vect(grid_sf), raster_med, field = "prop_irrigated")
prop_irrigated_percent <- prop_irrigated_raster * 100

prop_irrigated_percent <- mask(prop_irrigated_percent, vect(world_mask))
prop_irrigated_percent <- clamp(prop_irrigated_percent, lower = 0, upper = 100, values = TRUE)

raster_df_irrigated <- as.data.frame(prop_irrigated_percent, xy = TRUE, na.rm = FALSE)
colnames(raster_df_irrigated)[3] <- "irrigated_percent"

write.csv(
  raster_df_irrigated,
  file.path(dir_output, "irrigated_perc_raster_data.csv"),
  row.names = FALSE
)

raster_df_irrigated <- raster_df_irrigated %>%
  mutate(
    irrigated_percent_log = ifelse(
      is.na(irrigated_percent),
      NA_real_,
      log10(irrigated_percent + 1)
    )
  )

irrigated_plot <- ggplot(raster_df_irrigated) +
  geom_raster(aes(x = x, y = y, fill = irrigated_percent_log)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(
    option = "D",
    na.value = "gray90",
    name = "log10(% + 1)",
    limits = c(0, log10(101)),
    breaks = log10(c(1, 2, 5, 11, 51, 101)),
    labels = c("0", "1", "4", "10", "50", "100")
  ) +
  coord_sf(
    xlim = c(xmin(extent_med), xmax(extent_med)),
    ylim = c(ymin(extent_med), ymax(extent_med)),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(title = "Irrigated agriculture")

save_png_svg(
  irrigated_plot,
  file.path(dir_plots, "Mediterranean_irrigated_perc_cover.png"),
  file.path(dir_plots, "Mediterranean_irrigated_perc_cover.svg")
)

#---------------------------------------------------------------------------------------------------
# Export and plot land-cover diversity
#---------------------------------------------------------------------------------------------------

landcover_diversity_raster <- rasterize(
  vect(grid_sf),
  raster_med,
  field = "landcover_diversity"
)

landcover_diversity_raster <- mask(landcover_diversity_raster, vect(world_mask))
landcover_diversity_raster <- clamp(
  landcover_diversity_raster,
  lower = 0,
  upper = 1,
  values = TRUE
)

raster_df_landcover_diversity <- as.data.frame(
  landcover_diversity_raster,
  xy = TRUE,
  na.rm = FALSE
)

colnames(raster_df_landcover_diversity)[3] <- "landcover_diversity"

write.csv(
  raster_df_landcover_diversity,
  file.path(dir_output, "landcover_diversity_raster_data.csv"),
  row.names = FALSE
)

landcover_diversity_plot <- ggplot(raster_df_landcover_diversity) +
  geom_raster(aes(x = x, y = y, fill = landcover_diversity)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(
    option = "D",
    na.value = "gray90",
    limits = c(0, 1)
  ) +
  coord_sf(
    xlim = c(xmin(extent_med), xmax(extent_med)),
    ylim = c(ymin(extent_med), ymax(extent_med)),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Land-cover diversity",
    fill = "Diversity"
  )

save_png_svg(
  landcover_diversity_plot,
  file.path(dir_plots, "Mediterranean_landcover_diversity.png"),
  file.path(dir_plots, "Mediterranean_landcover_diversity.svg")
)

#---------------------------------------------------------------------------------------------------
# Grazing intensity / livestock density (cattle + sheep + goats; log transform)
#---------------------------------------------------------------------------------------------------
# The GLW 3 rasters store absolute numbers of animals per original 0.083333-degree pixel,
# not densities. Therefore, I first divide each species raster by its pixel area in km2 to obtain
# animals/km2. I then resample those density rasters to the shared 0.5-degree template grid,
# sum cattle + sheep + goats, and log-transform the combined livestock density.
#
# Output variables:
# - livestock_density_km2 = cattle + sheep + goats per km2
# - grazing_intensity_log = log(livestock_density_km2 + 1)

cattle_file <- file.path(dir_input, "Livestock_density_Gilbert_2018/cattle_density_2010_Da.tif")
sheep_file  <- file.path(dir_input, "Livestock_density_Gilbert_2018/sheep_density_2010_Da.tif")
goats_file  <- file.path(dir_input, "Livestock_density_Gilbert_2018/goat_density_2010_Da.tif")

graz_tif <- file.path(dir_output, "grazing_intensity_log_density_0p5deg.tif")
graz_csv <- file.path(dir_output, "grazing_intensity_raster_data_0p5deg.csv")
graz_png <- file.path(dir_plots,  "Mediterranean_grazing_intensity_0p5deg.png")
graz_svg <- file.path(dir_plots,  "Mediterranean_grazing_intensity_0p5deg.svg")

cattle <- rast(cattle_file)
sheep  <- rast(sheep_file)
goats  <- rast(goats_file)

# Crop first, then only project if needed.
cattle_med <- crop(cattle, extent_med)
if (!grepl("4326", crs(cattle_med))) {
  cattle_med <- project(cattle_med, "EPSG:4326", method = "bilinear")
  cattle_med <- crop(cattle_med, extent_med)
}

sheep_med <- crop(sheep, extent_med)
if (!grepl("4326", crs(sheep_med))) {
  sheep_med <- project(sheep_med, "EPSG:4326", method = "bilinear")
  sheep_med <- crop(sheep_med, extent_med)
}

goats_med <- crop(goats, extent_med)
if (!grepl("4326", crs(goats_med))) {
  goats_med <- project(goats_med, "EPSG:4326", method = "bilinear")
  goats_med <- crop(goats_med, extent_med)
}

# Calculate source-pixel area in km2.
# This is needed because GLW values are animals per pixel, not animals per km2.
# terra::cellSize accounts for the changing area of lon/lat pixels with latitude.
cattle_area_km2 <- cellSize(cattle_med, unit = "km")
sheep_area_km2  <- cellSize(sheep_med,  unit = "km")
goats_area_km2  <- cellSize(goats_med,  unit = "km")

# Convert each species from animals/source-pixel to animals/km2.
cattle_density_med <- cattle_med / cattle_area_km2
sheep_density_med  <- sheep_med  / sheep_area_km2
goats_density_med  <- goats_med  / goats_area_km2

names(cattle_density_med) <- "cattle_density_km2"
names(sheep_density_med)  <- "sheep_density_km2"
names(goats_density_med)  <- "goats_density_km2"

# Resample densities, not raw animal counts, to the shared 0.5-degree template.
cattle_density_0p5 <- resample(cattle_density_med, raster_med, method = "bilinear")
sheep_density_0p5  <- resample(sheep_density_med,  raster_med, method = "bilinear")
goats_density_0p5  <- resample(goats_density_med,  raster_med, method = "bilinear")

# Sum cattle + sheep + goats as an approximate combined livestock density.
livestock_density_km2 <- cattle_density_0p5 + sheep_density_0p5 + goats_density_0p5
names(livestock_density_km2) <- "livestock_density_km2"

# Mask to land within the Mediterranean study extent.
livestock_density_km2 <- mask(livestock_density_km2, vect(world_mask))

# Remove any negative artefacts if present.
livestock_density_km2 <- clamp(livestock_density_km2, lower = 0, values = TRUE)

# Log-transform for clustering.
grazing_log <- log(livestock_density_km2 + 1)
names(grazing_log) <- "grazing_intensity_log"

# Export log-transformed raster.
writeRaster(grazing_log, graz_tif, overwrite = TRUE)

# Export both the log-transformed variable and the interpretable original-density variable.
raster_df_graz <- as.data.frame(
  c(grazing_log, livestock_density_km2),
  xy = TRUE,
  na.rm = FALSE
)

write.csv(raster_df_graz, graz_csv, row.names = FALSE)

# Plot log-transformed livestock density.
graz_plot <- ggplot(raster_df_graz) +
  geom_raster(aes(x = x, y = y, fill = grazing_intensity_log)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D", na.value = "gray90") +
  coord_sf(
    xlim = c(xmin(extent_med), xmax(extent_med)),
    ylim = c(ymin(extent_med), ymax(extent_med)),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Livestock density",
    fill  = "log(animals/km² + 1)"
  )

save_png_svg(graz_plot, graz_png, graz_svg)

#===================================================================================================
# Soil Organic Carbon Stock, 0-30 cm (SoilGrids 250 m)
#===================================================================================================

# SoilGrids remote raster

# Use GDAL's /vsicurl/ virtual file system.
# This is often more reliable than passing the https URL directly to terra.
soc_url <- paste0(
  "/vsicurl?list_dir=no&max_retry=3&retry_delay=2&url=",
  "https://files.isric.org/soilgrids/latest/data/ocs/ocs_0-30cm_mean.vrt"
)

# Open remote SoilGrids VRT
soc_raw <- rast(soc_url)

# Optional checks
soc_raw
crs(soc_raw)
ext(soc_raw)

# Crop SoilGrids raster to the Mediterranean extent

# Create Mediterranean extent polygon in lon/lat
med_poly_ll <- as.polygons(extent_med, crs = "EPSG:4326")

# Project the Mediterranean polygon to the native SoilGrids CRS
med_poly_sg <- project(med_poly_ll, crs(soc_raw))

# Crop using the projected extent rather than the polygon object.
# This avoids some failures when reading remote VRT tiles.
soc_crop_raw <- crop(
  soc_raw,
  ext(med_poly_sg),
  snap = "out"
)

# Convert units and reproject to EPSG:4326

# SoilGrids stores ocs as integer values that need to be divided by 10.
# After division, values are kg C / m2.
soc_crop_kg_m2 <- soc_crop_raw / 10
names(soc_crop_kg_m2) <- "soc_0_30cm_kg_m2"

# Reproject cropped raster to lon/lat.
# Use bilinear interpolation because SOC stock is continuous.
soc_crop_ll <- project(
  soc_crop_kg_m2,
  "EPSG:4326",
  method = "bilinear"
)

# Crop exactly to the Mediterranean lon/lat extent
soc_crop_ll <- crop(soc_crop_ll, extent_med)

# Save cropped intermediate raster locally.
# This prevents repeated remote reads if later steps fail.
writeRaster(
  soc_crop_ll,
  file.path(dir_output, "soilgrids_soc_0_30cm_cropped_250m_kg_m2.tif"),
  overwrite = TRUE
)

# Aggregate to the 0.5° template grid

# Create 0.5° grid-cell polygons from the template raster
grid_05_poly <- as.polygons(raster_med, dissolve = FALSE)
grid_05_poly$cell_id <- seq_len(ncell(raster_med))

# Convert to sf for exactextractr
grid_05_sf <- st_as_sf(grid_05_poly)

# exactextractr works most reliably here with a RasterLayer object.
soc_crop_raster <- raster::raster(
  file.path(dir_output, "soilgrids_soc_0_30cm_cropped_250m_kg_m2.tif")
)

# Calculate mean SOC stock within each 0.5 degree grid cell
soc_mean_by_cell <- exactextractr::exact_extract(
  soc_crop_raster,
  grid_05_sf,
  "mean"
)

# Put values back onto the existing 0.5° template grid
soc_raster <- raster_med
values(soc_raster) <- soc_mean_by_cell
names(soc_raster) <- "soc_0_30cm_kg_m2"

# Mask to land within the Mediterranean region
soc_raster <- mask(soc_raster, vect(world_mask))

# Export raster data

# Convert raster to dataframe for export and plotting
raster_df_soc <- as.data.frame(soc_raster, xy = TRUE, na.rm = FALSE)
colnames(raster_df_soc)[3] <- "soc_0_30cm_kg_m2"

# Optional equivalent unit: tonnes C / ha
# 1 kg C / m2 = 10 tonnes C / ha
raster_df_soc <- raster_df_soc %>%
  mutate(soc_0_30cm_t_ha = soc_0_30cm_kg_m2 * 10)

write.csv(
  raster_df_soc,
  file.path(dir_output, "soil_organic_carbon_stock_0_30cm_raster_data.csv"),
  row.names = FALSE
)

# Plot

soc_plot <- ggplot(raster_df_soc) +
  geom_raster(aes(x = x, y = y, fill = soc_0_30cm_kg_m2)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D", na.value = "gray90") +
  coord_sf(
    xlim = c(xmin(extent_med), xmax(extent_med)),
    ylim = c(ymin(extent_med), ymax(extent_med)),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Soil organic carbon stock",
    fill = expression("SOC stock"~"(kg C m"^{-2}~")")
  )

save_png_svg(
  soc_plot,
  file.path(dir_plots, "Mediterranean_SOC_stock_0_30cm.png"),
  file.path(dir_plots, "Mediterranean_SOC_stock_0_30cm.svg")
)

#===================================================================================================
# CANDIDATE ADAPTIVE-CAPACITY INDICATORS
#===================================================================================================

#---------------------------------------------------------------------------------------------------
# Governance (World Bank) — vector plot + rasterised to 0.5°
#---------------------------------------------------------------------------------------------------
# Here I calculate an average governance percentile rank per country (most recent across indicators),
# then rasterise country values onto the 0.5° template grid.

governance_data <- data.table(
  read_excel(file.path(dir_input, "World_Bank_Governance_Indicator_data.xlsx"), sheet = "Sheet1")
)

governance_data_sorted <- governance_data %>%
  arrange(Country_Name, desc(Year))

governance_data_recent <- governance_data_sorted %>%
  group_by(Country_Name, Country_Code, indicator) %>%
  slice(1) %>%
  ungroup()

governance_data_recent$pctrank <- as.numeric(governance_data_recent$pctrank)

governance_data_summary <- governance_data_recent %>%
  group_by(Country_Name, Country_Code) %>%
  summarise(governance_pc_rank = mean(pctrank, na.rm = TRUE), .groups = "drop")

governance_data_summary$Country_Code <- toupper(governance_data_summary$Country_Code)

world_governance_data <- world %>%
  left_join(governance_data_summary, by = c("adm0_a3" = "Country_Code")) %>%
  st_transform(crs = 4326)

med_governance_data <- world_governance_data %>%
  filter(st_intersects(geometry, med_bbox, sparse = FALSE))

governance_med <- ggplot(med_governance_data) +
  geom_sf(aes(fill = governance_pc_rank), color = "black", size = 0.2) +
  scale_fill_viridis_c(option = "D", direction = 1, na.value = "gray90") +
  coord_sf(xlim = c(xmin(extent_med), xmax(extent_med)),
           ylim = c(ymin(extent_med), ymax(extent_med)),
           expand = FALSE) +
  theme_minimal() +
  labs(title = "Governance rank", fill = "Rank")

save_png_svg(
  governance_med,
  file.path(dir_plots, "Mediterranean_governance.png"),
  file.path(dir_plots, "Mediterranean_governance.svg")
)

raster_med_governance <- rasterize(vect(med_governance_data), raster_med,
                                   field = "governance_pc_rank", fun = mean, na.rm = TRUE)

# Mask to land within the Mediterranean study extent.
raster_med_governance <- mask(raster_med_governance, vect(world_mask))


raster_df_governance <- as.data.frame(raster_med_governance, xy = TRUE, na.rm = FALSE)
colnames(raster_df_governance)[3] <- "governance_perc_rank"

write.csv(raster_df_governance, file.path(dir_output, "Governance_raster_data.csv"), row.names = FALSE)

#---------------------------------------------------------------------------------------------------
# Road density (log transformed; candidate indicator retained for sensitivity analysis)
#---------------------------------------------------------------------------------------------------
# Here I resample road density to the 0.5° template grid, log-transform, and export outputs.

road_density <- rast(file.path(dir_input, "grip4_total_dens_m_km2.asc"))

# Crop first, then only project if needed.
road_density_clipped <- crop(road_density, extent_med)

if (!grepl("4326", crs(road_density_clipped))) {
  road_density_clipped <- project(road_density_clipped, "EPSG:4326", method = "bilinear")
  road_density_clipped <- crop(road_density_clipped, extent_med)
}

road_density_resampled <- resample(road_density_clipped, raster_med, method = "bilinear")
road_density_resampled_log <- log(road_density_resampled + 1)

# Mask to land within the Mediterranean study extent.
road_density_resampled_log <- mask(road_density_resampled_log, vect(world_mask))


raster_df_road_density <- as.data.frame(road_density_resampled_log, xy = TRUE, na.rm = FALSE)
colnames(raster_df_road_density)[3] <- "road_density_log"

write.csv(raster_df_road_density, file.path(dir_output, "Road_density_raster_data.csv"), row.names = FALSE)

road_density_med <- ggplot(raster_df_road_density) +
  geom_raster(aes(x = x, y = y, fill = road_density_log)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D", direction = 1, na.value = "gray90") +
  coord_sf(xlim = c(xmin(extent_med), xmax(extent_med)),
           ylim = c(ymin(extent_med), ymax(extent_med)),
           expand = FALSE) +
  theme_minimal() +
  labs(title = "Road density", fill = "Log Road Density")

save_png_svg(
  road_density_med,
  file.path(dir_plots, "Mediterranean_road_density.png"),
  file.path(dir_plots, "Mediterranean_road_density.svg")
)

#---------------------------------------------------------------------------------------------------
# Poverty proxy: Sub-national HDI (SHDI) from Global Data Lab
#---------------------------------------------------------------------------------------------------
# Here I rasterise SHDI values from subnational polygons directly onto the 0.5° template grid.

shdi_csv <- file.path(dir_input, "GDL_Subnational_HDI-data.csv")
gdl_shp  <- file.path(dir_input, "GDL_shapefiles/GDL Shapefiles V6.6 large.shp")

shdi_tif <- file.path(dir_output, "poverty_proxy_SHDI_2022_0p5deg.tif")
shdi_out_csv <- file.path(dir_output, "poverty_proxy_SHDI_2022_raster_data_0p5deg.csv")
shdi_png <- file.path(dir_plots, "Mediterranean_poverty_proxy_SHDI_2022_0p5deg.png")
shdi_svg <- file.path(dir_plots, "Mediterranean_poverty_proxy_SHDI_2022_0p5deg.svg")

shdi <- read.csv(shdi_csv, stringsAsFactors = FALSE)
stopifnot("GDLCODE" %in% names(shdi))

shdi_2022 <- shdi %>%
  dplyr::select(GDLCODE, shdi_2022 = `X2022`) %>%
  mutate(GDLCODE = as.character(GDLCODE),
         shdi_2022 = suppressWarnings(as.numeric(shdi_2022)))

gdl_sf <- sf::st_read(gdl_shp, quiet = TRUE)
stopifnot("gdlcode" %in% names(gdl_sf))
gdl_sf$GDLCODE <- as.character(gdl_sf$gdlcode)

# Ensure polygons are in the same CRS as the template before rasterisation.
gdl_joined <- gdl_sf %>%
  left_join(shdi_2022, by = "GDLCODE") %>%
  st_transform(4326)

gdl_v <- terra::vect(gdl_joined)
shdi_raster <- terra::rasterize(gdl_v, raster_med, field = "shdi_2022")

shdi_raster_med <- terra::crop(shdi_raster, extent_med)
shdi_raster_med <- terra::mask(shdi_raster_med, terra::vect(world_mask))
shdi_raster_med <- clamp(shdi_raster_med, lower = 0, upper = 1, values = TRUE)
names(shdi_raster_med) <- "shdi_2022"

writeRaster(shdi_raster_med, shdi_tif, overwrite = TRUE)

raster_df_shdi <- as.data.frame(shdi_raster_med, xy = TRUE, na.rm = FALSE)
colnames(raster_df_shdi)[3] <- "shdi_2022"
write.csv(raster_df_shdi, shdi_out_csv, row.names = FALSE)

shdi_plot <- ggplot(raster_df_shdi) +
  geom_raster(aes(x = x, y = y, fill = shdi_2022)) +
  geom_sf(data = world_outline, fill = NA, color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "D", na.value = "gray90") +
  coord_sf(xlim = c(xmin(extent_med), xmax(extent_med)),
           ylim = c(ymin(extent_med), ymax(extent_med)),
           expand = FALSE) +
  theme_minimal() +
  labs(title = "Human Development Index",
       fill = "SHDI")

save_png_svg(shdi_plot, shdi_png, shdi_svg)

#===================================================================================================
# Multi-panel plotting — all current indicators in one figure
#===================================================================================================
# Here I pull together the ggplot objects created in the indicator blocks above and arrange them
# into a single multi-panel figure.
#
# The aridity index is deliberately excluded from the plot because it will be used later for clipping
# rather than as an archetype indicator. However, it is included in the final combined CSV below.

# Helper to remove legends/keys and reduce margins from each individual plot.
no_key <- function(p) {
  p +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
}

all_indicator_panels <- (
  no_key(iav_plot) + no_key(water_stress_plot) + no_key(pop_density_med)
) / (
  no_key(crop_med_plot) + no_key(irrigated_plot) + no_key(landcover_diversity_plot)
) / (
  no_key(graz_plot) + no_key(soc_plot) + no_key(governance_med)
) / (
  no_key(road_density_med) + no_key(shdi_plot) + no_key(ppt_plot)
)

save_png_svg(
  all_indicator_panels,
  file.path(dir_plots, "all_indicators_multipanel.png"),
  file.path(dir_plots, "all_indicators_multipanel.svg"),
  width  = 12,
  height = 8.5
)

#===================================================================================================
# FINAL: Build one combined dataframe on the shared 0.5° template grid
#===================================================================================================
# Join all candidate indicators on the shared 0.5° grid. The aridity index is included as a
# support variable so Script 2 can restrict the analysis to drylands (AI < 0.65).

csv_files <- c(
  "interannual_variability_raster_data.csv",
  "water_stress_raster_data.csv",
  "population_raster_data.csv",
  "crop_perc_raster_data.csv",
  "irrigated_perc_raster_data.csv",
  "landcover_diversity_raster_data.csv",
  "grazing_intensity_raster_data_0p5deg.csv",
  "soil_organic_carbon_stock_0_30cm_raster_data.csv",
  "Governance_raster_data.csv",
  "Road_density_raster_data.csv",
  "poverty_proxy_SHDI_2022_raster_data_0p5deg.csv",
  "mean_annual_precipitation_WorldClim_BIO12_raster_data.csv",
  "Aridity_Index_raster_data.csv"
)

# I read each CSV and join by (x, y). Each file should be aligned to the same 0.5° grid.
df_list <- csv_files %>%
  purrr::set_names(nm = tools::file_path_sans_ext(csv_files)) %>%
  purrr::map(~ read_csv(file.path(dir_output, .x)))


combined_df <- reduce(df_list, full_join, by = c("x", "y"))

glimpse(combined_df)

write_csv(
  combined_df,
  file.path(dir_output, "med_05deg_all_indicators_combined.csv")
)


#===================================================================================================
# SUPPORTING INFORMATION INDICATOR MAPS FROM EXISTING CSV
#===================================================================================================


#---------------------------------------------------------------------------------------------------
# Load existing processed data
#---------------------------------------------------------------------------------------------------

dat <- read_csv(
  "output_data/Script1/med_05deg_all_indicators_combined.csv",
  show_col_types = FALSE
)

# Restrict to the same dryland analysis domain used in Script 2
dat <- dat %>%
  filter(
    y >= 28,
    y <= 45,
    !is.na(aridity_index),
    aridity_index < 0.65
  )

#---------------------------------------------------------------------------------------------------
# Common formatting
#---------------------------------------------------------------------------------------------------

map_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.2),
    panel.grid = element_blank(),
    
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 2)
    ),
    
    legend.position = "right",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    legend.key.height = unit(1.0, "cm"),
    legend.key.width = unit(0.35, "cm"),
    
    # keep margins very small so rows sit close together
    plot.margin = margin(1, 2, 1, 2)
  )

map_coord <- coord_sf(
  xlim = c(xmin(extent_med), xmax(extent_med)),
  ylim = c(ymin(extent_med), ymax(extent_med)),
  expand = FALSE
)

add_boundaries <- function() {
  geom_sf(
    data = world_outline,
    fill = NA,
    colour = "grey30",
    linewidth = 0.2,
    inherit.aes = FALSE
  )
}

#===================================================================================================
# FIGURE: EXPOSURE INDICATORS
#===================================================================================================

p1 <- ggplot(dat, aes(x, y, fill = mean_annual_precipitation)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    direction = -1,
    name = "Annual precipitation\n(mm/year)",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Climatic water limitation")

p2 <- ggplot(dat, aes(x, y, fill = interannual_variability)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    limits = c(0, 5),
    name = "Aqueduct score\n(0–5)",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Water-supply variability")

p3 <- ggplot(dat, aes(x, y, fill = pop_density_log)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    name = "log(people/km² + 1)",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Population pressure")

exposure_fig <- p1 / p2 / p3 +
  plot_layout(ncol = 1, heights = c(1, 1, 1))

ggsave(
  "plots_script1/SI_exposure_indicators.png",
  exposure_fig,
  width = 5,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  "plots_script1/SI_exposure_indicators.svg",
  exposure_fig,
  width = 5,
  height = 7,
  bg = "white"
)



#===================================================================================================
# FIGURE: SENSITIVITY INDICATORS
#===================================================================================================

p4 <- ggplot(dat, aes(x, y, fill = cropland_percent)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    limits = c(0, 100),
    name = "Cropland extent\n(%)",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Agricultural land dependence")

p5 <- ggplot(dat, aes(x, y, fill = irrigated_percent)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = "Irrigated cropland\n(%)",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Irrigation dependence")

p6 <- ggplot(dat, aes(x, y, fill = livestock_density_km2)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = "Cattle + sheep + goats\nkm⁻²",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Livestock density")

p7 <- ggplot(dat, aes(x, y, fill = soc_0_30cm_kg_m2)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    name = "SOC stock\n(kg C m⁻²)",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Soil condition")

sensitivity_fig <- p4 / p5 / p6 / p7 +
  plot_layout(ncol = 1, heights = c(1, 1, 1, 1))

ggsave(
  "plots_script1/SI_sensitivity_indicators.png",
  sensitivity_fig,
  width = 6.8,
  height = 11.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  "plots_script1/SI_sensitivity_indicators.svg",
  sensitivity_fig,
  width = 6.8,
  height = 11.5,
  bg = "white"
)

#===================================================================================================
# FIGURE: ADAPTIVE CAPACITY INDICATORS
#===================================================================================================

p8 <- ggplot(dat, aes(x, y, fill = governance_perc_rank)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    limits = c(0, 100),
    name = "Governance percentile\nrank",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Institutional capacity")

p9 <- ggplot(dat, aes(x, y, fill = shdi_2022)) +
  geom_tile(width = 0.5, height = 0.5) +
  add_boundaries() +
  scale_fill_viridis_c(
    limits = c(0, 1),
    name = "SHDI",
    na.value = "white"
  ) +
  map_coord + map_theme +
  labs(title = "Socio-economic capacity")

adaptive_fig <- p8 / p9 +
  plot_layout(ncol = 1, heights = c(1, 1))

ggsave(
  "plots_script1/SI_adaptive_capacity_indicators.png",
  adaptive_fig,
  width = 6.8,
  height = 5.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  "plots_script1/SI_adaptive_capacity_indicators.svg",
  adaptive_fig,
  width = 6.8,
  height = 5.8,
  bg = "white"
)
