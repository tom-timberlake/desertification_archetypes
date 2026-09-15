####################################################################################################
# ARCHETYPE ANALYSIS SCRIPT 2 — PREPARING INDICATORS FOR CLUSTERING
#
# This script prepares the candidate indicator set generated in Script 1 for the archetype analysis
# described in the manuscript.
#
# The workflow follows the Methods and Supporting Information:
#   1. Restrict the analysis to terrestrial dryland grid cells (aridity index < 0.65)
#   2. Check missingness, distributions and pairwise correlations among candidate indicators
#   3. Median-impute missing indicator values
#   4. Apply the selected transformation to irrigated area
#   5. Orient indicators so that higher processed values consistently indicate greater
#      vulnerability, pressure or constraint
#   6. Min-max scale all candidate indicators to 0–1
#   7. Export the processed candidate set for indicator selection and clustering in Script 3
#
# Population density, livestock density and road density were log-transformed in Script 1.
# Irrigated area is Yeo–Johnson transformed here because of its strongly right-skewed distribution.
# Aridity index is retained only to define the dryland analysis domain and is not used as a
# clustering indicator.
#
# Script 2 prepares all candidate indicators. The final nine-indicator set used in the main analysis
# is selected explicitly in Script 3, where alternative indicator sets are also tested.
####################################################################################################


# NOTE ON INPUT DATA
# --------------------------------------------------------------------------------------------------
# This script can be run independently of the large raw spatial datasets used in Script 1.
# The repository includes the analysis-ready output from Script 1:
#
#   output_data/Script1/med_05deg_all_indicators_combined.csv
#
# This contains the harmonised 0.5° candidate indicator dataset and is the starting point for all
# subsequent preprocessing and clustering analyses.


rm(list = ls())

#===================================================================================================
# 1. PACKAGES
#===================================================================================================

pkgs <- c(
  "readr", "dplyr", "tidyr", "purrr", "janitor",
  "moments", "recipes", "ggplot2", "tibble"
)

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

#===================================================================================================
# 2. INPUTS / OUTPUTS
#===================================================================================================

infile <- "output_data/Script1/med_05deg_all_indicators_combined.csv"

out_dir  <- "output_data/Script2"
plot_dir <- "plots_script2"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Main files used by Script 3.
out_features_csv <- file.path(out_dir, "indicators_transformed_scaled.csv")
out_xy_features  <- file.path(out_dir, "indicators_transformed_scaled_with_xy.csv")

# Original values are retained for interpretation and checking.
out_xy_original <- file.path(out_dir, "Indicator_data_clean_xy_original.csv")

#===================================================================================================
# 3. DATA CLEANING AND DIAGNOSTICS
#===================================================================================================

# Keep x/y as ID columns throughout.
xy_cols <- c("x", "y")

# Candidate indicator columns produced by Script 1.
# Aridity_index is deliberately not included here because it is used as a filtering/clipping variable.
indicator_cols <- c(
  "interannual_variability", # higher = more variable water supply / more vulnerable
  "water_stress",            # higher = greater baseline water stress / more vulnerable
  "mean_annual_precipitation", # higher = wetter conditions / less vulnerable; reversed later
  "pop_density_log",         # already log-transformed; higher = more pressure/vulnerability
  "cropland_percent",        # higher = more agricultural pressure/vulnerability
  "irrigated_percent",       # higher = more irrigation dependence / water-demand pressure
  "landcover_diversity",     # higher = more diverse landscape; reversed later
  "grazing_intensity_log",   # already log-transformed; higher = more grazing pressure/vulnerability
  "soc_0_30cm_kg_m2",        # higher = greater SOC stock; reversed later
  "governance_perc_rank",    # higher = less vulnerable; reversed later
  "road_density_log",        # higher = less isolated; reversed later
  "shdi_2022"                # higher = less vulnerable; reversed later
)

# Support variable used to define the dryland domain; not used as a clustering indicator.
support_cols <- c(
  "aridity_index"
)

#---------------------------------------------------------------------------------------------------
# 3.1 Load combined data
#---------------------------------------------------------------------------------------------------

df <- readr::read_csv(infile, show_col_types = FALSE) %>%
  janitor::clean_names()

# Hard spatial safety filter.
# This guarantees that no pixels below 28°N are included, even if an upstream layer
# accidentally retained values outside the intended study extent.
df <- df %>%
  filter(y >= 28, y <= 45)

cat(
  "Latitude range after hard spatial filter:",
  min(df$y, na.rm = TRUE), "to", max(df$y, na.rm = TRUE), "\n"
)

# Check that expected columns are present.
missing_expected <- setdiff(c(xy_cols, indicator_cols, support_cols), names(df))

if (length(missing_expected)) {
  stop(
    "These expected columns are missing from the combined CSV:\n",
    paste0(" - ", missing_expected, collapse = "\n"),
    "\n\nTip: check Script 1 output filenames and column names."
  )
}

# Keep only x/y, current indicators, and support variables.
df <- df %>%
  dplyr::select(all_of(c(xy_cols, indicator_cols, support_cols)))

#---------------------------------------------------------------------------------------------------
# 3.2 Remove non-land pixels
#---------------------------------------------------------------------------------------------------
# Some indicators can contain zeros over sea rather than NA.
# To avoid retaining ocean pixels, define land using indicators that are reliably NA over water.

land_check_cols <- c(
  "cropland_percent",
  "landcover_diversity",
  "soc_0_30cm_kg_m2",
  "shdi_2022"
)

# Keep rows where at least one land-defining indicator is non-NA.
land <- df %>%
  filter(if_any(all_of(land_check_cols), ~ !is.na(.x)))

cat("Kept land rows:", nrow(land), "of", nrow(df), "\n")

#---------------------------------------------------------------------------------------------------
# 3.3 Restrict the analysis to dryland pixels
#---------------------------------------------------------------------------------------------------
# Restrict the analysis to drylands using the aridity-index threshold described in the manuscript.
# The aridity index is used only to define the analysis domain.

dryland <- land %>%
  filter(!is.na(aridity_index) & aridity_index < 0.65)

cat("Kept dryland rows:", nrow(dryland), "of", nrow(land), "\n")

# Keep an untouched copy of original values for interpretation and manuscript summaries.
dryland_original <- dryland

# Use x/y plus candidate indicators for preprocessing.
dryland_analysis <- dryland %>%
  dplyr::select(all_of(c(xy_cols, indicator_cols)))

#---------------------------------------------------------------------------------------------------
# 3.4 Missingness diagnostics
#---------------------------------------------------------------------------------------------------

missing_tbl <- dryland %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "na_count") %>%
  mutate(
    total_rows = nrow(dryland),
    na_pct = round(100 * na_count / total_rows, 2)
  ) %>%
  arrange(desc(na_pct))

missing_tbl

#---------------------------------------------------------------------------------------------------
# 3.5 Distribution diagnostics
#---------------------------------------------------------------------------------------------------
# These summaries were used to identify strongly skewed indicators requiring transformation.
# Normality is not assumed by the clustering procedure, so I report skewness and kurtosis rather
# than formal normality tests.

num_df <- dryland_analysis %>%
  dplyr::select(-all_of(xy_cols)) %>%
  dplyr::select(where(is.numeric))

distribution_tbl <- purrr::map_dfr(names(num_df), function(v) {
  x <- num_df[[v]]
  x <- x[is.finite(x)]

  tibble(
    variable = v,
    n_non_na = length(x),
    skewness = if (length(x) >= 3 && sd(x) > 0) moments::skewness(x) else NA_real_,
    kurtosis_excess = if (length(x) >= 4 && sd(x) > 0) moments::kurtosis(x) - 3 else NA_real_
  )
})

print(distribution_tbl)

readr::write_csv(
  missing_tbl,
  file.path(out_dir, "indicator_missingness.csv")
)

readr::write_csv(
  distribution_tbl,
  file.path(out_dir, "indicator_distribution_diagnostics.csv")
)

#===================================================================================================
# 4. PREPROCESSING — READY FOR CLUSTERING
#===================================================================================================
# The preprocessing choices here correspond to the Methods and Appendix S1.
#
# Transformation strategy:
# - irrigated_percent: Yeo–Johnson transformation because of strong right-skew
# - pop_density_log, grazing_intensity_log and road_density_log were already log-transformed
#   in Script 1
# - remaining indicators are left untransformed
#
# Workflow:
# - median-impute missing indicator values
# - apply the irrigated-area transformation
# - remove zero-variance predictors
# - orient all indicators so that higher values indicate greater vulnerability, pressure or constraint
# - min-max scale each indicator to the range 0–1
#
# A previous working version z-standardised indicators before min-max scaling. That step was removed
# because it is mathematically redundant: min-max scaling the z-scores gives the same final 0–1
# values as min-max scaling the directionally aligned variables directly.

rec <- recipe(~ ., data = dryland_analysis) %>%
  update_role(all_of(xy_cols), new_role = "id") %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_YeoJohnson(irrigated_percent) %>%
  step_zv(all_predictors())

prepped <- prep(rec)
dryland_done <- bake(prepped, new_data = NULL)

#---------------------------------------------------------------------------------------------------
# 4.1 Directional alignment — ensure all indicators point toward vulnerability
#---------------------------------------------------------------------------------------------------
# After this step:
# - higher values consistently represent higher desertification vulnerability

dryland_vuln <- dryland_done %>%
  mutate(
    # already vulnerability-oriented pressures
    interannual_variability = interannual_variability,
    water_stress = water_stress,
    pop_density_log = pop_density_log,
    cropland_percent = cropland_percent,
    irrigated_percent = irrigated_percent,
    grazing_intensity_log = grazing_intensity_log,
    
    # reverse protective / resilience / capacity indicators
    landcover_diversity = -landcover_diversity,      # lower diversity = more vulnerable
    mean_annual_precipitation = -mean_annual_precipitation, # lower precipitation = more vulnerable
    soc_0_30cm_kg_m2 = -soc_0_30cm_kg_m2,            # lower SOC stock = more vulnerable
    governance_perc_rank = -governance_perc_rank,    # weaker governance = more vulnerable
    road_density_log = -road_density_log,            # lower accessibility = more isolated/vulnerable
    shdi_2022 = -shdi_2022                           # lower human development = more vulnerable
  )

#---------------------------------------------------------------------------------------------------
# 4.2 Min-max scale all candidate indicators to 0–1
#---------------------------------------------------------------------------------------------------

minmax_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

features_only <- dryland_vuln %>%
  dplyr::select(-all_of(xy_cols)) %>%
  mutate(across(everything(), minmax_01))

#---------------------------------------------------------------------------------------------------
# 4.3 Export processed and original-value datasets
#---------------------------------------------------------------------------------------------------

# Final 0–1 candidate indicator matrix with coordinates.
dryland_with_xy <- bind_cols(
  dryland_vuln %>% dplyr::select(all_of(xy_cols)),
  features_only
)

readr::write_csv(features_only, out_features_csv)
readr::write_csv(dryland_with_xy, out_xy_features)

# Original values are retained for interpretation and downstream summaries.
readr::write_csv(
  dryland_original %>%
    dplyr::select(all_of(c(xy_cols, indicator_cols, support_cols))),
  out_xy_original
)

# Basic checks before clustering.
na_total_scaled <- sum(is.na(features_only))
range01_ok <- all(
  vapply(
    features_only,
    function(x) min(x, na.rm = TRUE) >= 0 && max(x, na.rm = TRUE) <= 1,
    logical(1)
  )
)

cat(
  "\n=== Preprocessing checks ===\n",
  "NAs remaining in scaled indicators: ", na_total_scaled, "\n",
  "All indicators within [0,1]: ", range01_ok, "\n",
  sep = ""
)

#===================================================================================================
# 5. CORRELATION DIAGNOSTICS
#===================================================================================================
# Pairwise correlations are calculated across the processed candidate indicators. Because the
# transformations used for directional alignment and min-max scaling are monotonic linear operations
# (apart from the Yeo–Johnson transformation already applied to irrigation), these diagnostics are
# used to identify strongly overlapping candidate indicators before the final set is chosen in Script 3.

#---------------------------------------------------------------------------------------------------
# 5.1 Correlation matrix
#---------------------------------------------------------------------------------------------------

cor_mat <- cor(features_only, method = "pearson")

# Mean absolute off-diagonal correlation
mean_abs_cor <- cor_mat |>
  as.data.frame() |>
  tibble::rownames_to_column("var1") |>
  tidyr::pivot_longer(-var1, names_to = "var2", values_to = "r") |>
  dplyr::filter(var1 != var2) |>
  dplyr::mutate(abs_r = abs(r)) |>
  dplyr::summarise(mean_abs_r = mean(abs_r)) |>
  dplyr::pull(mean_abs_r)

cat(sprintf("\nMean absolute Pearson correlation (off-diagonal): %.2f\n", mean_abs_cor))

#---------------------------------------------------------------------------------------------------
# 5.2 Identify strongest pairwise correlations
#---------------------------------------------------------------------------------------------------

cor_pairs <- cor_mat |>
  as.data.frame() |>
  tibble::rownames_to_column("var1") |>
  tidyr::pivot_longer(-var1, names_to = "var2", values_to = "r") |>
  dplyr::filter(var1 < var2) |>
  dplyr::mutate(abs_r = abs(r)) |>
  dplyr::arrange(desc(abs_r))

cat("\nTop correlated indicator pairs (by |r|):\n")
print(head(cor_pairs, 10), row.names = FALSE)

high_cor_threshold <- 0.60

cat(sprintf("\nIndicator pairs with |r| ≥ %.2f:\n", high_cor_threshold))
print(
  dplyr::filter(cor_pairs, abs_r >= high_cor_threshold),
  row.names = FALSE
)

#---------------------------------------------------------------------------------------------------
# 5.3 Correlation matrix plot
#---------------------------------------------------------------------------------------------------

cor_long <- cor_mat |>
  as.data.frame() |>
  tibble::rownames_to_column("var1") |>
  tidyr::pivot_longer(-var1, names_to = "var2", values_to = "r")

# Order variables by hierarchical clustering on correlation structure.
var_order <- hclust(dist(1 - cor_mat))$order
ordered_vars <- colnames(cor_mat)[var_order]

cor_long$var1 <- factor(cor_long$var1, levels = ordered_vars)
cor_long$var2 <- factor(cor_long$var2, levels = ordered_vars)

# Keep lower triangle including diagonal.
cor_long_tri <- cor_long %>%
  mutate(
    var1 = factor(var1, levels = unique(var1)),
    var2 = factor(var2, levels = unique(var1))
  ) %>%
  filter(as.numeric(var1) >= as.numeric(var2))

# Dataset for diagonal tiles.
cor_diag <- cor_long_tri %>%
  filter(var1 == var2)

cor_plot <- ggplot() +
  
  # Triangle correlations, excluding diagonal.
  geom_tile(
    data = cor_long_tri %>% filter(var1 != var2),
    aes(x = var1, y = var2, fill = r),
    color = "white",
    linewidth = 0.2
  ) +
  
  # Grey diagonal.
  geom_tile(
    data = cor_diag,
    aes(x = var1, y = var2),
    fill = "grey85",
    color = "white",
    linewidth = 0.2
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "",
    x = NULL,
    y = NULL
  )

print(cor_plot)

ggsave(
  filename = file.path(plot_dir, "indicator_correlation_matrix.png"),
  plot = cor_plot,
  width = 6,
  height = 6,
  dpi = 600,
  bg = "white"
)

#---------------------------------------------------------------------------------------------------
# 5.4 Interpretation of correlation diagnostics
#---------------------------------------------------------------------------------------------------
# These diagnostics informed the final indicator selection reported in the manuscript. Land-cover
# diversity and road density showed substantial overlap with other candidate indicators, whereas
# livestock density was retained despite some correlation with population density because it represents
# a distinct land-use pressure. The main indicator set itself is specified explicitly in Script 3.

#===================================================================================================
# 6. FINAL CANDIDATE SET FOR SCRIPT 3
#===================================================================================================
# Script 3 makes the final inclusion/exclusion decisions. All 12 processed candidate indicators are
# therefore retained here so that the main analysis and alternative indicator-set sensitivity tests
# use exactly the same preprocessing.

final_candidate_set <- c(
  "mean_annual_precipitation",
  "interannual_variability",
  "water_stress",
  "pop_density_log",
  "cropland_percent",
  "irrigated_percent",
  "landcover_diversity",
  "grazing_intensity_log",
  "soc_0_30cm_kg_m2",
  "governance_perc_rank",
  "road_density_log",
  "shdi_2022"
)

missing_final_vars <- setdiff(final_candidate_set, names(dryland_with_xy))

if (length(missing_final_vars)) {
  stop(
    "These candidate indicators are missing from the processed dataset:\n",
    paste0(" - ", missing_final_vars, collapse = "\n")
  )
}

# Reorder explicitly before the final export used by Script 3.
features_only_final <- features_only %>%
  dplyr::select(all_of(final_candidate_set))

dryland_with_xy_final <- dryland_with_xy %>%
  dplyr::select(all_of(c(xy_cols, final_candidate_set)))

readr::write_csv(
  features_only_final,
  file.path(out_dir, "indicators_transformed_scaled.csv")
)

readr::write_csv(
  dryland_with_xy_final,
  file.path(out_dir, "indicators_transformed_scaled_with_xy.csv")
)

# Save correlation summaries used to document the indicator-selection checks.
readr::write_csv(
  cor_pairs,
  file.path(out_dir, "indicator_correlation_pairs_sorted.csv")
)

high_pairs <- cor_pairs %>%
  dplyr::filter(abs_r >= high_cor_threshold)

readr::write_csv(
  high_pairs,
  file.path(out_dir, "indicator_high_correlation_pairs.csv")
)

cat(
  "\nSaved processed candidate indicators for Script 3:\n",
  "- ", file.path(out_dir, "indicators_transformed_scaled.csv"), "\n",
  "- ", file.path(out_dir, "indicators_transformed_scaled_with_xy.csv"), "\n",
  "- ", out_xy_original, "\n",
  sep = ""
)

cat("\nScript 2 complete. Final indicator selection and clustering are carried out in Script 3.\n")
