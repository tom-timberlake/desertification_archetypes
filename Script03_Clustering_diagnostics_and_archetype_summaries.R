####################################################################################################
# ARCHETYPE ANALYSIS SCRIPT 3 — CLUSTERING, DIAGNOSTICS AND ARCHETYPE SUMMARIES
#
# This script runs the final archetype analysis described in the manuscript and Supporting
# Information. It takes the processed 0–1 candidate indicators from Script 2 and:
#   1. Selects the nine indicators used in the main analysis
#   2. Evaluates candidate solutions from k = 2–10
#   3. Fits the final six-archetype solution using hierarchical initialisation followed by k-means
#   4. Produces the main archetype map and within-archetype indicator distributions
#   5. Summarises archetypes in original/interpretable units for Table 2
#   6. Assesses variable influence by indicator blinding
#   7. Tests robustness to alternative indicator sets
#
# The clustering procedure follows the Methods and Appendix S1. For each k, Ward.D2 hierarchical
# clustering is first applied to a random subset of 100 grid cells. Mean profiles of those groups are
# then used as starting centres for k-means clustering of the full dataset. A fixed random seed of
# 123 is used throughout for reproducibility.
#
# All analyses reported in the manuscript/SI are enabled by default below. Individual blocks can be
# switched off when re-running selected parts of the workflow.
####################################################################################################

rm(list = ls())

# ==================================================================================================
# PURPOSE
# ==================================================================================================
# This script runs the archetype clustering workflow with explicit run switches at the top.
# Set each run_* option to TRUE or FALSE depending on which parts of the workflow you want to run.


# ==================================================================================================

# ==================================================================================================
# PACKAGES
# ==================================================================================================

pkgs <- c(
  "readr", "dplyr", "tidyr", "purrr", "janitor", "tibble",
  "cluster", "ggplot2", "RColorBrewer", "mclust", "forcats",
  "rnaturalearth", "sf", "scales", "grid"
)

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

# ==================================================================================================
# USER SETTINGS
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# 1) RUN SWITCHES
# --------------------------------------------------------------------------------------------------
# These correspond to analyses reported in the manuscript or Supporting Information.

run_basic_k_diagnostics   <- TRUE   # Fig. S4; silhouette is controlled separately below (Fig. S5)
run_stability_analysis    <- TRUE   # Fig. S6
run_branching_diagram     <- TRUE   # Fig. S7
run_final_solution        <- TRUE   # final six-archetype solution
run_cluster_map           <- TRUE   # main archetype map
run_cluster_profiles      <- TRUE   # mean scaled profiles retained for interpretation
run_boxplots              <- TRUE   # within-archetype distributions
run_original_summaries    <- TRUE   # Table 2 and supporting summaries
run_variable_importance   <- TRUE   # Table S2 / Fig. S8
run_sensitivity_analysis  <- TRUE   # Table S1

# --------------------------------------------------------------------------------------------------
# 2) INPUT FILES
# --------------------------------------------------------------------------------------------------
# Script 2 provides the transformed/min-max scaled candidate indicators. Script 1 provides the
# corresponding original-unit values used for interpretation and Table 2.

in_with_xy      <- "output_data/Script2/indicators_transformed_scaled_with_xy.csv"
infile_original <- "output_data/Script1/med_05deg_all_indicators_combined.csv"

# --------------------------------------------------------------------------------------------------
# 3) OUTPUT FOLDERS
# --------------------------------------------------------------------------------------------------
dir.create("output_data/Script3", showWarnings = FALSE, recursive = TRUE)

plot_dir <- "plots_script3"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------------------------------
# 4) ALL AVAILABLE INDICATORS
# --------------------------------------------------------------------------------------------------
# All candidate indicators exported by Script 2. After preprocessing, 1 represents greater
# vulnerability/pressure/constraint and 0 represents the opposite.

all_available_indicators <- c(
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

# --------------------------------------------------------------------------------------------------
# 5) MANUAL INDICATOR SELECTION  *** DECISION POINT ***
# --------------------------------------------------------------------------------------------------
# TRUE identifies the nine indicators used for the main archetype solution. Excluded candidates
# are retained for the sensitivity analysis below.

indicator_keep <- c(
  mean_annual_precipitation = TRUE,
  interannual_variability   = TRUE,
  water_stress              = FALSE,
  pop_density_log           = TRUE,
  cropland_percent          = TRUE,
  irrigated_percent         = TRUE,
  landcover_diversity       = FALSE,
  grazing_intensity_log     = TRUE,
  soc_0_30cm_kg_m2          = TRUE,
  governance_perc_rank      = TRUE,
  road_density_log          = FALSE,
  shdi_2022                 = TRUE
)

selected_indicators <- names(indicator_keep)[indicator_keep]

cat("Indicators retained for main analysis:\n")
print(selected_indicators)

cat("\nIndicators excluded from main analysis:\n")
print(setdiff(all_available_indicators, selected_indicators))

# --------------------------------------------------------------------------------------------------
# 6) K SETTINGS  *** DECISION POINT ***
# --------------------------------------------------------------------------------------------------
# k_grid is only needed if diagnostics, stability or branching are being run.
# chosen_k is the final solution used for interpretation and downstream outputs.

k_grid   <- 2:10
chosen_k <- 6

# --------------------------------------------------------------------------------------------------
# 7) DIAGNOSTIC SETTINGS  *** DECISION POINT ***
# --------------------------------------------------------------------------------------------------
# Full silhouette calculation can be very slow because it requires a full pairwise distance matrix.
# Keep calculate_silhouette = FALSE unless you specifically need this diagnostic.
# If TRUE, the script uses an optional sample rather than the full dataset.

calculate_silhouette <- TRUE
silhouette_sample_n  <- 5000

# --------------------------------------------------------------------------------------------------
# 8) STABILITY SETTINGS  *** DECISION POINT ***
# --------------------------------------------------------------------------------------------------
# Increase B for a more robust analysis, but it will take longer.
# prop_subsample = fraction of rows used in each resampled run for stability assessment.

B_stability     <- 40
prop_subsample  <- 0.80
ward_subset_n   <- 100
random_seed     <- 123

# --------------------------------------------------------------------------------------------------
# 9) BRANCHING DIAGRAM SETTINGS  *** DECISION POINT ***
# --------------------------------------------------------------------------------------------------
# Use a short k sequence for the diagram to keep it readable.

branching_diagram_k_grid <- 2:6
branching_diagram_coloured_k <- 6

# --------------------------------------------------------------------------------------------------
# 10) OPTIONAL SENSITIVITY SCENARIOS  *** DECISION POINT ***
# --------------------------------------------------------------------------------------------------
# Define additional indicator sets to test robustness of the final clustering.

sensitivity_indicator_sets <- list(
  main = selected_indicators,
  all_indicators = all_available_indicators,
  no_governance = setdiff(selected_indicators, "governance_perc_rank"),
  no_irrigation = setdiff(selected_indicators, "irrigated_percent"),
  no_interannual_variability = setdiff(selected_indicators, "interannual_variability"),
  no_governance_irrigation = setdiff(
    selected_indicators,
    c("governance_perc_rank", "irrigated_percent")
  )
)

# ==================================================================================================
# OUTPUT FILE NAMES
# ==================================================================================================

out_dir <- "output_data/Script3"

out_features_csv          <- file.path(out_dir, "indicators_main_scaled.csv")
out_k_diagnostics_csv     <- file.path(out_dir, "k_diagnostics.csv")
out_stability_reps_csv    <- file.path(out_dir, "k_stability_replicates.csv")
out_stability_summary_csv <- file.path(out_dir, "k_stability_summary.csv")

out_clusters_xy           <- file.path(out_dir, "archetype_assignments_xy.csv")
out_xy_feats_k            <- file.path(out_dir, "indicator_data_scaled_with_archetypes.csv")
out_profiles              <- file.path(out_dir, "archetype_profiles_scaled_means.csv")
out_profiles_long         <- file.path(out_dir, "archetype_profiles_scaled_means_long.csv")
out_medians               <- file.path(out_dir, "archetype_medians_original_units.csv")
out_means_orig            <- file.path(out_dir, "archetype_means_original_units.csv")
out_cluster_sizes         <- file.path(out_dir, "archetype_sizes.csv")
out_variable_importance   <- file.path(out_dir, "indicator_blinding_importance.csv")
out_boxplot_stats         <- file.path(out_dir, "archetype_boxplot_stats.csv")
out_sensitivity           <- file.path(out_dir, "indicator_set_sensitivity.csv")
out_summary_table         <- file.path(out_dir, "archetype_summary_table_wide.csv")

# ==================================================================================================
# READABLE INDICATOR LABELS
# ==================================================================================================
# Used consistently in heatmaps, profile plots, boxplots and variable-importance plots.

indicator_titles <- c(
  mean_annual_precipitation = "Annual precipitation",
  interannual_variability = "Interannual rainfall variability",
  water_stress            = "Water stress",
  pop_density_log         = "Population density",
  cropland_percent        = "Cropland cover",
  irrigated_percent       = "Irrigated area",
  landcover_diversity     = "Land-cover diversity",
  grazing_intensity_log   = "Livestock density",
  soc_0_30cm_kg_m2        = "Soil organic carbon",
  governance_perc_rank    = "Governance",
  road_density_log        = "Road density",
  shdi_2022               = "Human Development Index"
)


# ==================================================================================================
# HELPER FUNCTIONS
# ==================================================================================================

check_feature_matrix <- function(df, name = "feature data") {
  if (sum(is.na(df)) > 0) {
    stop("NA values present in ", name, ".")
  }
  if (sum(!is.finite(as.matrix(df))) > 0) {
    stop("Inf/-Inf values present in ", name, ".")
  }
  invisible(TRUE)
}


ward_seeded_kmeans <- function(X, k, subset_n = 1000, iter.max = 200, seed = 123) {
  set.seed(seed)
  
  n <- nrow(X)
  p <- ncol(X)
  
  idx <- sample(seq_len(n), size = min(subset_n, n), replace = FALSE)
  X_sub <- X[idx, , drop = FALSE]
  
  hc <- hclust(dist(X_sub), method = "ward.D2")
  sub_cl <- cutree(hc, k = k)
  
  init_centers <- matrix(NA_real_, nrow = k, ncol = p)
  for (j in seq_len(k)) {
    members <- X_sub[sub_cl == j, , drop = FALSE]
    
    if (nrow(members) == 0) {
      init_centers[j, ] <- X_sub[sample(seq_len(nrow(X_sub)), 1), ]
    } else {
      init_centers[j, ] <- colMeans(members)
    }
  }
  
  colnames(init_centers) <- colnames(X)
  
  km <- kmeans(X, centers = init_centers, iter.max = iter.max)
  return(km)
}

mean_silhouette_sample <- function(X, cluster_vec, sample_n = 5000, seed = 123) {
  set.seed(seed)
  n <- nrow(X)
  idx <- sample(seq_len(n), size = min(sample_n, n), replace = FALSE)
  dist_obj <- dist(X[idx, , drop = FALSE])
  sil <- cluster::silhouette(cluster_vec[idx], dist_obj)
  mean(sil[, 3])
}

stable_points_fraction <- function(cl1, cl2, n_pairs = 50000, seed = 123) {
  stopifnot(length(cl1) == length(cl2))
  n <- length(cl1)
  if (n < 2) return(NA_real_)
  
  set.seed(seed)
  i1 <- sample.int(n, size = n_pairs, replace = TRUE)
  i2 <- sample.int(n, size = n_pairs, replace = TRUE)
  
  same1 <- cl1[i1] == cl1[i2]
  same2 <- cl2[i1] == cl2[i2]
  
  mean(same1 == same2)
}

stability_for_k <- function(X, k, B = 40, prop = 0.80, subset_n = 5000, seed = 123) {
  set.seed(seed)
  n <- nrow(X)
  out <- vector("list", B)
  
  for (b in seq_len(B)) {
    idx1 <- sort(sample(seq_len(n), size = floor(prop * n), replace = FALSE))
    idx2 <- sort(sample(seq_len(n), size = floor(prop * n), replace = FALSE))
    
    km1 <- ward_seeded_kmeans(X[idx1, , drop = FALSE], k = k, subset_n = subset_n, seed = seed + b)
    km2 <- ward_seeded_kmeans(X[idx2, , drop = FALSE], k = k, subset_n = subset_n, seed = seed + 1000 + b)
    
    overlap <- intersect(idx1, idx2)
    pos1 <- match(overlap, idx1)
    pos2 <- match(overlap, idx2)
    
    cl1 <- km1$cluster[pos1]
    cl2 <- km2$cluster[pos2]
    
    ari_val <- mclust::adjustedRandIndex(cl1, cl2)
    stable_val <- stable_points_fraction(
      cl1,
      cl2,
      n_pairs = min(50000, length(overlap)^2),
      seed = seed + 2000 + b
    )
    
    out[[b]] <- tibble(
      k = k,
      rep = b,
      n_overlap = length(overlap),
      adjusted_rand = ari_val,
      stable_fraction = stable_val
    )
  }
  
  bind_rows(out)
}

summarise_original_units <- function(dat_join, selected_vars) {
  cluster_sizes <- dat_join %>% count(cluster, name = "n_cells")
  
  medians_tbl <- dat_join %>%
    group_by(cluster) %>%
    summarise(across(all_of(selected_vars), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
    left_join(cluster_sizes, by = "cluster") %>%
    relocate(n_cells, .after = cluster)
  
  means_tbl <- dat_join %>%
    group_by(cluster) %>%
    summarise(across(all_of(selected_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    left_join(cluster_sizes, by = "cluster") %>%
    relocate(n_cells, .after = cluster)
  
  list(medians = medians_tbl, means = means_tbl, sizes = cluster_sizes)
}

grid_cell_area_km2 <- function(lat, cell_size_deg = 0.5, earth_radius_km = 6371) {
  lat_min <- lat - cell_size_deg / 2
  lat_max <- lat + cell_size_deg / 2
  lon_width_rad <- cell_size_deg * pi / 180

  earth_radius_km^2 * lon_width_rad *
    abs(sin(lat_max * pi / 180) - sin(lat_min * pi / 180))
}

indicator_blinding_importance <- function(X, chosen_k, subset_n = 5000, seed = 123) {
  full_km <- ward_seeded_kmeans(X, k = chosen_k, subset_n = subset_n, seed = seed)
  full_cl <- full_km$cluster
  
  purrr::map_dfr(seq_len(ncol(X)), function(j) {
    X_blind <- X
    X_blind[, j] <- mean(X_blind[, j], na.rm = TRUE)
    
    blind_km <- ward_seeded_kmeans(X_blind, k = chosen_k, subset_n = subset_n, seed = seed)
    blind_cl <- blind_km$cluster
    
    ari_val <- mclust::adjustedRandIndex(full_cl, blind_cl)
    
    tibble(
      indicator    = colnames(X)[j],
      blinded_ari  = ari_val,
      importance   = 1 - ari_val
    )
  }) %>%
    arrange(desc(importance))
}

compute_boxplot_stats <- function(df_long) {
  df_long %>%
    group_by(cluster, indicator) %>%
    summarise(
      mean   = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q25    = quantile(value, 0.25, na.rm = TRUE),
      q75    = quantile(value, 0.75, na.rm = TRUE),
      p05    = quantile(value, 0.05, na.rm = TRUE),
      p95    = quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

fit_solution_for_indicator_set <- function(data_with_xy, indicator_vec, k, subset_n = 1000, seed = 123) {
  feat <- data_with_xy %>% dplyr::select(all_of(indicator_vec))
  check_feature_matrix(feat, name = paste("feature matrix for", paste(indicator_vec, collapse = ", ")))
  
  
  X_tmp <- as.matrix(feat)
  km <- ward_seeded_kmeans(X_tmp, k = k, subset_n = subset_n, seed = seed)
  
  tibble(
    row_id = seq_len(nrow(data_with_xy)),
    cluster_num = km$cluster
  )
}

require_final_solution <- function() {
  if (!exists("cl", inherits = TRUE)) {
    stop(
      "This block requires the final solution object `cl`. ",
      "Set run_final_solution <- TRUE before running downstream output blocks."
    )
  }
}

# ==================================================================================================
# LOAD DATA
# ==================================================================================================
cat("\nLoading data...\n")

dryland_with_xy <- readr::read_csv(in_with_xy, show_col_types = FALSE)
orig_raw        <- readr::read_csv(infile_original, show_col_types = FALSE) %>% janitor::clean_names()

xy_cols <- c("x", "y")

missing_xy_main <- setdiff(xy_cols, names(dryland_with_xy))
if (length(missing_xy_main)) stop("Missing x/y columns in scaled input: ", paste(missing_xy_main, collapse = ", "))

missing_xy_orig <- setdiff(xy_cols, names(orig_raw))
if (length(missing_xy_orig)) stop("Missing x/y columns in original input: ", paste(missing_xy_orig, collapse = ", "))

missing_inds <- setdiff(selected_indicators, names(dryland_with_xy))
if (length(missing_inds)) stop("Missing selected indicators in scaled input: ", paste(missing_inds, collapse = ", "))

missing_inds_orig <- setdiff(selected_indicators, names(orig_raw))
if (length(missing_inds_orig)) {
  warning(
    "These selected indicators are missing from original-value file and will not appear in original-unit summaries: ",
    paste(missing_inds_orig, collapse = ", ")
  )
}

dryland_with_xy <- dryland_with_xy %>%
  mutate(across(all_of(xy_cols), as.double))

orig_raw <- orig_raw %>%
  mutate(across(all_of(xy_cols), as.double))

# ==================================================================================================
# BUILD FEATURE MATRIX FOR MAIN ANALYSIS
# ==================================================================================================
cat("\nBuilding feature matrix...\n")

features_only <- dryland_with_xy %>%
  dplyr::select(dplyr::all_of(selected_indicators))

check_feature_matrix(features_only, name = "main selected feature data")


readr::write_csv(features_only, out_features_csv)

X <- as.matrix(features_only)
set.seed(random_seed)

# ==================================================================================================
# STEP 1: BASIC DIAGNOSTICS ACROSS K
# ==================================================================================================
if (run_basic_k_diagnostics) {
  cat("\nRunning basic k diagnostics...\n")
  
  k_diag_tbl <- purrr::map_dfr(k_grid, function(k) {
    cat("  k = ", k, "\n", sep = "")
    km <- ward_seeded_kmeans(X, k = k, subset_n = ward_subset_n, seed = random_seed + k)
    
    out <- tibble(
      k = k,
      tot_withinss = km$tot.withinss,
      betweenss    = km$betweenss
    )
    
    if (calculate_silhouette) {
      out <- out %>%
        mutate(
          mean_silhouette = mean_silhouette_sample(
            X,
            km$cluster,
            sample_n = silhouette_sample_n,
            seed = random_seed + k
          )
        )
    }
    
    out
  })
  
  readr::write_csv(k_diag_tbl, out_k_diagnostics_csv)
  print(k_diag_tbl)
  
  p_elbow <- ggplot(k_diag_tbl, aes(x = k, y = tot_withinss)) +
    geom_line() +
    geom_point() +
    theme_minimal(base_size = 12) +
    labs(
      title = "Within-cluster variation across candidate solutions",
      x = "Number of clusters (k)",
      y = "Total within-cluster sum of squares"
    )
  
  ggsave(file.path(plot_dir, "SI_k_within_cluster_sum_of_squares.png"), p_elbow, width = 7, height = 5, dpi = 600, bg = "white")
  
  if (calculate_silhouette) {
    p_sil <- ggplot(k_diag_tbl, aes(x = k, y = mean_silhouette)) +
      geom_line() +
      geom_point() +
      theme_minimal(base_size = 12) +
      labs(
        title = "Cluster separation across candidate solutions",
        x = "Number of clusters (k)",
        y = "Mean silhouette width"
      )
    
    ggsave(file.path(plot_dir, "SI_k_mean_silhouette.png"), p_sil, width = 7, height = 5, dpi = 600, bg = "white")
  }
} else {
  cat("\nSkipping basic k diagnostics.\n")
}

# ==================================================================================================
# STEP 2: STABILITY / CONSISTENCY GRAPH ACROSS K
# ==================================================================================================
if (run_stability_analysis) {
  cat("\nRunning stability analysis across k...\n")
  
  stability_reps <- purrr::map_dfr(
    k_grid,
    ~ stability_for_k(
      X = X,
      k = .x,
      B = B_stability,
      prop = prop_subsample,
      subset_n = ward_subset_n,
      seed = random_seed + .x * 100
    )
  )
  
  stability_summary <- stability_reps %>%
    group_by(k) %>%
    summarise(
      mean_adjusted_rand = mean(adjusted_rand, na.rm = TRUE),
      sd_adjusted_rand   = sd(adjusted_rand, na.rm = TRUE),
      mean_stable_fraction = mean(stable_fraction, na.rm = TRUE),
      sd_stable_fraction   = sd(stable_fraction, na.rm = TRUE),
      .groups = "drop"
    )
  
  readr::write_csv(stability_reps, out_stability_reps_csv)
  readr::write_csv(stability_summary, out_stability_summary_csv)
  print(stability_summary)
  
  p_consistency_ari <- ggplot(stability_summary, aes(x = k, y = mean_adjusted_rand)) +
    geom_line() +
    geom_point() +
    geom_errorbar(
      aes(
        ymin = mean_adjusted_rand - sd_adjusted_rand,
        ymax = mean_adjusted_rand + sd_adjusted_rand
      ),
      width = 0.15
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Stability of candidate cluster solutions",
      x = "Number of clusters (k)",
      y = "Mean adjusted Rand index"
    )
  
  ggsave(file.path(plot_dir, "SI_k_stability_adjusted_rand.png"), p_consistency_ari,
         width = 7, height = 5, dpi = 600, bg = "white")
  
  p_consistency_stable <- ggplot(stability_summary, aes(x = k, y = mean_stable_fraction)) +
    geom_line() +
    geom_point() +
    geom_errorbar(
      aes(
        ymin = mean_stable_fraction - sd_stable_fraction,
        ymax = mean_stable_fraction + sd_stable_fraction
      ),
      width = 0.15
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Pairwise consistency across candidate cluster solutions",
      x = "Number of clusters (k)",
      y = "Mean stable-points fraction"
    )
  
  ggsave(file.path(plot_dir, "k_stability_pairwise_consistency.png"), p_consistency_stable,
         width = 7, height = 5, dpi = 600, bg = "white")
} else {
  cat("\nSkipping stability analysis.\n")
}

# ==================================================================================================
# STEP 3: BRANCHING STRUCTURE ACROSS K (FIG. S7)
# ==================================================================================================
if (run_branching_diagram) {
  cat("\nBuilding representative solutions for the branching diagram...\n")

  k_grid_branch <- branching_diagram_k_grid
  final_k <- branching_diagram_coloured_k

  solution_list <- purrr::map(
    k_grid_branch,
    ~ ward_seeded_kmeans(X, k = .x, subset_n = ward_subset_n, seed = random_seed + .x)
  )
  names(solution_list) <- paste0("k", k_grid_branch)

  membership_tbl <- purrr::map_dfc(seq_along(k_grid_branch), function(i) {
    tibble(!!paste0("k", k_grid_branch[i]) := solution_list[[i]]$cluster)
  })

  branching_edges <- purrr::map_dfr(seq_along(k_grid_branch)[-length(k_grid_branch)], function(i) {
    k_from <- k_grid_branch[i]
    k_to   <- k_grid_branch[i + 1]
    
    cl_from <- membership_tbl[[paste0("k", k_from)]]
    cl_to   <- membership_tbl[[paste0("k", k_to)]]
    
    tibble(
      cluster_from_num = cl_from,
      cluster_to_num   = cl_to
    ) %>%
      count(cluster_from_num, cluster_to_num, name = "n") %>%
      group_by(cluster_from_num) %>%
      mutate(prop_from = n / sum(n)) %>%
      ungroup() %>%
      mutate(
        k_from = k_from,
        k_to   = k_to,
        cluster_from = paste0("k", k_from, "_C", cluster_from_num),
        cluster_to   = paste0("k", k_to, "_C", cluster_to_num)
      )
  })
  
  node_tbl <- purrr::map_dfr(k_grid_branch, function(k) {
    cl_vec <- membership_tbl[[paste0("k", k)]]
    
    tibble(cluster_num = sort(unique(cl_vec))) %>%
      mutate(
        k = k,
        cluster_id = paste0("k", k, "_C", cluster_num)
      )
  })
  
  node_tbl <- node_tbl %>%
    group_by(k) %>%
    arrange(cluster_num, .by_group = TRUE) %>%
    mutate(order = row_number()) %>%
    ungroup()
  
  for (i in seq_along(k_grid_branch)[-1]) {
    k_prev <- k_grid_branch[i - 1]
    k_curr <- k_grid_branch[i]
    
    prev_order <- node_tbl %>%
      filter(k == k_prev) %>%
      dplyr::select(cluster_id, order)
    
    curr_edges <- branching_edges %>%
      filter(k_from == k_prev, k_to == k_curr) %>%
      left_join(prev_order, by = c("cluster_from" = "cluster_id")) %>%
      group_by(cluster_to) %>%
      summarise(
        barycenter = weighted.mean(order, w = prop_from, na.rm = TRUE),
        .groups = "drop"
      )
    
    curr_nodes <- node_tbl %>%
      filter(k == k_curr) %>%
      left_join(curr_edges, by = c("cluster_id" = "cluster_to")) %>%
      mutate(barycenter = ifelse(is.na(barycenter), Inf, barycenter)) %>%
      arrange(barycenter, cluster_num) %>%
      mutate(order = row_number()) %>%
      dplyr::select(cluster_id, order)
    
    node_tbl <- node_tbl %>%
      rows_update(curr_nodes, by = "cluster_id")
  }
  
  x_positions <- tibble(
    k = k_grid_branch,
    x = seq_along(k_grid_branch)
  )
  
  node_tbl <- node_tbl %>%
    left_join(x_positions, by = "k") %>%
    group_by(k) %>%
    mutate(
      n_k = n(),
      y = rev(seq_len(n_k)),
      y_plot = y
    ) %>%
    ungroup()
  
  edge_plot_tbl <- branching_edges %>%
    left_join(
      node_tbl %>% dplyr::select(cluster_id, x_from = x, y_from = y_plot),
      by = c("cluster_from" = "cluster_id")
    ) %>%
    left_join(
      node_tbl %>% dplyr::select(cluster_id, x_to = x, y_to = y_plot),
      by = c("cluster_to" = "cluster_id")
    ) %>%
    mutate(
      strength_cat = case_when(
        prop_from > 0.75 ~ "> 75 %",
        prop_from > 0.50 ~ "> 50 %",
        prop_from > 0.25 ~ "> 25 %",
        prop_from > 0.05 ~ "> 5 %",
        TRUE ~ NA_character_
      ),
      strength_cat = factor(strength_cat, levels = c("> 75 %", "> 50 %", "> 25 %", "> 5 %")),
      line_size = case_when(
        prop_from > 0.75 ~ 1.2,
        prop_from > 0.50 ~ 0.9,
        prop_from > 0.25 ~ 0.6,
        prop_from > 0.05 ~ 0.35,
        TRUE ~ 0.2
      ),
      line_alpha = case_when(
        prop_from > 0.75 ~ 1.0,
        prop_from > 0.50 ~ 0.8,
        prop_from > 0.25 ~ 0.5,
        prop_from > 0.05 ~ 0.25,
        TRUE ~ 0.0
      ),
      line_col = case_when(
        prop_from > 0.75 ~ "black",
        prop_from > 0.50 ~ "grey35",
        prop_from > 0.25 ~ "grey60",
        prop_from > 0.05 ~ "grey80",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(strength_cat))
  
  final_palette <- c(
    "#1f5aa6", "#b792d2", "#ef3b2c", "#f0c419",
    "#33a02c", "#9e9e9e", "#000000", "#1f9acb",
    "#fb9a99", "#6a3d9a", "#b15928", "#a6cee3"
  )
  
  node_tbl <- node_tbl %>%
    mutate(
      fill_col = ifelse(k == final_k, final_palette[cluster_num], "grey80")
    )
  
  arrow_spec <- grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  
  box_width  <- 0.34
  box_height <- 0.26
  
  legend_x_start <- max(x_positions$x) + 0.8
  legend_x_end   <- legend_x_start + 0.6
  legend_y_top   <- max(node_tbl$y_plot)
  
  legend_df <- tibble(
    x = legend_x_start,
    xend = legend_x_end,
    y = legend_y_top + c(-0.5, -1.2, -1.9, -2.6),
    yend = legend_y_top + c(-0.5, -1.2, -1.9, -2.6),
    label = c("> 75 %", "> 50 %", "> 25 %", "> 5 %"),
    line_size = c(1.2, 0.9, 0.6, 0.35),
    line_col = c("black", "grey35", "grey60", "grey80"),
    line_alpha = c(1.0, 0.8, 0.5, 0.25)
  )
  
  p_branching_kok <- ggplot() +
    geom_segment(
      data = edge_plot_tbl,
      aes(
        x = x_from + box_width / 2,
        y = y_from,
        xend = x_to - box_width / 2,
        yend = y_to
      ),
      linewidth = edge_plot_tbl$line_size,
      colour = edge_plot_tbl$line_col,
      alpha = edge_plot_tbl$line_alpha,
      arrow = arrow_spec,
      lineend = "round"
    ) +
    geom_rect(
      data = node_tbl,
      aes(
        xmin = x - box_width / 2,
        xmax = x + box_width / 2,
        ymin = y_plot - box_height / 2,
        ymax = y_plot + box_height / 2
      ),
      fill = node_tbl$fill_col,
      colour = "grey45",
      linewidth = 0.3
    ) +
    geom_text(
      data = x_positions,
      aes(x = x, y = max(node_tbl$y_plot) + 0.8, label = paste0("K=", k)),
      fontface = "bold",
      size = 4
    ) +
    geom_segment(
      data = legend_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = legend_df$line_size,
      colour = legend_df$line_col,
      alpha = legend_df$line_alpha,
      arrow = arrow_spec
    ) +
    geom_text(
      data = legend_df,
      aes(x = xend + 0.15, y = y, label = label),
      hjust = 0,
      size = 3.8
    ) +
    annotate(
      "text",
      x = legend_x_start,
      y = legend_y_top + 0.2,
      label = "Proportion from k",
      hjust = 0,
      fontface = "bold",
      size = 4
    ) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(
      breaks = x_positions$x,
      labels = NULL,
      limits = c(min(x_positions$x) - 0.3, legend_x_end + 1.5),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_y_continuous(
      breaks = NULL,
      expand = expansion(mult = c(0.08, 0.12))
    ) +
    theme_void() +
    theme(
      plot.margin = margin(15, 15, 15, 15),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    labs(title = "Branching structure of candidate cluster solutions")
  
  ggsave(
    filename = file.path(plot_dir, "SI_k_branching_structure.png"),
    plot = p_branching_kok,
    width = 6,
    height = 5,
    dpi = 600,
    bg = "white"
  )
} else {
  cat("\nSkipping branching diagram.\n")
}

# ==================================================================================================
# STEP 4: FINAL K CHECK
# ==================================================================================================
if (is.na(chosen_k)) {
  stop("chosen_k is NA. Review the diagnostics and set chosen_k manually, then rerun.")
}

if (!chosen_k %in% k_grid) {
  stop("chosen_k must be one of: ", paste(k_grid, collapse = ", "))
}

# ==================================================================================================
# STEP 5: FIT FINAL CHOSEN SOLUTION AND RELABEL AS A1-A6
# ==================================================================================================
if (run_final_solution) {
  cat("\nFitting final clustering solution for chosen_k = ", chosen_k, "...\n", sep = "")
  
  km_final <- ward_seeded_kmeans(X, k = chosen_k, subset_n = ward_subset_n, seed = random_seed)
  
  archetype_lookup <- tibble::tibble(
    cluster_num_raw = c(1, 5, 2, 6, 3, 4),
    archetype_num   = seq_len(chosen_k),
    cluster         = factor(paste0("A", archetype_num), levels = paste0("A", seq_len(chosen_k)))
  )
  
  cl <- dryland_with_xy %>%
    mutate(cluster_num_raw = km_final$cluster) %>%
    left_join(archetype_lookup, by = "cluster_num_raw")
  
  if (any(is.na(cl$cluster))) {
    stop("Some final cluster numbers were not found in archetype_lookup. Check the relabelling table.")
  }
  
  cl <- cl %>%
    mutate(
      cluster_num = archetype_num,
      cluster = factor(cluster, levels = paste0("A", seq_len(chosen_k)))
    )
  
  assignments_xy <- cl %>%
    dplyr::select(all_of(xy_cols), cluster, cluster_num, cluster_num_raw)
  
  readr::write_csv(assignments_xy, out_clusters_xy)
  readr::write_csv(cl, out_xy_feats_k)
} else {
  cat("\nSkipping final clustering solution.\n")
}

# ==================================================================================================
# STEP 6: MAP FINAL CLUSTERS
# ==================================================================================================
if (run_cluster_map) {
  require_final_solution()
  cat("\nMapping final clusters...\n")
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(4326)
  
  cluster_plot <- ggplot(cl, aes(x = x, y = y, fill = cluster)) +
    geom_tile(width = 0.5, height = 0.5) +
    geom_sf(data = world, color = "grey30", fill = NA, linewidth = 0.2, inherit.aes = FALSE) +
    coord_sf(xlim = c(-10, 40), ylim = c(28, 45), expand = FALSE) +
    scale_fill_brewer(palette = "Set3", drop = FALSE) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("Mediterranean archetype clusters (k = ", chosen_k, ")"),
      fill = "Archetype"
    )
  
  ggsave(
    filename = file.path(plot_dir, paste0("Mediterranean_clusters_k", chosen_k, ".png")),
    plot = cluster_plot,
    width = 8,
    height = 6,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0("Mediterranean_clusters_k", chosen_k, ".svg")),
    plot = cluster_plot,
    width = 8,
    height = 6,
    bg = "white"
  )
  
} else {
  cat("\nSkipping final cluster map.\n")
}

# ==================================================================================================
# STEP 7: MEAN SCALED ARCHETYPE PROFILES
# ==================================================================================================
# Mean scaled profiles are retained for interpretation of the archetypes, as described in the
# manuscript, but no additional profile figures are generated here.

if (run_cluster_profiles) {
  require_final_solution()
  cat("\nCalculating mean scaled archetype profiles...\n")

  profiles <- cl %>%
    group_by(cluster) %>%
    summarise(across(all_of(selected_indicators), mean), .groups = "drop") %>%
    arrange(cluster)

  readr::write_csv(profiles, out_profiles)

  profiles_long_tbl <- profiles %>%
    pivot_longer(-cluster, names_to = "indicator", values_to = "score") %>%
    mutate(
      indicator = factor(indicator, levels = selected_indicators),
      indicator_label = indicator_titles[as.character(indicator)]
    )

  readr::write_csv(profiles_long_tbl, out_profiles_long)

} else {
  cat("\nSkipping mean scaled archetype profiles.\n")
}

# ==================================================================================================
# STEP 8: BOXPLOTS / WITHIN-CLUSTER SPREAD
# ==================================================================================================
if (run_boxplots) {
  require_final_solution()
  cat("\nCreating boxplots...\n")
  
  # Order of indicators in the boxplot.
  boxplot_indicator_order <- c(
    "mean_annual_precipitation",
    "interannual_variability",
    "pop_density_log",
    "cropland_percent",
    "irrigated_percent",
    "grazing_intensity_log",
    "soc_0_30cm_kg_m2",
    "governance_perc_rank",
    "shdi_2022"
  )
  
  # Keep only indicators that are actually included in the selected clustering set.
  boxplot_indicator_order <- intersect(boxplot_indicator_order, selected_indicators)
  
  boxplot_long <- cl %>%
    dplyr::select(cluster, all_of(boxplot_indicator_order)) %>%
    pivot_longer(-cluster, names_to = "indicator", values_to = "value") %>%
    mutate(
      # For the boxplot display:
      # - annual precipitation is shown in its natural direction:
      #   0 = lowest precipitation, 1 = highest precipitation
      # - governance, SHDI and SOC are also shown in their natural direction:
      #   higher = more governance / human development / SOC
      value = dplyr::case_when(
        indicator == "mean_annual_precipitation" ~ 1 - value,
        indicator %in% c("governance_perc_rank", "shdi_2022", "soc_0_30cm_kg_m2") ~ 1 - value,
        TRUE ~ value
      ),
      indicator = factor(indicator, levels = boxplot_indicator_order),
      indicator_label = indicator_titles[as.character(indicator)],
      indicator_label = factor(
        indicator_label,
        levels = indicator_titles[boxplot_indicator_order]
      )
    )
  
  boxplot_stats <- compute_boxplot_stats(boxplot_long)
  readr::write_csv(boxplot_stats, out_boxplot_stats)
  
  p_box <- ggplot(boxplot_long, aes(x = cluster, y = value)) +
    geom_boxplot(outlier.size = 0.2) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~ indicator_label, scales = "fixed") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("Within-cluster distributions by indicator (k = ", chosen_k, ")"),
      x = "Archetype",
      y = "Scaled value"
    )
  
  ggsave(
    filename = file.path(plot_dir, "archetype_indicator_boxplots.png"),
    plot = p_box,
    width = 8,
    height = 8,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(plot_dir, "archetype_indicator_boxplots.svg"),
    plot = p_box,
    width = 8,
    height = 8,
    bg = "white"
  )
} else {
  cat("\nSkipping boxplots.\n")
}

# ==================================================================================================
# STEP 9: ORIGINAL-UNIT SUMMARIES
# ==================================================================================================
if (run_original_summaries) {
  require_final_solution()
  cat("\nCreating original-unit summaries...\n")
  
  if (!exists("assignments_xy", inherits = TRUE)) {
    assignments_xy <- cl %>%
      dplyr::select(all_of(xy_cols), cluster, cluster_num, cluster_num_raw)
  }
  
  orig_join <- orig_raw %>%
    inner_join(assignments_xy, by = setNames(xy_cols, xy_cols))
  
  # Standard original-unit summaries for the variables used directly in clustering.
  # This keeps grazing_intensity_log in the standard summary outputs, so the rest of the workflow
  # remains consistent with the selected_indicators vector.
  orig_selected_vars <- intersect(selected_indicators, names(orig_join))
  orig_summaries <- summarise_original_units(orig_join, orig_selected_vars)
  
  readr::write_csv(orig_summaries$medians, out_medians)
  readr::write_csv(orig_summaries$means, out_means_orig)
  readr::write_csv(orig_summaries$sizes, out_cluster_sizes)
  
  print(orig_summaries$medians)
  print(orig_summaries$means)
} else {
  cat("\nSkipping original-unit summaries.\n")
}

# ==================================================================================================
# STEP 10: SUMMARY TABLE FOR THE MANUSCRIPT
# ==================================================================================================
# This creates one wide table with archetypes as columns and summary metrics/indicator medians as rows.
#
# For the manuscript table, values are reported in interpretable original units:
# - population density is back-transformed from pop_density_log to people/km2
# - livestock density uses livestock_density_km2 from Script 1, which is cattle + sheep + goats per km2
# - other indicators are reported in the original units exported from Script 1

if (run_original_summaries) {
  require_final_solution()
  cat("\nCreating wide archetype summary table...\n")
  
  archetype_names <- c(
    A1 = "High-capacity agriculture with greater irrigation",
    A2 = "High-capacity rainfed systems",
    A3 = "Mid-capacity agriculture with greater irrigation",
    A4 = "Lower-capacity high-pressure rainfed systems",
    A5 = "Low-capacity arid rangelands",
    A6 = "Low-capacity hyper-arid rangelands"
  )
  
  if (!exists("assignments_xy", inherits = TRUE)) {
    assignments_xy <- cl %>%
      dplyr::select(all_of(xy_cols), cluster, cluster_num, cluster_num_raw)
  }
  
  orig_join <- orig_raw %>%
    inner_join(assignments_xy, by = setNames(xy_cols, xy_cols))
  
  # -----------------------------------------------------------------------------------------------
  # Indicator metadata
  # -----------------------------------------------------------------------------------------------
  # These are the variables and units used in the final paper summary table.
  # Note that grazing_intensity_log remains the clustering variable, but the table reports
  # livestock_density_km2, which was exported from Script 1 as cattle + sheep + goats per km2.
  
  indicator_summary_metadata <- tibble::tribble(
    ~indicator,                    ~row_label,                              ~units,
    "mean_annual_precipitation",   "Annual precipitation",                  "mm/year",
    "interannual_variability",     "Interannual rainfall variability",      "Aqueduct score, 0-5",
    "pop_density_people_km2",      "Population density",                    "people/km²",
    "cropland_percent",            "Cropland cover",                        "% of grid-cell land area",
    "irrigated_percent",           "Irrigated area",                        "% of grid-cell land area",
    "irrigated_share_cropland",    "Irrigated share of cropland",           "% of classified cropland",
    "livestock_density_km2",       "Livestock density",                     "cattle + sheep + goats per km²",
    "soc_0_30cm_kg_m2",            "Soil organic carbon",                   "kg C/m², 0-30 cm",
    "governance_perc_rank",        "Governance",                            "percentile rank, 0-100",
    "shdi_2022",                   "Human Development Index",               "index, 0-1"
  )
  
  # Keep only indicators included in the final selected clustering set.
  # Where the clustering variable is transformed, map it to a clearer original-scale table variable.
  selected_summary_indicators <- c(
    if ("mean_annual_precipitation" %in% selected_indicators) "mean_annual_precipitation",
    if ("interannual_variability" %in% selected_indicators) "interannual_variability",
    if ("pop_density_log" %in% selected_indicators) "pop_density_people_km2",
    if ("cropland_percent" %in% selected_indicators) "cropland_percent",
    if ("irrigated_percent" %in% selected_indicators) "irrigated_percent",
    if ("irrigated_percent" %in% selected_indicators) "irrigated_share_cropland",
    if ("grazing_intensity_log" %in% selected_indicators) "livestock_density_km2",
    if ("soc_0_30cm_kg_m2" %in% selected_indicators) "soc_0_30cm_kg_m2",
    if ("governance_perc_rank" %in% selected_indicators) "governance_perc_rank",
    if ("shdi_2022" %in% selected_indicators) "shdi_2022"
  )
  
  indicator_summary_metadata <- indicator_summary_metadata %>%
    filter(indicator %in% selected_summary_indicators) %>%
    mutate(indicator = factor(indicator, levels = selected_summary_indicators)) %>%
    arrange(indicator)
  
  # -----------------------------------------------------------------------------------------------
  # Unit conversions for table
  # -----------------------------------------------------------------------------------------------
  # pop_density_log was created in Script 1 as log(population density people/km2 + 1).
  # livestock_density_km2 is now created directly in Script 1 as cattle + sheep + goats per km2.
  
  orig_join_for_summary <- orig_join %>%
    mutate(
      pop_density_people_km2 = exp(pop_density_log) - 1,
      
      # Percentage of classified cropland that is irrigated.
      # Keep as NA where there is no classified cropland to avoid division by zero.
      irrigated_share_cropland = ifelse(
        cropland_percent > 0,
        100 * irrigated_percent / cropland_percent,
        NA_real_
      )
    )
  
  missing_summary_vars <- setdiff(selected_summary_indicators, names(orig_join_for_summary))
  
  if (length(missing_summary_vars)) {
    stop(
      "The following variables needed for the summary table are missing from the original-unit file: ",
      paste(missing_summary_vars, collapse = ", "),
      "\nIf livestock_density_km2 is missing, rerun Script 1 including the final combined CSV step."
    )
  }
  
  # -----------------------------------------------------------------------------------------------
  # Area percentage
  # -----------------------------------------------------------------------------------------------
  # Grid cells are equally weighted observations in the clustering, but 0.5° longitude-latitude
  # cells vary in physical area with latitude. Therefore, archetype area percentages are calculated
  # using the approximate physical area of each retained grid cell.
  
  area_pct_tbl <- orig_join_for_summary %>%
    mutate(cell_area_km2 = grid_cell_area_km2(y)) %>%
    group_by(cluster) %>%
    summarise(
      area_km2 = sum(cell_area_km2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      value = 100 * area_km2 / sum(area_km2)
    ) %>%
    select(cluster, value)
  
  # -----------------------------------------------------------------------------------------------
  # Population percentage
  # -----------------------------------------------------------------------------------------------
  # pop_density_people_km2 is people per km2. To estimate population per 0.5-degree grid cell, multiply
  # by approximate cell area in km2.
  
  population_pct_tbl <- orig_join_for_summary %>%
    mutate(
      cell_area_km2 = grid_cell_area_km2(y),
      population_est = pop_density_people_km2 * cell_area_km2
    ) %>%
    group_by(cluster) %>%
    summarise(population_est = sum(population_est, na.rm = TRUE), .groups = "drop") %>%
    mutate(value = 100 * population_est / sum(population_est, na.rm = TRUE)) %>%
    dplyr::select(cluster, value)
  
  # -----------------------------------------------------------------------------------------------
  # Median values in original/interpretable units
  # -----------------------------------------------------------------------------------------------
  
  indicator_medians_long <- orig_join_for_summary %>%
    group_by(cluster) %>%
    summarise(
      across(all_of(selected_summary_indicators), ~ median(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = -cluster,
      names_to = "indicator",
      values_to = "value"
    ) %>%
    left_join(indicator_summary_metadata, by = "indicator") %>%
    dplyr::select(row_label, units, cluster, value)
  
  # -----------------------------------------------------------------------------------------------
  # Combine metadata rows, percentage rows and indicator median rows
  # -----------------------------------------------------------------------------------------------
  
  summary_long <- dplyr::bind_rows(
    tibble(
      row_label = "Archetype name",
      units = "",
      cluster = factor(names(archetype_names), levels = paste0("A", seq_len(chosen_k))),
      value = as.character(unname(archetype_names))
    ),
    area_pct_tbl %>%
      mutate(
        row_label = "Area",
        units = "%",
        value = as.character(round(value, 2))
      ) %>%
      dplyr::select(row_label, units, cluster, value),
    population_pct_tbl %>%
      mutate(
        row_label = "Population",
        units = "%",
        value = as.character(round(value, 2))
      ) %>%
      dplyr::select(row_label, units, cluster, value),
    indicator_medians_long %>%
      mutate(
        value = as.character(round(value, 3))
      )
  )
  
  # Ensure rows appear in the requested order.
  summary_row_order <- c(
    "Archetype name",
    "Area",
    "Population",
    indicator_summary_metadata$row_label
  )
  
  archetype_summary_wide <- summary_long %>%
    mutate(
      row_label = factor(row_label, levels = summary_row_order),
      cluster = factor(cluster, levels = paste0("A", seq_len(chosen_k)))
    ) %>%
    arrange(row_label, cluster) %>%
    pivot_wider(
      names_from = cluster,
      values_from = value
    ) %>%
    mutate(row_label = as.character(row_label)) %>%
    rename(
      Variable = row_label,
      Units = units
    )
  
  readr::write_csv(
    archetype_summary_wide,
    out_summary_table
  )
  
  print(archetype_summary_wide)
}


# ==================================================================================================
# STEP 11: INDICATOR INFLUENCE (VARIABLE BLINDING; TABLE S2 / FIG. S8)
# ==================================================================================================
if (run_variable_importance) {
  cat("\nRunning indicator-blinding analysis...\n")
  
  importance_tbl <- indicator_blinding_importance(
    X = X,
    chosen_k = chosen_k,
    subset_n = ward_subset_n,
    seed = random_seed
  )
  
  readr::write_csv(importance_tbl, out_variable_importance)
  print(importance_tbl)
  
  importance_tbl <- importance_tbl %>%
    mutate(
      indicator_label = indicator_titles[indicator],
      indicator_label = factor(indicator_label, levels = indicator_titles[selected_indicators])
    )
  
  p_importance <- ggplot(
    importance_tbl,
    aes(x = reorder(indicator_label, importance), y = importance)
  ) +
    geom_col() +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("Influence of individual indicators (k = ", chosen_k, ")"),
      x = NULL,
      y = "Importance (cluster change after blinding)"
    )
  
  ggsave(
    filename = file.path(plot_dir, "indicator_blinding_importance.png"),
    plot = p_importance,
    width = 5,
    height = 5,
    dpi = 600,
    bg = "white"
  )
} else {
  cat("\nSkipping indicator-blinding analysis.\n")
}

# ==================================================================================================
# STEP 12: SENSITIVITY TO ALTERNATIVE INDICATOR SETS (TABLE S1)
# ==================================================================================================
if (run_sensitivity_analysis) {
  cat("\nRunning sensitivity analysis across alternative indicator sets...\n")
  
  main_solution <- fit_solution_for_indicator_set(
    data_with_xy = dryland_with_xy,
    indicator_vec = sensitivity_indicator_sets$main,
    k = chosen_k,
    subset_n = ward_subset_n,
    seed = random_seed
  )
  
  sensitivity_tbl <- purrr::imap_dfr(sensitivity_indicator_sets, function(ind_set, set_name) {
    sol <- fit_solution_for_indicator_set(
      data_with_xy = dryland_with_xy,
      indicator_vec = ind_set,
      k = chosen_k,
      subset_n = ward_subset_n,
      seed = random_seed
    )
    
    ari_val <- mclust::adjustedRandIndex(main_solution$cluster_num, sol$cluster_num)
    
    tibble(
      scenario = set_name,
      n_indicators = length(ind_set),
      adjusted_rand_vs_main = ari_val,
      indicators = paste(ind_set, collapse = "; ")
    )
  })
  
  readr::write_csv(sensitivity_tbl, out_sensitivity)
  print(sensitivity_tbl)
  
  p_sensitivity <- ggplot(sensitivity_tbl, aes(x = reorder(scenario, adjusted_rand_vs_main), y = adjusted_rand_vs_main)) +
    geom_col() +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("Sensitivity of chosen solution (k = ", chosen_k, ") to indicator set"),
      x = NULL,
      y = "Adjusted Rand vs main solution"
    )
  
  ggsave(
    filename = file.path(plot_dir, "indicator_set_sensitivity.png"),
    plot = p_sensitivity,
    width = 8,
    height = 5,
    dpi = 600,
    bg = "white"
  )
} else {
  cat("\nSkipping sensitivity analysis.\n")
}

# ==================================================================================================
# STEP 13: COMPLETION MESSAGE AND SESSION INFORMATION
# ==================================================================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "session_info.txt")
)

cat("\n================================================================================\n")
cat("ARCHETYPE ANALYSIS COMPLETE.\n")
cat("Blocks run:\n")
cat("- k diagnostics: ", run_basic_k_diagnostics, "\n", sep = "")
cat("- stability analysis: ", run_stability_analysis, "\n", sep = "")
cat("- branching diagram: ", run_branching_diagram, "\n", sep = "")
cat("- final solution: ", run_final_solution, "\n", sep = "")
cat("- cluster map: ", run_cluster_map, "\n", sep = "")
cat("- mean scaled profiles: ", run_cluster_profiles, "\n", sep = "")
cat("- within-archetype boxplots: ", run_boxplots, "\n", sep = "")
cat("- original-unit summaries: ", run_original_summaries, "\n", sep = "")
cat("- indicator blinding: ", run_variable_importance, "\n", sep = "")
cat("- indicator-set sensitivity: ", run_sensitivity_analysis, "\n", sep = "")
cat("================================================================================\n")
