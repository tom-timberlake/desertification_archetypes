# Archetypes of desertification vulnerability — analysis code

This repository contains the R code used to construct and analyse spatial archetypes of desertification vulnerability across Mediterranean drylands.

The workflow is organised into three scripts.

## Scripts

### `Script01_Importing_and_harmonising_indicator_layers_GitHub.R`

Documents the construction of the spatial indicator dataset from publicly available source layers.

The script:

- defines the Mediterranean study extent and common 0.5° grid;
- imports and harmonises candidate exposure, sensitivity and adaptive-capacity indicators;
- prepares the aridity index used to define the dryland analysis domain;
- exports individual processed layers and Supporting Information indicator maps; and
- combines all candidate indicators into a single grid-cell dataset.

The original geospatial source files are not included in this repository because several are large third-party datasets that are publicly available from their original providers.

Instead, the repository includes the combined output from Script 1:

`output_data/Script1/med_05deg_all_indicators_combined.csv`

This is the direct input to Script 2 and allows the full analytical workflow to be reproduced without downloading and reprocessing all source spatial datasets.

Users wishing to reproduce the complete data-construction stage should obtain the original datasets from the sources cited in the manuscript and retain the filenames and folder structure expected by Script 1.

### `Script02_Preparing_indicators_for_clustering_GitHub.R`

Prepares the combined grid-cell dataset for clustering.

The script:

- restricts the dataset to terrestrial dryland grid cells using aridity index < 0.65;
- checks missingness, distributions and pairwise correlations among candidate indicators;
- median-imputes missing values;
- applies the selected transformation to irrigated area;
- orients indicators so that higher processed values represent greater vulnerability, pressure or constraint;
- scales all candidate indicators to a common 0–1 range; and
- exports the processed candidate indicator dataset used by Script 3.

The main output is:

`output_data/Script2/indicators_transformed_scaled_with_xy.csv`

### `Script03_Clustering_diagnostics_and_archetype_summaries_GitHub.R`

Runs the final archetype analysis and produces the principal clustering diagnostics, robustness checks and manuscript summaries.

The script:

- selects the nine indicators used in the main analysis;
- compares candidate cluster solutions from k = 2–10;
- fits the final six-archetype solution using Ward.D2 hierarchical clustering to initialise k-means;
- uses a fixed random seed of 123 for reproducibility;
- evaluates cluster stability and branching structure;
- produces archetype assignments and within-archetype indicator summaries;
- generates original-unit summaries used for the manuscript table;
- assesses individual indicator influence using variable blinding; and
- tests robustness to alternative indicator sets.

The final six-archetype solution uses a hierarchical subset of 100 grid cells. Stability is assessed using 40 replicate comparisons of independent 80% subsamples. Mean silhouette width is calculated from a random sample of 5,000 grid cells.

## Running the analysis

To reproduce the analysis from the processed spatial dataset:

1. Place `med_05deg_all_indicators_combined.csv` in:
  
  `output_data/Script1/`
  
2. Run Script 2.
  
3. Run Script 3.
  

Script 1 only needs to be run if reconstructing the harmonised grid-cell dataset from the original public spatial sources.

The scripts create their required output and plot directories automatically.

## Data availability

The processed grid-cell dataset required to reproduce Scripts 2 and 3 is included with the archived code.

The original source layers used in Script 1 are not redistributed because they are publicly available third-party spatial datasets and, in several cases, are too large for practical inclusion in the repository.

Workshop data used to contextualise the archetypes are openly available from the University of York Research Data Repository:

https://doi.org/10.15124/623a29f4-f080-4d6e-9526-063f51ffae16

## Software

The analysis was conducted in R. Required packages are loaded within each script. Script 3 also writes `session_info.txt` to record the R and package versions used when the analysis is run.

## Citation

If using this code or processed dataset, please cite the associated manuscript:

**Timberlake et al.** *Archetypes of desertification vulnerability in Mediterranean drylands and their implications for socio-ecological resilience.*
