#!/usr/bin/env Rscript

find_project_root <- function() {
  candidates <- unique(normalizePath(
    c(getwd(), dirname(dirname(commandArgs(trailingOnly = FALSE)[1]))),
    mustWork = FALSE
  ))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "config", "route_l.R"))) return(candidate)
  }
  stop("Run this script from the CommunityCarbonMap_v3 project root.")
}

project_root <- find_project_root()
setwd(project_root)
source("R/covariates/external_carbon_model.R")

inputs <- prepare_external_model_inputs(
  external_csv = file.path(
    project_root, "data", "derived", "covariate_analysis",
    "external_core_covariates_riley_canpeat.csv"
  ),
  community_csv = file.path(
    project_root, "data", "derived", "covariate_analysis",
    "community_core_covariates.csv"
  ),
  band_inventory_csv = file.path(
    project_root, "outputs", "covariate_analysis", "part2_qa",
    "band_inventory.csv"
  ),
  block_size_m = 100000,
  seed = 470316L
)
output_dir <- file.path(
  project_root, "outputs", "covariate_analysis", "part5_external_model"
)
write_external_model_preparation(inputs, output_dir)

cv <- run_nested_external_cv(
  inputs,
  output_dir = output_dir,
  inner_folds = 4L,
  inner_trees = 300L,
  outer_trees = 1000L,
  seed = 470316L,
  num_threads = 2L
)
final <- fit_final_external_models(
  inputs,
  output_dir = output_dir,
  inner_trees = 500L,
  final_trees = 1500L,
  seed = 470316L,
  num_threads = 2L
)
plot_external_model_results(cv, final, output_dir)

cat("\nExternal cross-validation metrics\n")
print(cv$metrics[cv$metrics$source_scope == "all", ], row.names = FALSE)
cat("\nCommunity final transfer metrics\n")
print(final$community_metrics, row.names = FALSE)
message("Part 5 outputs written to outputs/covariate_analysis/part5_external_model")

