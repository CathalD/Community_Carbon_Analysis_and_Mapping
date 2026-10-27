# Invalidate only named generated files, preserving learner notes and inputs.
.stage_files <- list(
  field = file.path("outputs/part1", c("common_depth_stock.png", "core_stocks.csv", "cumulative_stock.png", "depth_profiles.png", "descriptive_statistics.csv", "layer_stocks.csv", "qc_flags.csv")),
  comparison = file.path("outputs/part1/comparison", c("R_session.txt", "community_matched_depth_values.csv", "community_vs_reference_groups.csv", "concentration_by_dataset.png", "density_by_dataset.png", "figure_reference_groups.csv", "reference_availability.csv", "reference_input_checksums.csv", "reference_summary_by_dataset.csv", "reference_summary_by_ecosystem.csv", "reference_summary_by_general_region.csv", "reference_summary_by_source_region.csv", "reference_summary_by_study_group.csv", "stock_by_dataset.png", "stock_by_ecosystem.png", "stock_by_region.png", "stock_by_source_region.png", "stock_by_study.png")),
  area = file.path("outputs/part2/area", c("location_map.png", "area_analysis.rds", "area_totals.csv", "ecosystem_map.png", "grid_values.csv", "group_stock_0_30.tif", "group_stock_map.png", "illustrative_ecosystem_cells.gpkg", "method_metadata.csv", "study_area.png")),
  spatial = file.path("outputs/part2/spatial", c("prior_comparison.csv", "distance_map.png", "grid_values.csv", "holdout_metrics.csv", "holdout_predictions.csv", "idw_stock_0_30.tif", "idw_stock_map.png", "method_metadata.csv", "method_totals.csv", "model_coefficients.csv", "model_stock_0_30.tif", "model_stock_map.png", "nearest_core_distance_m.tif", "observed_predicted.png", "simulated_elevation_m.tif", "simulated_wetness.tif", "simulation_training.csv", "spatial_analysis.rds")),
  bayes = file.path("outputs/part3", c("bayesian_analysis.rds", "bayesian_grid.csv", "bayesian_maps.png", "bayesian_summary.csv", "change_kg_m2.tif", "method_metadata.csv", "new_mean_kg_m2.tif", "new_sd_kg_m2.tif", "posterior_lower95_kg_m2.tif", "posterior_mean_kg_m2.tif", "posterior_sd_kg_m2.tif", "posterior_upper95_kg_m2.tif", "prior_mean_kg_m2.tif", "prior_sd_kg_m2.tif", "uncertainty.png")),
  kriging = file.path("outputs/part2/kriging", c("simulation.csv", "holdout.csv", "metrics.csv", "training_variogram.csv", "fitted_variogram.csv", "generating_variogram.csv", "grid_values.csv", "variogram.png", "kriging_maps.png", "kriging_mean_kg_m2.tif", "kriging_sd_kg_m2.tif", "kriging_analysis.rds", "method_metadata.csv"))
)
.stage_dependencies <- list(field=names(.stage_files), comparison="comparison",
 area=c("area","spatial","kriging","bayes"), spatial=c("spatial","bayes"),
 kriging="kriging", bayes="bayes")
.stage_marker <- "outputs/run_complete.txt"
invalidate_stage <- function(stage) {
  paths <- c(unlist(.stage_files[.stage_dependencies[[stage]]],use.names=FALSE),
    .stage_marker)
  unlink(paths)
  if(any(file.exists(paths))) stop("Close generated output files before rerunning.")
}
run_stage <- function(stage, code) {
  invalidate_stage(stage)
  tryCatch(eval(substitute(code), envir=parent.frame()), error=function(e) {
    invalidate_stage(stage)
    stop("Stage ",stage," failed; its outputs and dependent results were invalidated. ",
      conditionMessage(e), call.=FALSE)
  })
}
