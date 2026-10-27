# The short route: field measurements -> stocks -> sample map and IDW map.
run_course <- function(workbook_path, aoi_path, output_dir = "outputs/my_data",
  density_override_J7 = NA_real_, data_origin = "User data",
  grid_cellsize_m = NULL) {
  if (!file.exists("Community_Carbon_Analysis_and_Mapping.Rproj")) {
    stop("Open Community_Carbon_Analysis_and_Mapping.Rproj, then run this code in its Console.")
  }
  # Clear only our named results, so a failed rerun cannot leave an old map
  # beside new measurements. Preserve input files and anything the user adds.
  generated <- c(file.path(output_dir, "input_files.csv"),
    file.path(output_dir, "part1", c("layer_stocks.csv", "core_stocks.csv",
      "descriptive_statistics.csv", "qc_flags.csv", "depth_profiles.png",
      "cumulative_stock.png", "common_depth_stock.png")),
    file.path(output_dir, "maps", c("location_map.png", "sample_map.png", "idw_stock_map.png",
      "idw_stock_0_30.tif", "grid_values.csv", "mapping_metadata.csv", "excluded_cores.csv")))
  # The short route may share its output directory with the optional local route.
  # Invalidate those downstream products too, including after a failed field import.
  source("R/pipeline.R", local=TRUE)
  generated <- c(generated,
    file.path(output_dir, "part1/comparison", basename(.stage_files$comparison)),
    file.path(output_dir, c("local_grid_template.csv", "local_grid.gpkg", "local_stage_status.csv",
      "local_area_totals.csv", "local_group_map.png", "local_model_map.png", "local_holdout.csv",
      "local_model_metrics.csv", "local_model_coefficients.csv", "local_bayesian.csv", "local_bayesian_map.png")))
  input_paths <- normalizePath(c(workbook_path, aoi_path), mustWork = FALSE)
  if (any(normalizePath(generated, mustWork = FALSE) %in% input_paths)) {
    stop("Keep input files outside the generated output paths. Use data/input/ and data/aoi/.")
  }
  unlink(generated)
  if (any(file.exists(generated))) stop("Close the old output files so this run can replace them.")
  if (!file.exists(workbook_path)) {
    stop("Field-data file not found: ", workbook_path,
      ". Copy the template to data/input/my_field_data.csv, fill it in, and save it.")
  }
  if (!file.exists(aoi_path)) {
    stop("Area boundary not found: ", aoi_path,
      ". Save your polygon as data/aoi/my_aoi.geojson. See docs/start_here.md.")
  }
  source("R/field.R", local = TRUE)
  source("R/spatial.R", local = TRUE)
  input <- import_core_workbook(workbook_path,
    density_override_J7 = density_override_J7)
  field_out <- file.path(output_dir, "part1")
  field <- analyse_core_data(input, out = field_out, data_origin = data_origin)
  message("Field results: ", field_out, ". Review qc_flags.csv before using the maps.")
  config <- list(aoi_path = aoi_path, illustrative = FALSE, data_origin = data_origin)
  maps <- run_simple_mapping(config,
    core_path = file.path(field_out, "core_stocks.csv"),
    out = file.path(output_dir, "maps"), cellsize_m = grid_cellsize_m)
  manifest <- data.frame(input = c("Field data", "Area boundary"),
    file = c(workbook_path, aoi_path),
    md5 = unname(tools::md5sum(c(workbook_path, aoi_path))),
    data_origin = data_origin)
  utils::write.csv(manifest, file.path(output_dir, "input_files.csv"), row.names = FALSE)
  print(field$cores[c("core_id", "ecosystem", "stock_0_30_kg_m2")], row.names = FALSE)
  message("Finished. Open ", output_dir,
    "/part1 for stocks and profiles, and ", output_dir, "/maps for the maps.")
  invisible(list(field = field, maps = maps, output_dir = output_dir))
}
