source(file.path("R", "external_profiles", "harmonization_utils.R"))
source(file.path("R", "external_profiles", "ingest_canpeat.R"))

source_root <- Sys.getenv("CARBON_RESOURCES_DIR", unset = file.path("data", "external", "Open_Core_Data"))
output_dir <- Sys.getenv(
  "EXTERNAL_CORE_OUTPUT_DIR",
  unset = file.path("data", "derived", "external_carbon_profiles")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result <- ingest_canpeat(source_root)
write_csv(result$layers, file.path(output_dir, "canpeat_layers_harmonized.csv"), na = "")
write_csv(build_core_table(result$layers, result$metadata),
          file.path(output_dir, "canpeat_cores_harmonized.csv"), na = "")
