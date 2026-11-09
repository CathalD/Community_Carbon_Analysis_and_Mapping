suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

source(file.path("R", "external_profiles", "harmonization_utils.R"))
source(file.path("R", "external_profiles", "ingest_riley.R"))
source(file.path("R", "external_profiles", "ingest_wosis.R"))
source(file.path("R", "external_profiles", "ingest_canpeat.R"))
source(file.path("R", "external_profiles", "ingest_janousek.R"))
source(file.path("R", "external_profiles", "ingest_npdb.R"))

source_root <- Sys.getenv("CARBON_RESOURCES_DIR", unset = file.path("data", "external", "Open_Core_Data"))
output_dir <- Sys.getenv(
  "EXTERNAL_CORE_OUTPUT_DIR",
  unset = file.path("data", "derived", "external_carbon_profiles")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ingesters <- list(
  RILEY = ingest_riley,
  WOSIS = ingest_wosis,
  CANPEAT = ingest_canpeat,
  JANOUSEK = ingest_janousek,
  NPDB = ingest_npdb
)

results <- imap(ingesters, function(ingester, code) {
  message("Ingesting ", code, " ...")
  result <- ingester(source_root)
  validate_ingest_result(result, code)
  core_table <- build_core_table(result$layers, result$metadata)
  write_csv(result$layers, file.path(output_dir, paste0(tolower(code), "_layers_harmonized.csv")), na = "")
  write_csv(core_table, file.path(output_dir, paste0(tolower(code), "_cores_harmonized.csv")), na = "")
  list(layers = result$layers, metadata = result$metadata, cores = core_table)
})

combined <- bind_rows(map(results, "cores")) %>%
  mutate(
    source_ecosystem = ecosystem,
    broad_ecosystem_group = case_when(
      database_code == "CANPEAT" & tolower(source_ecosystem) == "bog" ~ "peatland_bog",
      database_code == "CANPEAT" & tolower(source_ecosystem) == "fen" ~ "peatland_fen",
      database_code == "CANPEAT" & tolower(source_ecosystem) == "swamp" ~ "peatland_swamp",
      database_code == "CANPEAT" & tolower(source_ecosystem) == "marsh" ~ "peatland_marsh",
      database_code %in% c("CANPEAT", "RILEY") ~ "peatland_other_or_unspecified",
      database_code == "JANOUSEK" ~ "coastal_wetland",
      database_code == "WOSIS" & tolower(source_ecosystem) == "organic_surface" ~
        "organic_surface_profile",
      database_code == "WOSIS" ~ "soil_profile_unspecified",
      database_code == "NPDB" ~ "terrestrial_soil_profile",
      TRUE ~ "unclassified"
    )
  ) %>%
  arrange(database_code, core_uid)
if (anyDuplicated(combined$core_uid)) stop("Combined table has duplicate core_uid values", call. = FALSE)

gee_table <- combined %>%
  filter(gee_eligible) %>%
  mutate(
    # Earth Engine table import is most reliable with simple scalar fields.
    system_index = str_replace_all(core_uid, "[^A-Za-z0-9_-]", "_")
  ) %>%
  relocate(system_index, core_uid, database_code, source_core_id, latitude, longitude)

canada_gee_table <- gee_table %>%
  filter(tolower(country) == "canada") %>%
  arrange(broad_ecosystem_group, database_code, core_uid)

canada_ecosystem_table <- canada_gee_table %>%
  arrange(broad_ecosystem_group, source_ecosystem, database_code, core_uid)

canada_ecosystem_summary <- canada_ecosystem_table %>%
  group_by(broad_ecosystem_group, database_code, source_ecosystem) %>%
  summarise(
    n_cores = n(),
    n_with_observed_stock = sum(has_any_observed_stock, na.rm = TRUE),
    n_with_complete_interval = sum(has_any_complete_interval, na.rm = TRUE),
    n_with_overlapping_layers = sum(has_overlapping_layers %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(broad_ecosystem_group, database_code, desc(n_cores))

qa <- combined %>%
  group_by(database, database_code) %>%
  summarise(
    n_cores = n(),
    n_coordinate_valid = sum(coordinate_valid, na.rm = TRUE),
    n_with_observed_stock = sum(has_any_observed_stock, na.rm = TRUE),
    n_with_complete_interval = sum(has_any_complete_interval, na.rm = TRUE),
    n_with_overlapping_layers = sum(has_overlapping_layers %in% TRUE, na.rm = TRUE),
    n_complete_0_15cm = sum(carbon_complete_0_15cm %in% TRUE, na.rm = TRUE),
    n_complete_15_30cm = sum(carbon_complete_15_30cm %in% TRUE, na.rm = TRUE),
    n_complete_30_50cm = sum(carbon_complete_30_50cm %in% TRUE, na.rm = TRUE),
    n_complete_50_100cm = sum(carbon_complete_50_100cm %in% TRUE, na.rm = TRUE),
    n_complete_100_200cm = sum(carbon_complete_100_200cm %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

manifest <- tibble::tribble(
  ~database_code, ~database, ~local_inputs, ~carbon_unit, ~notes,
  "RILEY", "Riley Ontario peat dataset", "RileyData_raw.csv", "kg C/m2", "Organic C preferred; total C fallback is flagged",
  "WOSIS", "WoSIS 2023", "WOSIS_BDOD_DATA;WOSIS_ORGC_DATA;WOSIS_TOTC_DATA", "kg C/m2", "Correct g/kg conversion uses division by 100, not 10",
  "CANPEAT", "Canadian Peatland Database (CanPeat)", "CANPEAT_CORES_DATA;CANPEAT_PROFILES_DATA", "kg C/m2", "Measured/estimated flags retained",
  "JANOUSEK", "Janousek global coastal wetland carbon cores", "Global_Core_Locations.csv;Global_Core_Samples.csv", "kg C/m2", "Core ID is the location join key; source StudyID mismatches are flagged",
  "NPDB", "Agriculture and Agri-Food Canada National Pedon Database", "npdb_carbon_cores.csv;npdb_carbon_samples.csv", "kg C/m2", "Uses source-prepared horizon tables and retains source flags"
)

write_csv(combined, file.path(output_dir, "external_carbon_cores_all.csv"), na = "")
write_csv(gee_table, file.path(output_dir, "external_carbon_cores_for_gee.csv"), na = "")
write_csv(canada_gee_table, file.path(output_dir, "external_carbon_cores_canada_for_gee.csv"), na = "")
write_csv(canada_ecosystem_table, file.path(output_dir, "external_carbon_cores_canada_by_ecosystem.csv"), na = "")
write_csv(canada_ecosystem_summary, file.path(output_dir, "external_carbon_cores_canada_ecosystem_summary.csv"), na = "")
write_csv(qa, file.path(output_dir, "external_carbon_cores_qa_summary.csv"), na = "")
write_csv(manifest, file.path(output_dir, "external_carbon_source_manifest.csv"), na = "")

message("Wrote harmonized outputs to: ", normalizePath(output_dir, mustWork = TRUE))
print(qa, n = nrow(qa))
