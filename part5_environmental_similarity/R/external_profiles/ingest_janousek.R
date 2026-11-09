ingest_janousek <- function(source_root) {
  root <- file.path(source_root, "JANOUSEK_DATA")
  location_path <- file.path(root, "Global_Core_Locations.csv")
  sample_path <- file.path(root, "Global_Core_Samples.csv")
  locations <- read_csv(
    location_path,
    show_col_types = FALSE,
    progress = FALSE,
    col_types = cols(core_id = col_character(), study_id = col_character(), .default = col_guess())
  ) %>%
    filter(core_id != "SampID")
  samples <- read_csv(
    sample_path,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA"),
    col_types = cols(
      core_id = col_character(), StudyID = col_character(), soc_percent = col_character(),
      perc_C_OM = col_character(), perc_C_C = col_character(), .default = col_guess()
    )
  ) %>%
    mutate(source_row = row_number())
  assert_columns(locations, c("study_id", "core_id", "latitude", "longitude"), "Global_Core_Locations.csv")
  assert_columns(
    samples,
    c("StudyID", "core_id", "SubSampID", "depth_min", "depth_max", "sample_length",
      "bulk_density", "soc_percent", "perc_C_C", "perc_C_OM"),
    "Global_Core_Samples.csv"
  )

  metadata <- locations %>%
    transmute(
      database = "Janousek global coastal wetland carbon cores",
      database_code = "JANOUSEK",
      source_core_id = core_id,
      core_uid = str_c("JANOUSEK:", source_core_id),
      study_id = study_id,
      latitude = suppressWarnings(as.numeric(latitude)),
      longitude = suppressWarnings(as.numeric(longitude)),
      country = NA_character_,
      region = as.character(state),
      sampling_year = NA_integer_,
      ecosystem = as.character(ecosystem),
      source_files = "Global_Core_Locations.csv;Global_Core_Samples.csv",
      source_license = "Not supplied in local extract; verify before redistribution",
      source_notes = "Location core_id is globally unique; sample StudyID mismatches are retained as QA flags"
    ) %>%
    select(all_of(STANDARD_METADATA_COLUMNS))

  core_slim <- metadata %>%
    select(core_uid, source_core_id, location_study_id = study_id, latitude, longitude,
           country, region, sampling_year, ecosystem)
  layers <- samples %>%
    mutate(source_core_id = core_id) %>%
    left_join(core_slim, by = "source_core_id") %>%
    mutate(
      soc_numeric = suppressWarnings(as.numeric(soc_percent)),
      direct_carbon_numeric = suppressWarnings(as.numeric(perc_C_C)),
      om_carbon_numeric = suppressWarnings(as.numeric(perc_C_OM)),
      carbon_pct = coalesce(soc_numeric, direct_carbon_numeric, om_carbon_numeric),
      carbon_basis = case_when(
        is.finite(soc_numeric) ~ "source_soc_percent",
        is.finite(direct_carbon_numeric) ~ "direct_carbon_percent_fallback",
        is.finite(om_carbon_numeric) ~ "organic_matter_derived_carbon_percent_fallback",
        TRUE ~ NA_character_
      ),
      database = "Janousek global coastal wetland carbon cores",
      database_code = "JANOUSEK",
      study_id = StudyID,
      layer_id = str_c("JANOUSEK:", source_core_id, ":", SubSampID, ":", source_row),
      top_cm = depth_min,
      bottom_cm = depth_max,
      thickness_cm = bottom_cm - top_cm,
      bulk_density_g_cm3 = bulk_density,
      stock_kg_m2 = if_else(
        is.finite(carbon_pct) & is.finite(bulk_density_g_cm3) & thickness_cm > 0,
        stock_from_percent(carbon_pct, bulk_density_g_cm3, thickness_cm),
        NA_real_
      ),
      stock_method = "carbon_percent_x_bulk_density_x_depth_span",
      measurement_status = "source_measurement_or_source_SOC_composite",
      layer_quality_flag = case_when(
        is.na(core_uid) ~ "location_missing",
        !is.na(location_study_id) & !is.na(StudyID) & location_study_id != StudyID &
          !str_detect(StudyID, fixed(location_study_id)) ~ "sample_location_study_id_mismatch",
        is.finite(sample_length) & is.finite(thickness_cm) & abs(sample_length - thickness_cm) > 0.11 ~
          "sample_length_differs_from_depth_span",
        !is.finite(carbon_pct) ~ "carbon_missing",
        !is.finite(bulk_density_g_cm3) ~ "bulk_density_missing",
        thickness_cm <= 0 | !is.finite(thickness_cm) ~ "invalid_depth_interval",
        TRUE ~ NA_character_
      ),
      source_file = "Global_Core_Samples.csv"
    ) %>%
    filter(!is.na(core_uid)) %>%
    select(all_of(STANDARD_LAYER_COLUMNS))

  list(layers = layers, metadata = metadata)
}
