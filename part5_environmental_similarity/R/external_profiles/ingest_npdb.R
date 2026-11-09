ingest_npdb <- function(source_root) {
  root <- file.path(source_root, "AGCAN_DATA")
  core_path <- file.path(root, "npdb_carbon_cores.csv")
  sample_path <- file.path(root, "npdb_carbon_samples.csv")
  cores <- read_csv(
    core_path, show_col_types = FALSE, progress = FALSE,
    col_types = cols(pedon_id = col_character(), .default = col_guess())
  )
  samples <- read_csv(
    sample_path, show_col_types = FALSE, progress = FALSE,
    col_types = cols(pedon_id = col_character(), .default = col_guess())
  ) %>%
    mutate(source_row = row_number())
  assert_columns(cores, c("pedon_id", "latitude", "longitude", "cal_year"), "npdb_carbon_cores.csv")
  assert_columns(
    samples,
    c("pedon_id", "layer_id", "top_depth_cm", "bot_depth_cm", "carbon_pct",
      "bulk_density", "stock_Mgha", "stock_computable"),
    "npdb_carbon_samples.csv"
  )

  metadata <- cores %>%
    transmute(
      database = "Agriculture and Agri-Food Canada National Pedon Database",
      database_code = "NPDB",
      source_core_id = pedon_id,
      core_uid = str_c("NPDB:", source_core_id),
      study_id = NA_character_,
      latitude = latitude,
      longitude = longitude,
      country = "Canada",
      region = as.character(province),
      sampling_year = as.integer(cal_year),
      ecosystem = coalesce(as.character(drainage), NA_character_),
      source_files = "npdb_carbon_cores.csv;npdb_carbon_samples.csv",
      source_license = "Not supplied in local extract; verify before redistribution",
      source_notes = "Uses source-prepared pedon and horizon tables; carbon stock is recalculated from concentration, bulk density and thickness"
    ) %>%
    select(all_of(STANDARD_METADATA_COLUMNS))

  core_slim <- metadata %>%
    select(core_uid, source_core_id, study_id, latitude, longitude, country, region,
           sampling_year, ecosystem)
  layers <- samples %>%
    mutate(source_core_id = pedon_id) %>%
    # npdb_carbon_samples.csv repeats the core-level latitude/longitude on
    # every horizon row; drop them so the join below doesn't create
    # latitude.x/latitude.y and leave the plain columns missing.
    select(-any_of(c("latitude", "longitude"))) %>%
    left_join(core_slim, by = "source_core_id") %>%
    mutate(
      database = "Agriculture and Agri-Food Canada National Pedon Database",
      database_code = "NPDB",
      layer_id = str_c("NPDB:", source_core_id, ":", layer_id, ":", source_row),
      top_cm = top_depth_cm,
      bottom_cm = bot_depth_cm,
      thickness_cm = bottom_cm - top_cm,
      bulk_density_g_cm3 = bulk_density,
      carbon_basis = if_else(is.finite(carbon_pct), "source_carbon_percent", NA_character_),
      stock_kg_m2 = if_else(
        stock_computable %in% TRUE & is.finite(carbon_pct) &
          is.finite(bulk_density_g_cm3) & thickness_cm > 0,
        stock_from_percent(carbon_pct, bulk_density_g_cm3, thickness_cm),
        NA_real_
      ),
      stock_method = "carbon_percent_x_bulk_density_x_thickness",
      measurement_status = "source_prepared_NPDB_horizon",
      layer_quality_flag = case_when(
        flag_thick_invalid %in% TRUE | thickness_cm <= 0 | !is.finite(thickness_cm) ~ "invalid_depth_interval",
        flag_carbon_missing %in% TRUE | !is.finite(carbon_pct) ~ "carbon_missing",
        flag_bd_missing %in% TRUE | !is.finite(bulk_density_g_cm3) ~ "bulk_density_missing",
        flag_bd_low %in% TRUE ~ "bulk_density_low",
        flag_bd_high %in% TRUE ~ "bulk_density_high",
        flag_open_bottom %in% TRUE ~ "open_bottom_horizon",
        flag_negative_top %in% TRUE ~ "negative_top_depth",
        flag_carbon_high %in% TRUE ~ "carbon_high",
        TRUE ~ NA_character_
      ),
      source_file = "npdb_carbon_samples.csv"
    ) %>%
    select(all_of(STANDARD_LAYER_COLUMNS))

  list(layers = layers, metadata = metadata)
}
