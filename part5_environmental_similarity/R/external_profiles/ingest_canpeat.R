ingest_canpeat <- function(source_root) {
  root <- file.path(source_root, "CANPEAT_DATA")
  core_path <- file.path(root, "CANPEAT_CORES_DATA")
  profile_path <- file.path(root, "CANPEAT_PROFILES_DATA")
  # CANPEAT_CORES_DATA carries a malformed byte where a UTF-8 BOM should be
  # (0xEF with no continuation bytes), which is invalid UTF-8 and makes
  # base R's regex engine fail on that column name. Strip non-ASCII bytes
  # from the header only (data values are untouched) so `matches()` works.
  sanitize_names <- function(data) {
    names(data) <- iconv(names(data), to = "ASCII", sub = "")
    data
  }
  cores <- read_csv(core_path, show_col_types = FALSE, progress = FALSE, name_repair = "unique") %>%
    sanitize_names() %>%
    select(-any_of("...1")) %>%
    rename_with(~ "CORE_ID", matches("CORE_ID"))
  samples <- read_csv(profile_path, show_col_types = FALSE, progress = FALSE, name_repair = "unique") %>%
    sanitize_names() %>%
    select(-any_of("...1")) %>%
    rename_with(~ "CORE_ID", matches("CORE_ID")) %>%
    mutate(source_row = row_number())
  assert_columns(cores, c("CORE_ID", "LATITUDE", "LONGITUDE", "SAMPLING_YR"), "CANPEAT_CORES_DATA")
  assert_columns(
    samples,
    c("CORE_ID", "SAMPLE_NO", "UPPER_SAMP_DEPTH", "SAMP_THICK", "BULK_DENSITY",
      "C_TOT_PCT", "C_ORG_PCT", "SAMP_CARB_MGHA"),
    "CANPEAT_PROFILES_DATA"
  )

  metadata <- cores %>%
    transmute(
      database = "Canadian Peatland Database (CanPeat)",
      database_code = "CANPEAT",
      source_core_id = as.character(CORE_ID),
      core_uid = str_c("CANPEAT:", source_core_id),
      study_id = as.character(SOURCE_SITEID),
      latitude = LATITUDE,
      longitude = LONGITUDE,
      country = "Canada",
      region = as.character(PROV_TERR),
      sampling_year = as.integer(SAMPLING_YR),
      ecosystem = coalesce(as.character(CWCS_CLASS), as.character(CSSC_ORDER)),
      source_files = "CANPEAT_CORES_DATA;CANPEAT_PROFILES_DATA",
      source_license = "Not supplied in local extract; verify before redistribution",
      source_notes = "Organic carbon and bulk density include measured/estimated source flags"
    ) %>%
    select(all_of(STANDARD_METADATA_COLUMNS))

  core_slim <- metadata %>%
    select(core_uid, source_core_id, study_id, latitude, longitude, country, region,
           sampling_year, ecosystem)
  layers <- samples %>%
    mutate(source_core_id = as.character(CORE_ID)) %>%
    left_join(core_slim, by = "source_core_id") %>%
    mutate(
      carbon_pct = coalesce(C_ORG_PCT, C_TOT_PCT),
      carbon_basis = case_when(
        is.finite(C_ORG_PCT) ~ "organic_carbon_percent",
        is.finite(C_TOT_PCT) ~ "total_carbon_percent_fallback",
        is.finite(SAMP_CARB_MGHA) ~ "source_sample_carbon_stock_fallback",
        TRUE ~ NA_character_
      ),
      database = "Canadian Peatland Database (CanPeat)",
      database_code = "CANPEAT",
      layer_id = str_c("CANPEAT:", source_core_id, ":", SAMPLE_NO, ":", source_row),
      top_cm = UPPER_SAMP_DEPTH,
      bottom_cm = UPPER_SAMP_DEPTH + SAMP_THICK,
      thickness_cm = bottom_cm - top_cm,
      bulk_density_g_cm3 = BULK_DENSITY,
      stock_kg_m2 = case_when(
        is.finite(carbon_pct) & is.finite(bulk_density_g_cm3) & thickness_cm > 0 ~
          stock_from_percent(carbon_pct, bulk_density_g_cm3, thickness_cm),
        is.finite(SAMP_CARB_MGHA) & thickness_cm > 0 ~ SAMP_CARB_MGHA * 0.1,
        TRUE ~ NA_real_
      ),
      stock_method = case_when(
        is.finite(carbon_pct) & is.finite(bulk_density_g_cm3) & thickness_cm > 0 ~
          "carbon_percent_x_bulk_density_x_thickness",
        is.finite(SAMP_CARB_MGHA) ~ "source_sample_stock_Mg_ha_converted",
        TRUE ~ NA_character_
      ),
      measurement_status = pmap_chr(
        list(BD_MEAS_EST, C_ORG_MEAS_EST, C_TOT_MEAS_EST),
        ~ collapse_nonmissing(c(...))
      ),
      layer_quality_flag = case_when(
        !is.finite(carbon_pct) & !is.finite(SAMP_CARB_MGHA) ~ "carbon_missing",
        !is.finite(bulk_density_g_cm3) & !is.finite(SAMP_CARB_MGHA) ~ "bulk_density_missing",
        thickness_cm <= 0 | !is.finite(thickness_cm) ~ "invalid_depth_interval",
        TRUE ~ NA_character_
      ),
      source_file = "CANPEAT_PROFILES_DATA"
    ) %>%
    select(all_of(STANDARD_LAYER_COLUMNS))

  list(layers = layers, metadata = metadata)
}
