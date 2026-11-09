ingest_riley <- function(source_root) {
  path <- file.path(source_root, "OGS_DATA", "RileyData_raw.csv")
  raw <- read_csv(
    path,
    skip = 1,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA", "N/A", "null")
  )
  assert_columns(
    raw,
    c("ID", "NTS", "Peatla_No.", "Sample_Point_Location", "Lat", "Long",
      "Top_Int_(CM)", "Bot_Int_(CM)", "Dry_Bulk_Dens", "C_TTL_pct", "C_Org_pct"),
    "RileyData_raw.csv"
  )

  layers <- raw %>%
    mutate(
      source_core_id = str_c(NTS, Peatla_No., Sample_Point_Location, sep = "|"),
      core_uid = str_c("RILEY", source_core_id, sep = ":"),
      carbon_pct = coalesce(C_Org_pct, C_TTL_pct),
      carbon_basis = case_when(
        is.finite(C_Org_pct) ~ "organic_carbon_percent",
        is.finite(C_TTL_pct) ~ "total_carbon_percent_fallback",
        TRUE ~ NA_character_
      ),
      database = "Riley Ontario peat dataset",
      database_code = "RILEY",
      study_id = as.character(NTS),
      layer_id = str_c("RILEY:", ID),
      latitude = Lat,
      longitude = Long,
      country = "Canada",
      region = "Ontario",
      sampling_year = NA_integer_,
      # Peat_Type_S/C/L/X are layer composition scores that sum to ten, not
      # ecosystem labels. Use the defensible broad class for grouping.
      ecosystem = "peatland",
      top_cm = `Top_Int_(CM)`,
      bottom_cm = `Bot_Int_(CM)`,
      thickness_cm = bottom_cm - top_cm,
      bulk_density_g_cm3 = Dry_Bulk_Dens,
      stock_kg_m2 = if_else(
        is.finite(carbon_pct) & is.finite(bulk_density_g_cm3) & thickness_cm > 0,
        stock_from_percent(carbon_pct, bulk_density_g_cm3, thickness_cm),
        NA_real_
      ),
      stock_method = "carbon_percent_x_dry_bulk_density_x_thickness",
      measurement_status = "source_measurement",
      layer_quality_flag = case_when(
        is.finite(C_Org_pct) & is.finite(C_TTL_pct) & C_Org_pct > C_TTL_pct + 0.1 ~
          "organic_carbon_exceeds_total_carbon",
        !is.finite(carbon_pct) ~ "carbon_missing",
        !is.finite(bulk_density_g_cm3) ~ "bulk_density_missing",
        thickness_cm <= 0 | !is.finite(thickness_cm) ~ "invalid_depth_interval",
        TRUE ~ NA_character_
      ),
      source_file = "RileyData_raw.csv"
    ) %>%
    select(all_of(STANDARD_LAYER_COLUMNS))

  metadata <- layers %>%
    group_by(database, database_code, core_uid, source_core_id, study_id) %>%
    summarise(
      latitude = first_nonmissing(latitude),
      longitude = first_nonmissing(longitude),
      country = first_nonmissing(country),
      region = first_nonmissing(region),
      sampling_year = first_nonmissing(sampling_year),
      ecosystem = collapse_nonmissing(ecosystem),
      source_files = "RileyData_raw.csv",
      source_license = "Not supplied in local extract; verify before redistribution",
      source_notes = "Core key combines NTS, peatland number and sample-point location",
      .groups = "drop"
    ) %>%
    select(all_of(STANDARD_METADATA_COLUMNS))

  list(layers = layers, metadata = metadata)
}
