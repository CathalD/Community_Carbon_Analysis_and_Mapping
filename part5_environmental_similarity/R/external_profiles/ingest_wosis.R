ingest_wosis <- function(source_root) {
  root <- file.path(source_root, "WOSIS_DATA")
  paths <- c(
    bd = file.path(root, "WOSIS_BDOD_DATA"),
    org = file.path(root, "WOSIS_ORGC_DATA"),
    total = file.path(root, "WOSIS_TOTC_DATA")
  )
  common <- c(
    "profile_id", "layer_id", "profile_code", "layer_name", "upper_depth",
    "lower_depth", "organic_surface", "dataset_id", "country_name", "longitude",
    "latitude", "positional_uncertainty", "region", "continent", "date", "licence"
  )
  bd <- read_csv(paths[["bd"]], show_col_types = FALSE, progress = FALSE) %>%
    select(any_of(common), bulk_density_g_cm3 = value_avg)
  org <- read_csv(paths[["org"]], show_col_types = FALSE, progress = FALSE) %>%
    select(layer_id, organic_carbon_g_kg = value_avg)
  total <- read_csv(paths[["total"]], show_col_types = FALSE, progress = FALSE) %>%
    select(layer_id, total_carbon_g_kg = value_avg)
  assert_columns(bd, c(common, "bulk_density_g_cm3"), "WOSIS_BDOD_DATA")

  layers <- bd %>%
    left_join(org, by = "layer_id") %>%
    left_join(total, by = "layer_id") %>%
    mutate(
      carbon_g_kg = coalesce(organic_carbon_g_kg, total_carbon_g_kg),
      carbon_pct = carbon_g_kg / 10,
      carbon_basis = case_when(
        is.finite(organic_carbon_g_kg) ~ "organic_carbon_g_kg",
        is.finite(total_carbon_g_kg) ~ "total_carbon_g_kg_fallback",
        TRUE ~ NA_character_
      ),
      database = "WoSIS 2023",
      database_code = "WOSIS",
      source_core_id = as.character(profile_id),
      core_uid = str_c("WOSIS:", source_core_id),
      study_id = dataset_id,
      layer_id = str_c("WOSIS:", as.character(layer_id)),
      country = country_name,
      sampling_year = as_year(date),
      ecosystem = if_else(organic_surface %in% TRUE, "organic_surface", NA_character_),
      top_cm = upper_depth,
      bottom_cm = lower_depth,
      thickness_cm = bottom_cm - top_cm,
      stock_kg_m2 = if_else(
        is.finite(carbon_g_kg) & is.finite(bulk_density_g_cm3) & thickness_cm > 0,
        stock_from_g_per_kg(carbon_g_kg, bulk_density_g_cm3, thickness_cm),
        NA_real_
      ),
      stock_method = "carbon_g_kg_x_bulk_density_x_thickness",
      measurement_status = "WoSIS standardized source value",
      layer_quality_flag = case_when(
        !is.finite(carbon_g_kg) ~ "carbon_missing",
        !is.finite(bulk_density_g_cm3) ~ "bulk_density_missing",
        thickness_cm <= 0 | !is.finite(thickness_cm) ~ "invalid_depth_interval",
        TRUE ~ NA_character_
      ),
      source_file = "WOSIS_BDOD_DATA;WOSIS_ORGC_DATA;WOSIS_TOTC_DATA"
    ) %>%
    select(all_of(STANDARD_LAYER_COLUMNS))

  metadata <- bd %>%
    mutate(
      database = "WoSIS 2023",
      database_code = "WOSIS",
      source_core_id = as.character(profile_id),
      core_uid = str_c("WOSIS:", source_core_id),
      study_id = dataset_id,
      country = country_name,
      sampling_year = as_year(date),
      ecosystem = if_else(organic_surface %in% TRUE, "organic_surface", NA_character_)
    ) %>%
    group_by(database, database_code, core_uid, source_core_id) %>%
    summarise(
      study_id = collapse_nonmissing(study_id),
      latitude = first_nonmissing(latitude),
      longitude = first_nonmissing(longitude),
      country = first_nonmissing(country),
      region = first_nonmissing(region),
      sampling_year = first_nonmissing(sampling_year),
      ecosystem = collapse_nonmissing(ecosystem),
      source_files = "WOSIS_BDOD_DATA;WOSIS_ORGC_DATA;WOSIS_TOTC_DATA",
      source_license = collapse_nonmissing(licence),
      source_notes = "Only profiles with a WoSIS bulk-density layer are retained; organic carbon is preferred over total carbon",
      .groups = "drop"
    ) %>%
    select(all_of(STANDARD_METADATA_COLUMNS))

  list(layers = layers, metadata = metadata)
}
