suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

CARBON_INTERVALS <- tibble::tribble(
  ~interval,    ~top_cm, ~bottom_cm,
  "0_15cm",          0,         15,
  "15_30cm",        15,         30,
  "30_50cm",        30,         50,
  "50_100cm",       50,        100,
  "100_200cm",     100,        200
)

STANDARD_LAYER_COLUMNS <- c(
  "database", "database_code", "core_uid", "source_core_id", "study_id",
  "layer_id", "latitude", "longitude", "country", "region",
  "sampling_year", "ecosystem", "top_cm", "bottom_cm", "thickness_cm",
  "bulk_density_g_cm3", "carbon_pct", "carbon_basis", "stock_kg_m2",
  "stock_method", "measurement_status", "layer_quality_flag", "source_file"
)

STANDARD_METADATA_COLUMNS <- c(
  "database", "database_code", "core_uid", "source_core_id", "study_id",
  "latitude", "longitude", "country", "region", "sampling_year",
  "ecosystem", "source_files", "source_license", "source_notes"
)

assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s", label, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(data)
}

first_nonmissing <- function(x) {
  keep <- !is.na(x)
  if (is.character(x)) keep <- keep & nzchar(str_trim(x))
  if (!any(keep)) return(x[NA_integer_][1])
  x[which(keep)[1]]
}

collapse_nonmissing <- function(x, sep = ";") {
  x <- unique(str_trim(as.character(x[!is.na(x)])))
  x <- x[nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  paste(sort(x), collapse = sep)
}

safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

as_year <- function(x) {
  value <- suppressWarnings(as.integer(str_extract(as.character(x), "[12][0-9]{3}")))
  ifelse(value >= 1800L & value <= 2100L, value, NA_integer_)
}

valid_coordinate <- function(latitude, longitude) {
  !is.na(latitude) & !is.na(longitude) &
    latitude >= -90 & latitude <= 90 & longitude >= -180 & longitude <= 180
}

stock_from_percent <- function(carbon_pct, bulk_density_g_cm3, thickness_cm) {
  # kg C/m2 = C(%) * bulk density(g/cm3) * thickness(cm) / 10.
  carbon_pct * bulk_density_g_cm3 * thickness_cm / 10
}

stock_from_g_per_kg <- function(carbon_g_kg, bulk_density_g_cm3, thickness_cm) {
  # kg C/m2 = C(g/kg) * bulk density(g/cm3) * thickness(cm) / 100.
  # This is one tenth of the formula in the legacy draft WoSIS script.
  carbon_g_kg * bulk_density_g_cm3 * thickness_cm / 100
}

normalize_layers <- function(layers) {
  assert_columns(layers, STANDARD_LAYER_COLUMNS, "harmonized layer table")
  layers %>%
    mutate(
      across(all_of(c("database", "database_code", "core_uid", "source_core_id",
                      "study_id", "layer_id", "country", "region", "ecosystem",
                      "carbon_basis", "stock_method", "measurement_status",
                      "layer_quality_flag", "source_file")), as.character),
      across(all_of(c("latitude", "longitude", "top_cm", "bottom_cm",
                      "thickness_cm", "bulk_density_g_cm3", "carbon_pct",
                      "stock_kg_m2")), as.numeric),
      sampling_year = as.integer(sampling_year)
    ) %>%
    select(all_of(STANDARD_LAYER_COLUMNS))
}

normalize_metadata <- function(metadata) {
  assert_columns(metadata, STANDARD_METADATA_COLUMNS, "core metadata table")
  metadata <- metadata %>%
    mutate(
      across(all_of(c("database", "database_code", "core_uid", "source_core_id",
                      "study_id", "country", "region", "ecosystem", "source_files",
                      "source_license", "source_notes")), as.character),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude),
      sampling_year = as.integer(sampling_year)
    ) %>%
    select(all_of(STANDARD_METADATA_COLUMNS))
  if (anyDuplicated(metadata$core_uid)) {
    duplicate_ids <- unique(metadata$core_uid[duplicated(metadata$core_uid)])
    stop(
      sprintf("Core metadata contains duplicate core_uid values: %s",
              paste(head(duplicate_ids, 10), collapse = ", ")),
      call. = FALSE
    )
  }
  metadata
}

covered_length_cm <- function(top, bottom, interval_top, interval_bottom, keep = NULL) {
  if (is.null(keep)) keep <- rep(TRUE, length(top))
  good <- keep & is.finite(top) & is.finite(bottom) & bottom > top
  if (!any(good)) return(0)
  starts <- pmax(top[good], interval_top)
  ends <- pmin(bottom[good], interval_bottom)
  use <- ends > starts
  if (!any(use)) return(0)
  segments <- tibble(start = starts[use], end = ends[use]) %>% arrange(start, end)
  total <- 0
  current_start <- segments$start[1]
  current_end <- segments$end[1]
  if (nrow(segments) > 1L) {
    for (i in 2:nrow(segments)) {
      if (segments$start[i] <= current_end) {
        current_end <- max(current_end, segments$end[i])
      } else {
        total <- total + current_end - current_start
        current_start <- segments$start[i]
        current_end <- segments$end[i]
      }
    }
  }
  total + current_end - current_start
}

has_layer_overlap <- function(top, bottom, tolerance = 1e-8) {
  good <- is.finite(top) & is.finite(bottom) & bottom > top
  if (sum(good) < 2L) return(FALSE)
  segments <- tibble(top = top[good], bottom = bottom[good]) %>% arrange(top, bottom)
  prior_max <- cummax(segments$bottom)
  any(segments$top[-1] < prior_max[-length(prior_max)] - tolerance)
}

summarize_one_core <- function(core_layers, intervals = CARBON_INTERVALS) {
  valid_depth <- with(core_layers, is.finite(top_cm) & is.finite(bottom_cm) & bottom_cm > top_cm)
  overlap_flag <- has_layer_overlap(core_layers$top_cm, core_layers$bottom_cm)
  profile_top <- safe_min(core_layers$top_cm[valid_depth])
  profile_bottom <- safe_max(core_layers$bottom_cm[valid_depth])

  result <- tibble(
    core_uid = first(core_layers$core_uid),
    profile_top_cm = profile_top,
    profile_bottom_cm = profile_bottom,
    n_layers = nrow(core_layers),
    n_valid_depth_layers = sum(valid_depth),
    n_stock_layers = sum(valid_depth & is.finite(core_layers$stock_kg_m2)),
    n_invalid_depth_layers = sum(!valid_depth),
    has_overlapping_layers = overlap_flag,
    carbon_basis_summary = collapse_nonmissing(core_layers$carbon_basis),
    stock_method_summary = collapse_nonmissing(core_layers$stock_method),
    measurement_status_summary = collapse_nonmissing(core_layers$measurement_status),
    layer_quality_flags = collapse_nonmissing(core_layers$layer_quality_flag)
  )

  for (i in seq_len(nrow(intervals))) {
    label <- intervals$interval[i]
    interval_top <- intervals$top_cm[i]
    interval_bottom <- intervals$bottom_cm[i]
    width <- interval_bottom - interval_top
    overlap_cm <- pmax(
      0,
      pmin(core_layers$bottom_cm, interval_bottom) -
        pmax(core_layers$top_cm, interval_top)
    )
    overlap_cm[!valid_depth] <- 0
    stock_ok <- valid_depth & is.finite(core_layers$stock_kg_m2)
    contributions <- core_layers$stock_kg_m2 * overlap_cm / core_layers$thickness_cm
    contributions[!stock_ok | overlap_cm <= 0] <- NA_real_
    observed_stock <- if (all(is.na(contributions))) NA_real_ else sum(contributions, na.rm = TRUE)
    depth_coverage <- covered_length_cm(
      core_layers$top_cm, core_layers$bottom_cm, interval_top, interval_bottom, valid_depth
    )
    stock_coverage <- covered_length_cm(
      core_layers$top_cm, core_layers$bottom_cm, interval_top, interval_bottom, stock_ok
    )
    stock_complete <- stock_coverage >= width - 1e-6 && !overlap_flag

    result[[paste0("carbon_stock_observed_", label, "_kg_m2")]] <- observed_stock
    result[[paste0("carbon_stock_", label, "_kg_m2")]] <- if (stock_complete) observed_stock else NA_real_
    result[[paste0("depth_coverage_", label, "_fraction")]] <- min(1, depth_coverage / width)
    result[[paste0("stock_coverage_", label, "_fraction")]] <- min(1, stock_coverage / width)
    result[[paste0("carbon_complete_", label)]] <- stock_complete
  }
  result
}

build_core_table <- function(layers, metadata, intervals = CARBON_INTERVALS) {
  layers <- normalize_layers(layers)
  metadata <- normalize_metadata(metadata)
  orphan_layers <- setdiff(unique(layers$core_uid), metadata$core_uid)
  if (length(orphan_layers) > 0L) {
    stop(
      sprintf("%d layer core IDs have no metadata row; first IDs: %s",
              length(orphan_layers), paste(head(orphan_layers, 10), collapse = ", ")),
      call. = FALSE
    )
  }

  layer_summary <- layers %>%
    group_by(core_uid) %>%
    group_split(.keep = TRUE) %>%
    map_dfr(summarize_one_core, intervals = intervals)

  metadata %>%
    left_join(layer_summary, by = "core_uid") %>%
    mutate(
      coordinate_valid = valid_coordinate(latitude, longitude),
      gee_eligible = coordinate_valid,
      has_any_observed_stock = if_any(starts_with("carbon_stock_observed_"), ~ is.finite(.x)),
      has_any_complete_interval = if_any(starts_with("carbon_complete_"), ~ .x %in% TRUE)
    )
}

validate_ingest_result <- function(result, label) {
  if (!is.list(result) || !all(c("layers", "metadata") %in% names(result))) {
    stop(sprintf("%s ingester must return list(layers=..., metadata=...)", label), call. = FALSE)
  }
  result$layers <- normalize_layers(result$layers)
  result$metadata <- normalize_metadata(result$metadata)
  invisible(result)
}
