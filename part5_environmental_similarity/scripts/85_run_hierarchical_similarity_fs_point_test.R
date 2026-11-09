#!/usr/bin/env Rscript

# One-command runner for the Fort Severn hierarchical-similarity point test.
# It audits the downloaded Earth Engine shards, runs the model self-test,
# predicts held-out community points, evaluates measured shallow intervals, and
# compares the result with simple equal-source donor baselines.

options(stringsAsFactors = FALSE)

DEFAULT_POINTS <- file.path("data", "external", "Hierarchical_Data")
DEFAULT_STOCKS <- file.path(
  "data", "derived", "external_carbon_profiles",
  "external_carbon_cores_all.csv"
)
DEFAULT_TRUTH <- file.path(
  "outputs", "covariate_analysis", "part7_vertical_profile",
  "community_depth_segments.csv"
)
DEFAULT_WEIGHTS <- file.path(
  "config", "hierarchical_similarity_depth_exponents.csv"
)
DEFAULT_OUTPUT <- file.path(
  "outputs", "covariate_analysis",
  "part16_hierarchical_similarity_micro_pilot", "fs_aoi_point_test_v1"
)

parse_cli <- function(args) {
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      out$help <- TRUE
    } else if (arg == "--skip-self-test") {
      out$skip_self_test <- TRUE
    } else if (startsWith(arg, "--") && grepl("=", arg, fixed = TRUE)) {
      bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
      out[[gsub("-", "_", bits[[1L]], fixed = TRUE)]] <-
        paste(bits[-1L], collapse = "=")
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  out
}

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "y", "1")
}

read_csv <- function(path) {
  utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "-9999", "-9999.0")
  )
}

write_checked <- function(x, path, overwrite) {
  if (file.exists(path) && !overwrite) {
    stop("Output exists; use --overwrite=true: ", path, call. = FALSE)
  }
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

run_rscript <- function(script, args = character()) {
  command <- file.path(R.home("bin"), "Rscript")
  status <- system2(command, c(shQuote(script), vapply(args, shQuote, character(1))))
  if (!identical(status, 0L)) {
    stop("R step failed (status ", status, "): ", script, call. = FALSE)
  }
}

audit_point_shards <- function(points_dir, output_dir, overwrite) {
  files <- sort(list.files(
    points_dir,
    pattern = "^HierarchicalSimilarity_DonorCovariates_.*[.]csv$",
    full.names = TRUE
  ))
  if (!length(files)) stop("No point shards found in: ", points_dir, call. = FALSE)
  tables <- lapply(files, read_csv)
  reference <- names(tables[[1L]])
  bad_schema <- which(!vapply(
    tables, function(x) identical(names(x), reference), logical(1)
  ))
  if (length(bad_schema)) {
    stop("Shard schema/order mismatch: ", files[[bad_schema[[1L]]]], call. = FALSE)
  }
  points <- do.call(rbind, tables)
  rownames(points) <- NULL
  required <- c("core_uid", "database_code", "record_role")
  missing <- setdiff(required, names(points))
  if (length(missing)) {
    stop("Point shards lack: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (any(!nzchar(trimws(points$core_uid))) || anyDuplicated(points$core_uid)) {
    stop("core_uid must be non-empty and unique across all shards.", call. = FALSE)
  }
  role <- tolower(trimws(points$record_role))
  if (!all(role %in% c("donor", "validation"))) {
    stop("Unexpected record_role values: ", paste(unique(role), collapse = ", "), call. = FALSE)
  }
  if (!any(role == "donor") || !any(role == "validation")) {
    stop("Both donor and validation rows are required.", call. = FALSE)
  }

  response_like <- grep(
    "(^carbon_|stock|coverage|complete_interval|profile_|target_kg)",
    names(points), value = TRUE, ignore.case = TRUE
  )
  if (length(response_like)) {
    stop(
      "Response-like fields leaked into Earth Engine output: ",
      paste(response_like, collapse = ", "), call. = FALSE
    )
  }

  block_names <- list(
    local = grep("^local_", names(points), value = TRUE),
    regional = grep("^regional_", names(points), value = TRUE),
    broad = setdiff(
      grep("^broad_", names(points), value = TRUE),
      "broad_ecosystem_group"
    ),
    qa = grep("^qa_", names(points), value = TRUE)
  )
  expected_counts <- c(local = 31L, regional = 7L, broad = 5L, qa = 4L)
  observed_counts <- lengths(block_names)
  if (!identical(observed_counts[names(expected_counts)], expected_counts)) {
    stop(
      "Unexpected v1 band counts. Expected ",
      paste(names(expected_counts), expected_counts, collapse = "; "),
      "; found ", paste(names(observed_counts), observed_counts, collapse = "; "),
      call. = FALSE
    )
  }
  bands <- unlist(block_names, use.names = FALSE)
  for (name in bands) points[[name]] <- suppressWarnings(as.numeric(points[[name]]))

  file_audit <- data.frame(
    file = basename(files),
    rows = vapply(tables, nrow, integer(1)),
    columns = vapply(tables, ncol, integer(1)),
    bytes = file.info(files)$size,
    stringsAsFactors = FALSE
  )
  role_database <- as.data.frame(
    table(record_role = role, database_code = points$database_code),
    stringsAsFactors = FALSE
  )
  role_database <- role_database[role_database$Freq > 0, , drop = FALSE]
  missingness <- data.frame(
    covariate = bands,
    scale_block = sub("_.*$", "", bands),
    missing_n = vapply(points[bands], function(x) sum(!is.finite(x)), integer(1)),
    missing_fraction = vapply(points[bands], function(x) mean(!is.finite(x)), numeric(1)),
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    diagnostic = c(
      "csv_shards", "total_points", "donor_points", "validation_points",
      "duplicate_core_uid", "schemas_identical", "response_like_fields",
      "local_bands", "regional_bands", "broad_bands", "qa_bands"
    ),
    value = c(
      length(files), nrow(points), sum(role == "donor"),
      sum(role == "validation"), sum(duplicated(points$core_uid)), TRUE,
      length(response_like), observed_counts[c("local", "regional", "broad", "qa")]
    ),
    stringsAsFactors = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_checked(file_audit, file.path(output_dir, "input_file_audit.csv"), overwrite)
  write_checked(role_database, file.path(output_dir, "input_role_database_counts.csv"), overwrite)
  write_checked(missingness, file.path(output_dir, "input_covariate_missingness.csv"), overwrite)
  write_checked(summary, file.path(output_dir, "input_audit_summary.csv"), overwrite)
  list(points = points, files = files, summary = summary)
}

make_baseline_comparison <- function(points, stocks_path, validation_path) {
  stocks <- read_csv(stocks_path)
  validation <- read_csv(validation_path)
  donor_ids <- points$core_uid[tolower(points$record_role) == "donor"]
  stocks <- stocks[stocks$core_uid %in% donor_ids, , drop = FALSE]
  intervals <- intersect(
    unique(validation$interval_code), c("0_15cm", "15_30cm")
  )
  groups <- c("ALL_COMMUNITY", "FS_ONLY")
  rows <- list()
  for (interval in intervals) {
    stock_name <- paste0("carbon_stock_", interval, "_kg_m2")
    complete_name <- paste0("carbon_complete_", interval)
    eligible <- stocks[
      truthy(stocks[[complete_name]]) & is.finite(stocks[[stock_name]]),
      , drop = FALSE
    ]
    source_mean <- tapply(
      eligible[[stock_name]], eligible$database_code, mean, na.rm = TRUE
    )
    source_median <- tapply(
      eligible[[stock_name]], eligible$database_code, median, na.rm = TRUE
    )
    sources <- intersect(c("CANPEAT", "RILEY"), names(source_mean))
    baseline_mean <- mean(source_mean[sources])
    baseline_median <- mean(source_median[sources])
    for (group in groups) {
      observed <- validation[validation$interval_code == interval, , drop = FALSE]
      if (group == "FS_ONLY") {
        observed <- observed[grepl("^FS", observed$plot_id, ignore.case = TRUE), ]
      }
      model_error <- observed$prediction_kg_m2 - observed$observed_kg_m2
      mean_error <- baseline_mean - observed$observed_kg_m2
      median_error <- baseline_median - observed$observed_kg_m2
      rows[[length(rows) + 1L]] <- data.frame(
        validation_group = group,
        interval_code = interval,
        n = nrow(observed),
        model_rmse_kg_m2 = sqrt(mean(model_error^2)),
        equal_source_mean_baseline_kg_m2 = baseline_mean,
        mean_baseline_rmse_kg_m2 = sqrt(mean(mean_error^2)),
        rmse_improvement_vs_mean_fraction =
          1 - sqrt(mean(model_error^2)) / sqrt(mean(mean_error^2)),
        equal_source_median_baseline_kg_m2 = baseline_median,
        median_baseline_rmse_kg_m2 = sqrt(mean(median_error^2)),
        rmse_improvement_vs_median_fraction =
          1 - sqrt(mean(model_error^2)) / sqrt(mean(median_error^2)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

main <- function() {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (isTRUE(cli$help)) {
    cat(paste0(
      "Fort Severn hierarchical-similarity point test\n\n",
      "Run from the project root:\n",
      "  Rscript scripts/85_run_hierarchical_similarity_fs_point_test.R ",
      "--overwrite=true\n\n",
      "Optional: --points=PATH --stocks=PATH --truth=PATH --weights=PATH ",
      "--output=PATH --skip-self-test\n"
    ))
    return(invisible(TRUE))
  }
  points_dir <- if (!is.null(cli$points)) cli$points else DEFAULT_POINTS
  stocks_path <- if (!is.null(cli$stocks)) cli$stocks else DEFAULT_STOCKS
  truth_path <- if (!is.null(cli$truth)) cli$truth else DEFAULT_TRUTH
  weights_path <- if (!is.null(cli$weights)) cli$weights else DEFAULT_WEIGHTS
  output_dir <- if (!is.null(cli$output)) cli$output else DEFAULT_OUTPUT
  overwrite <- !is.null(cli$overwrite) && truthy(cli$overwrite)

  for (path in c(points_dir, stocks_path, truth_path, weights_path)) {
    if (!file.exists(path)) stop("Required input does not exist: ", path, call. = FALSE)
  }
  audit <- audit_point_shards(points_dir, output_dir, overwrite)
  cat("Input audit passed: ", nrow(audit$points), " points in ",
      length(audit$files), " shards.\n", sep = "")

  model_script <- file.path("scripts", "83_run_hierarchical_similarity_micro_pilot.R")
  evaluation_script <- file.path("scripts", "84_evaluate_hierarchical_similarity_community.R")
  if (!isTRUE(cli$skip_self_test)) run_rscript(model_script, "--self-test")
  run_rscript(model_script, c(
    paste0("--points=", points_dir),
    paste0("--stocks=", stocks_path),
    paste0("--weights=", weights_path),
    paste0("--output=", output_dir),
    paste0("--overwrite=", tolower(overwrite))
  ))
  predictions <- file.path(output_dir, "point_predictions.csv")
  run_rscript(evaluation_script, c(
    paste0("--predictions=", predictions),
    paste0("--truth=", truth_path),
    paste0("--output=", output_dir),
    paste0("--overwrite=", tolower(overwrite))
  ))

  validation_path <- file.path(output_dir, "community_validation_rows.csv")
  baseline <- make_baseline_comparison(
    audit$points, stocks_path, validation_path
  )
  baseline_path <- file.path(
    output_dir, "community_validation_baseline_comparison.csv"
  )
  write_checked(baseline, baseline_path, overwrite)
  cat("\nSimilarity model versus simple donor baselines:\n")
  print(baseline, row.names = FALSE)
  cat("\nComplete run output: ", normalizePath(output_dir), "\n", sep = "")
  invisible(baseline)
}

if (sys.nframe() == 0L) main()
