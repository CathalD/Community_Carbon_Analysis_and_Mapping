#!/usr/bin/env Rscript

# Hierarchical, multi-scale similarity micro-pilot.
#
# This is an isolated research prototype. It reads covariates exported at points,
# joins complete five-depth donor stocks, and predicts each target independently.
# It never reads from or writes to outputs/route_l.

options(stringsAsFactors = FALSE)

DEPTH_SPEC <- data.frame(
  interval_code = c("0_15cm", "15_30cm", "30_50cm", "50_100cm", "100_200cm"),
  top_cm = c(0, 15, 30, 50, 100),
  bottom_cm = c(15, 30, 50, 100, 200),
  local = c(0.50, 0.40, 0.30, 0.20, 0.10),
  regional = c(0.30, 0.35, 0.40, 0.40, 0.35),
  broad = c(0.20, 0.25, 0.30, 0.40, 0.55),
  stringsAsFactors = FALSE
)

DEFAULT_OUTPUT_DIR <- file.path(
  "outputs", "covariate_analysis", "part16_hierarchical_similarity_micro_pilot"
)

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(is.finite(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "y", "1", "complete")
}

parse_boolean <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  value <- truthy(x)
  if (length(value) != 1L) stop("Expected one boolean value.", call. = FALSE)
  value
}

first_nonempty <- function(...) {
  values <- list(...)
  for (value in values) {
    if (!is.null(value) && length(value) && !is.na(value[[1L]]) && nzchar(value[[1L]])) {
      return(value[[1L]])
    }
  }
  NULL
}

parse_cli <- function(args) {
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      out$help <- TRUE
      next
    }
    if (arg == "--self-test") {
      out$self_test <- TRUE
      next
    }
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Unknown argument: ", arg, ". Use --help for the interface.", call. = FALSE)
    }
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", bits[[1L]], fixed = TRUE)
    out[[key]] <- paste(bits[-1L], collapse = "=")
  }
  out
}

usage <- function() {
  cat(paste0(
    "Hierarchical multi-scale similarity micro-pilot\n\n",
    "Required:\n",
    "  --points=PATH       One point CSV or a directory of CSV shards. Requires\n",
    "                      identical columns, unique core_uid, and numeric\n",
    "                      local_*, regional_* and broad_* fields.\n",
    "  --stocks=PATH       Harmonized stocks with core_uid, source, five stock\n",
    "                      fields, and explicit completeness flags.\n\n",
    "Optional:\n",
    "  --output=PATH       Default: ", DEFAULT_OUTPUT_DIR, "\n",
    "  --weights=PATH      CSV overriding depth exponents. Columns:\n",
    "                      interval_code, local, regional, broad.\n",
    "  --raster=PATH       Small multi-band GeoTIFF with the same covariates.\n",
    "  --id-column=NAME    Default: core_uid.\n",
    "  --role-column=NAME  Default: record_role when present.\n",
    "  --target-roles=LIST Default: target,community,validation,prediction.\n",
    "  --top-n=N           Top donors retained per source/depth; default 5.\n",
    "  --min-source-donors=N  Default 3.\n",
    "  --min-block-overlap-fraction=X  Pairwise variables required per block; default 0.50.\n",
    "  --support-warning-threshold=X  Default 0.10.\n",
    "  --allow-inferred-completeness=true  Opt out of strict completeness flags.\n",
    "  --overwrite=true    Replace earlier Part 16 files only; default false.\n",
    "  --self-test         Run a dependency-free synthetic test.\n\n",
    "The same options may be supplied with HMS_POINTS_CSV, HMS_STOCKS_CSV,\n",
    "HMS_OUTPUT_DIR, HMS_WEIGHTS_CSV and HMS_RASTER environment variables.\n"
  ))
}

read_csv <- function(path) {
  if (!file.exists(path)) stop("Input does not exist: ", path, call. = FALSE)
  utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "-9999", "-9999.0")
  )
}

read_point_csvs <- function(path) {
  if (!dir.exists(path)) {
    out <- read_csv(path)
    attr(out, "point_shard_count") <- 1L
    attr(out, "point_shard_manifest") <- data.frame(
      file = normalizePath(path, mustWork = TRUE), rows = nrow(out), stringsAsFactors = FALSE
    )
    return(out)
  }
  files <- sort(list.files(
    path, pattern = "^HierarchicalSimilarity_DonorCovariates_.*\\.csv$",
    full.names = TRUE, recursive = TRUE
  ))
  if (!length(files)) {
    stop(
      "Point directory contains no donor shards named ",
      "HierarchicalSimilarity_DonorCovariates_*.csv: ", path, call. = FALSE
    )
  }
  tables <- lapply(files, read_csv)
  reference_names <- names(tables[[1L]])
  mismatch <- which(!vapply(tables, function(x) identical(names(x), reference_names), logical(1)))
  if (length(mismatch)) {
    stop(
      "Point CSV shards do not have identical columns/order. First mismatch: ",
      files[[mismatch[[1L]]]], call. = FALSE
    )
  }
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  attr(out, "point_shard_count") <- length(files)
  attr(out, "point_shard_manifest") <- data.frame(
    file = normalizePath(files, mustWork = TRUE),
    rows = vapply(tables, nrow, integer(1)), stringsAsFactors = FALSE
  )
  out
}

canonical_interval <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("\\s", "", z)
  z <- gsub("-", "_", z, fixed = TRUE)
  z <- gsub("to", "_", z, fixed = TRUE)
  z <- gsub("cm", "", z, fixed = TRUE)
  z <- gsub("__+", "_", z)
  out <- paste0(z, "cm")
  out[out == "0_15cm"] <- "0_15cm"
  out[out == "15_30cm"] <- "15_30cm"
  out[out == "30_50cm"] <- "30_50cm"
  out[out == "50_100cm"] <- "50_100cm"
  out[out == "100_200cm"] <- "100_200cm"
  out
}

choose_column <- function(data, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit)) return(hit[[1L]])
  if (required) {
    stop(label, " column not found. Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  NULL
}

load_depth_weights <- function(path = NULL) {
  weights <- DEPTH_SPEC
  if (!is.null(path)) {
    custom <- read_csv(path)
    required <- c("interval_code", "local", "regional", "broad")
    if (!all(required %in% names(custom))) {
      stop("Weight CSV must contain: ", paste(required, collapse = ", "), call. = FALSE)
    }
    custom$interval_code <- canonical_interval(custom$interval_code)
    if (anyDuplicated(custom$interval_code)) stop("Repeated interval in weight CSV.", call. = FALSE)
    index <- match(weights$interval_code, custom$interval_code)
    if (anyNA(index)) stop("Weight CSV must contain all five depth intervals.", call. = FALSE)
    weights[, c("local", "regional", "broad")] <- custom[index, c("local", "regional", "broad")]
  }
  for (name in c("local", "regional", "broad")) {
    weights[[name]] <- suppressWarnings(as.numeric(weights[[name]]))
  }
  exponent_sum <- rowSums(weights[, c("local", "regional", "broad")])
  if (any(!is.finite(as.matrix(weights[, c("local", "regional", "broad")]))) ||
      any(as.matrix(weights[, c("local", "regional", "broad")]) < 0) ||
      any(!is.finite(exponent_sum)) || any(exponent_sum <= 0)) {
    stop("Depth exponents must be finite, non-negative, and sum to more than zero.", call. = FALSE)
  }
  weights[, c("local", "regional", "broad")] <-
    weights[, c("local", "regional", "broad")] / exponent_sum
  weights$exponent_sum <- rowSums(weights[, c("local", "regional", "broad")])
  weights
}

prepare_stock_long <- function(stocks, id_column, allow_inferred_completeness = FALSE) {
  if (!id_column %in% names(stocks)) stop("Stocks lack ID column: ", id_column, call. = FALSE)
  source_column <- choose_column(
    stocks, c("database_code", "donor_source", "source", "database"), "Donor source"
  )

  if ("interval_code" %in% names(stocks)) {
    response_column <- choose_column(
      stocks,
      c("target_kg_m2", "carbon_stock_kg_m2", "stock_kg_m2", "observed_kg_m2"),
      "Long-table stock"
    )
    complete_column <- choose_column(
      stocks, c("complete_interval", "carbon_complete", "stock_complete"),
      "Long-table completeness", required = !allow_inferred_completeness
    )
    complete <- if (is.null(complete_column)) {
      is.finite(suppressWarnings(as.numeric(stocks[[response_column]])))
    } else {
      truthy(stocks[[complete_column]])
    }
    out <- data.frame(
      core_uid = as.character(stocks[[id_column]]),
      donor_source = as.character(stocks[[source_column]]),
      interval_code = canonical_interval(stocks$interval_code),
      stock_kg_m2 = suppressWarnings(as.numeric(stocks[[response_column]])),
      complete_interval = complete,
      stringsAsFactors = FALSE
    )
  } else {
    rows <- vector("list", nrow(DEPTH_SPEC))
    for (i in seq_len(nrow(DEPTH_SPEC))) {
      code <- DEPTH_SPEC$interval_code[[i]]
      stock_column <- choose_column(
        stocks,
        c(
          paste0("carbon_stock_", code, "_kg_m2"),
          paste0("stock_", code, "_kg_m2"),
          paste0("target_", code, "_kg_m2"),
          paste0(code, "_kg_m2")
        ),
        paste0("Stock for ", code)
      )
      complete_column <- choose_column(
        stocks,
        c(
          paste0("carbon_complete_", code),
          paste0("complete_", code),
          paste0("stock_complete_", code)
        ),
        paste0("Completeness flag for ", code),
        required = !allow_inferred_completeness
      )
      values <- suppressWarnings(as.numeric(stocks[[stock_column]]))
      complete <- if (is.null(complete_column)) is.finite(values) else truthy(stocks[[complete_column]])
      rows[[i]] <- data.frame(
        core_uid = as.character(stocks[[id_column]]),
        donor_source = as.character(stocks[[source_column]]),
        interval_code = code,
        stock_kg_m2 = values,
        complete_interval = complete,
        stringsAsFactors = FALSE
      )
    }
    out <- do.call(rbind, rows)
  }

  out <- out[out$interval_code %in% DEPTH_SPEC$interval_code, , drop = FALSE]
  out$core_uid <- trimws(out$core_uid)
  out$donor_source <- trimws(out$donor_source)
  if (any(!nzchar(out$core_uid)) || any(!nzchar(out$donor_source))) {
    stop("Every donor stock row must have a non-empty ID and source.", call. = FALSE)
  }
  if (anyDuplicated(out[, c("core_uid", "interval_code")])) {
    stop("Stocks repeat a core/depth combination.", call. = FALSE)
  }
  source_count <- vapply(split(out$donor_source, out$core_uid), function(z) length(unique(z)), integer(1))
  if (any(source_count != 1L)) stop("A donor core is assigned to more than one source.", call. = FALSE)
  out$complete_interval <- out$complete_interval & is.finite(out$stock_kg_m2) & out$stock_kg_m2 >= 0
  out
}

prepare_points <- function(
    points, id_column, role_column = NULL, target_roles = NULL,
    allow_no_targets = FALSE) {
  if (!id_column %in% names(points)) stop("Points lack ID column: ", id_column, call. = FALSE)
  points[[id_column]] <- trimws(as.character(points[[id_column]]))
  if (any(!nzchar(points[[id_column]])) || anyDuplicated(points[[id_column]])) {
    stop("Point IDs must be non-empty and unique.", call. = FALSE)
  }
  blocks <- list(
    local = grep("^local_", names(points), value = TRUE),
    regional = grep("^regional_", names(points), value = TRUE),
    broad = grep("^broad_", names(points), value = TRUE)
  )
  # Prefixes define the scale blocks, but metadata can legitimately share a
  # prefix (notably broad_ecosystem_group). Only numeric or numeric-like fields
  # may enter continuous distance. This preserves optional numeric predictors
  # while preventing character metadata from appearing as an all-missing band.
  numeric_like <- function(x) {
    if (is.numeric(x) || is.integer(x) || is.logical(x)) return(TRUE)
    text <- trimws(as.character(x))
    present <- !is.na(text) & nzchar(text)
    if (!any(present)) return(TRUE)
    parsed <- suppressWarnings(as.numeric(text[present]))
    all(is.finite(parsed))
  }
  blocks <- lapply(blocks, function(candidates) {
    candidates[vapply(points[candidates], numeric_like, logical(1))]
  })
  if (any(lengths(blocks) == 0L)) {
    missing <- names(blocks)[lengths(blocks) == 0L]
    stop("No covariates found for block(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  covariates <- unlist(blocks, use.names = FALSE)
  categorical_codes <- covariates[grepl("_categorical_code$", covariates, ignore.case = TRUE)]
  if (length(categorical_codes)) {
    stop(
      "Raw categorical codes cannot enter continuous similarity distance: ",
      paste(categorical_codes, collapse = ", "),
      ". One-hot encode them or create an explicit exact-match similarity before running.",
      call. = FALSE
    )
  }
  for (name in covariates) points[[name]] <- suppressWarnings(as.numeric(points[[name]]))

  if (is.null(role_column) && "record_role" %in% names(points)) role_column <- "record_role"
  if (!is.null(role_column)) {
    if (!role_column %in% names(points)) stop("Role column not found: ", role_column, call. = FALSE)
    role <- tolower(trimws(as.character(points[[role_column]])))
    target <- role %in% tolower(target_roles)
    if (!any(target) && !allow_no_targets) {
      stop("Role column exists, but no rows match --target-roles.", call. = FALSE)
    }
  } else {
    target <- rep(TRUE, nrow(points))
  }
  list(data = points, blocks = blocks, covariates = covariates, target = target)
}

fit_robust_scaling <- function(donor_points, blocks) {
  rows <- list()
  dropped <- character()
  for (block in names(blocks)) {
    for (name in blocks[[block]]) {
      values <- donor_points[[name]]
      centre <- suppressWarnings(stats::median(values, na.rm = TRUE))
      mad_spread <- suppressWarnings(stats::mad(values, center = centre, constant = 1.4826, na.rm = TRUE))
      sd_spread <- suppressWarnings(stats::sd(values, na.rm = TRUE))
      mad_usable <- is.finite(mad_spread) && mad_spread > 1e-12
      sd_usable <- is.finite(sd_spread) && sd_spread > 1e-12
      spread <- if (mad_usable) mad_spread else if (sd_usable) sd_spread else NA_real_
      scale_method <- if (mad_usable) "MAD" else if (sd_usable) "SD_fallback" else "excluded"
      usable <- is.finite(centre) && is.finite(spread) && spread > 1e-12
      rows[[length(rows) + 1L]] <- data.frame(
        covariate = name, scale_block = block, donor_median = centre,
        donor_mad = mad_spread, donor_sd = sd_spread, donor_scale = spread,
        scale_method = scale_method,
        donor_nonmissing_n = sum(is.finite(values)),
        donor_missing_fraction = mean(!is.finite(values)),
        included = usable,
        exclusion_reason = if (usable) "" else "non_finite_or_zero_donor_MAD_and_SD",
        stringsAsFactors = FALSE
      )
      if (!usable) dropped <- c(dropped, name)
    }
  }
  table <- do.call(rbind, rows)
  kept_blocks <- lapply(blocks, setdiff, y = dropped)
  if (any(lengths(kept_blocks) == 0L)) {
    stop(
      "Robust scaling removed every covariate in block(s): ",
      paste(names(kept_blocks)[lengths(kept_blocks) == 0L], collapse = ", "),
      call. = FALSE
    )
  }
  list(table = table, blocks = kept_blocks)
}

scale_matrix <- function(data, variables, scaling_table) {
  values <- as.matrix(data[, variables, drop = FALSE])
  storage.mode(values) <- "double"
  index <- match(variables, scaling_table$covariate)
  values <- sweep(values, 2, scaling_table$donor_median[index], "-")
  sweep(values, 2, scaling_table$donor_scale[index], "/")
}

weighted_quantile <- function(values, weights, probabilities = c(0.1, 0.5, 0.9)) {
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(keep)) return(rep(NA_real_, length(probabilities)))
  values <- values[keep]
  weights <- weights[keep]
  order_index <- order(values)
  values <- values[order_index]
  weights <- weights[order_index] / sum(weights)
  cumulative <- cumsum(weights)
  vapply(probabilities, function(p) values[which(cumulative >= p)[1L]], numeric(1))
}

prepare_model <- function(points_info, stock_long, id_column, min_block_overlap_fraction = 0.50) {
  points <- points_info$data
  initial_complete_stock_ids <- unique(stock_long$core_uid[stock_long$complete_interval])
  unmatched_stock_ids <- setdiff(initial_complete_stock_ids, points[[id_column]])
  donor_ids <- initial_complete_stock_ids
  donor_index <- match(donor_ids, points[[id_column]])
  matched <- !is.na(donor_index)
  donor_ids <- donor_ids[matched]
  donor_index <- donor_index[matched]
  if (!length(donor_ids)) stop("No complete stock donors match point covariates.", call. = FALSE)

  donor_points <- points[donor_index, , drop = FALSE]
  scaling <- fit_robust_scaling(donor_points, points_info$blocks)
  blocks <- scaling$blocks
  covariates <- unlist(blocks, use.names = FALSE)
  if (nrow(donor_points) < 2L) {
    stop("Fewer than two donors have matched point covariates.", call. = FALSE)
  }
  minimum_overlap <- pmax(1L, ceiling(lengths(blocks) * min_block_overlap_fraction))
  names(minimum_overlap) <- names(blocks)

  donor_source <- vapply(donor_ids, function(id) {
    unique(stock_long$donor_source[stock_long$core_uid == id])[[1L]]
  }, character(1))
  donor_scaled <- lapply(blocks, function(vars) scale_matrix(donor_points, vars, scaling$table))
  stock_long <- stock_long[stock_long$core_uid %in% donor_ids, , drop = FALSE]

  list(
    donor_ids = donor_ids,
    donor_source = donor_source,
    donor_scaled = donor_scaled,
    donor_points = donor_points,
    stock_long = stock_long,
    scaling = scaling$table,
    blocks = blocks,
    minimum_overlap = minimum_overlap,
    min_block_overlap_fraction = min_block_overlap_fraction,
    covariates = covariates,
    initial_complete_stock_ids = initial_complete_stock_ids,
    unmatched_stock_ids = unmatched_stock_ids
  )
}

empty_prediction_row <- function(target_id, depth, status, missing_covariates = "") {
  data.frame(
    target_id = target_id,
    interval_code = depth$interval_code,
    top_cm = depth$top_cm,
    bottom_cm = depth$bottom_cm,
    prediction_kg_m2 = NA_real_,
    mixture_p10_kg_m2 = NA_real_,
    mixture_median_kg_m2 = NA_real_,
    mixture_p90_kg_m2 = NA_real_,
    predictive_spread_sd_kg_m2 = NA_real_,
    within_source_sd_kg_m2 = NA_real_,
    between_source_sd_kg_m2 = NA_real_,
    n_sources = 0L,
    n_donors = 0L,
    mean_effective_donors_per_source = NA_real_,
    local_support = NA_real_,
    regional_support = NA_real_,
    broad_support = NA_real_,
    local_overlap_fraction = NA_real_,
    regional_overlap_fraction = NA_real_,
    broad_overlap_fraction = NA_real_,
    combined_support = NA_real_,
    maximum_combined_similarity = NA_real_,
    local_exponent = depth$local,
    regional_exponent = depth$regional,
    broad_exponent = depth$broad,
    status = status,
    missing_covariates = missing_covariates,
    stringsAsFactors = FALSE
  )
}

predict_one_target <- function(
    target_row, target_id, model, depth_weights, min_source_donors = 3L,
    support_warning_threshold = 0.10, top_n = 5L, keep_details = TRUE) {
  target_finite <- vapply(
    target_row[model$covariates],
    function(x) is.finite(suppressWarnings(as.numeric(x[[1L]]))),
    logical(1)
  )
  missing <- model$covariates[!target_finite]
  target_block_n <- vapply(
    model$blocks, function(vars) sum(target_finite[vars]), integer(1)
  )
  insufficient_blocks <- names(target_block_n)[
    target_block_n < model$minimum_overlap[names(target_block_n)]
  ]
  if (length(insufficient_blocks)) {
    rows <- do.call(rbind, lapply(seq_len(nrow(depth_weights)), function(i) {
      empty_prediction_row(
        target_id, depth_weights[i, , drop = FALSE],
        paste0("INSUFFICIENT_TARGET_BLOCK_OVERLAP:", paste(insufficient_blocks, collapse = ";")),
        paste(missing, collapse = ";")
      )
    }))
    return(list(predictions = rows, sources = NULL, donors = NULL, blocks = NULL))
  }

  block_distance <- list()
  block_similarity <- list()
  block_overlap_n <- list()
  block_overlap_fraction <- list()
  for (block in names(model$blocks)) {
    variables <- model$blocks[[block]]
    target_scaled <- scale_matrix(target_row, variables, model$scaling)
    donor_values <- model$donor_scaled[[block]]
    target_values <- as.numeric(target_scaled)
    pairwise_valid <- is.finite(donor_values) & matrix(
      is.finite(target_values), nrow(donor_values), ncol(donor_values), byrow = TRUE
    )
    overlap_n <- rowSums(pairwise_valid)
    difference <- sweep(donor_values, 2, target_values, "-")
    squared_difference <- difference^2
    squared_difference[!pairwise_valid] <- 0
    distance <- rowSums(squared_difference) / pmax(overlap_n, 1L)
    distance[overlap_n < model$minimum_overlap[[block]]] <- NA_real_
    block_distance[[block]] <- distance
    block_similarity[[block]] <- exp(-0.5 * distance)
    block_overlap_n[[block]] <- overlap_n
    block_overlap_fraction[[block]] <- overlap_n / length(variables)
  }

  prediction_rows <- list()
  source_rows <- list()
  donor_rows <- list()
  support_rows <- list()
  for (d in seq_len(nrow(depth_weights))) {
    depth <- depth_weights[d, , drop = FALSE]
    stock <- model$stock_long[
      model$stock_long$interval_code == depth$interval_code & model$stock_long$complete_interval,
      , drop = FALSE
    ]
    donor_position <- match(stock$core_uid, model$donor_ids)
    keep <- !is.na(donor_position) & stock$core_uid != target_id
    stock <- stock[keep, , drop = FALSE]
    donor_position <- donor_position[keep]
    if (!nrow(stock)) {
      prediction_rows[[d]] <- empty_prediction_row(target_id, depth, "NO_COMPLETE_DONORS")
      next
    }

    local_similarity <- block_similarity$local[donor_position]
    regional_similarity <- block_similarity$regional[donor_position]
    broad_similarity <- block_similarity$broad[donor_position]
    local_distance <- block_distance$local[donor_position]
    regional_distance <- block_distance$regional[donor_position]
    broad_distance <- block_distance$broad[donor_position]
    local_overlap_n <- block_overlap_n$local[donor_position]
    regional_overlap_n <- block_overlap_n$regional[donor_position]
    broad_overlap_n <- block_overlap_n$broad[donor_position]
    local_overlap <- block_overlap_fraction$local[donor_position]
    regional_overlap <- block_overlap_fraction$regional[donor_position]
    broad_overlap <- block_overlap_fraction$broad[donor_position]
    valid_similarity <- is.finite(local_similarity) & is.finite(regional_similarity) &
      is.finite(broad_similarity)
    stock <- stock[valid_similarity, , drop = FALSE]
    donor_position <- donor_position[valid_similarity]
    local_similarity <- local_similarity[valid_similarity]
    regional_similarity <- regional_similarity[valid_similarity]
    broad_similarity <- broad_similarity[valid_similarity]
    local_distance <- local_distance[valid_similarity]
    regional_distance <- regional_distance[valid_similarity]
    broad_distance <- broad_distance[valid_similarity]
    local_overlap_n <- local_overlap_n[valid_similarity]
    regional_overlap_n <- regional_overlap_n[valid_similarity]
    broad_overlap_n <- broad_overlap_n[valid_similarity]
    local_overlap <- local_overlap[valid_similarity]
    regional_overlap <- regional_overlap[valid_similarity]
    broad_overlap <- broad_overlap[valid_similarity]
    if (!nrow(stock)) {
      prediction_rows[[d]] <- empty_prediction_row(target_id, depth, "NO_DONORS_MEET_BLOCK_OVERLAP")
      next
    }
    log_weight <- depth$local * log(pmax(local_similarity, .Machine$double.xmin)) +
      depth$regional * log(pmax(regional_similarity, .Machine$double.xmin)) +
      depth$broad * log(pmax(broad_similarity, .Machine$double.xmin))
    combined_similarity <- exp(log_weight)
    positive_weight <- which(is.finite(combined_similarity) & combined_similarity > 0)
    source_groups <- split(positive_weight, stock$donor_source[positive_weight])
    source_groups <- source_groups[lengths(source_groups) >= min_source_donors]
    if (!length(source_groups)) {
      prediction_rows[[d]] <- empty_prediction_row(target_id, depth, "NO_SOURCE_MEETS_MINIMUM_DONORS")
      next
    }

    this_sources <- list()
    this_donors <- list()
    this_support <- list()
    mixture_values <- numeric()
    mixture_weights <- numeric()
    source_means <- numeric()
    source_variances <- numeric()
    source_effective_n <- numeric()
    source_support <- numeric()
    source_maximum <- numeric()
    source_block_weighted <- list(local = numeric(), regional = numeric(), broad = numeric())
    source_overlap_weighted <- list(local = numeric(), regional = numeric(), broad = numeric())

    for (source in names(source_groups)) {
      group <- source_groups[[source]]
      raw_weight <- combined_similarity[group]
      if (!any(is.finite(raw_weight) & raw_weight > 0)) next
      normalized <- raw_weight / sum(raw_weight)
      response <- stock$stock_kg_m2[group]
      estimate <- sum(normalized * response)
      variance <- sum(normalized * (response - estimate)^2)
      effective_n <- 1 / sum(normalized^2)
      quantiles <- weighted_quantile(response, normalized)
      weighted_block <- c(
        local = sum(normalized * local_similarity[group]),
        regional = sum(normalized * regional_similarity[group]),
        broad = sum(normalized * broad_similarity[group])
      )
      weighted_overlap <- c(
        local = sum(normalized * local_overlap[group]),
        regional = sum(normalized * regional_overlap[group]),
        broad = sum(normalized * broad_overlap[group])
      )

      source_means[source] <- estimate
      source_variances[source] <- variance
      source_effective_n[source] <- effective_n
      source_support[source] <- mean(raw_weight)
      source_maximum[source] <- max(raw_weight)
      for (block in names(source_block_weighted)) {
        source_block_weighted[[block]][source] <- weighted_block[[block]]
        source_overlap_weighted[[block]][source] <- weighted_overlap[[block]]
      }
      mixture_values <- c(mixture_values, response)
      mixture_weights <- c(mixture_weights, normalized)

      this_sources[[source]] <- data.frame(
        target_id = target_id, interval_code = depth$interval_code,
        donor_source = source, prediction_kg_m2 = estimate,
        donor_weighted_sd_kg_m2 = sqrt(variance),
        donor_p10_kg_m2 = quantiles[[1L]], donor_median_kg_m2 = quantiles[[2L]],
        donor_p90_kg_m2 = quantiles[[3L]], n_donors = length(group),
        effective_donors = effective_n,
        combined_similarity_mean = mean(raw_weight),
        combined_similarity_max = max(raw_weight),
        local_similarity_weighted = weighted_block[["local"]],
        regional_similarity_weighted = weighted_block[["regional"]],
        broad_similarity_weighted = weighted_block[["broad"]],
        local_overlap_fraction_weighted = weighted_overlap[["local"]],
        regional_overlap_fraction_weighted = weighted_overlap[["regional"]],
        broad_overlap_fraction_weighted = weighted_overlap[["broad"]],
        stringsAsFactors = FALSE
      )

      if (keep_details) {
        selected <- head(order(raw_weight, decreasing = TRUE), top_n)
        this_donors[[source]] <- data.frame(
          target_id = target_id, interval_code = depth$interval_code,
          donor_source = source,
          donor_rank_within_source = seq_along(selected),
          donor_core_uid = stock$core_uid[group[selected]],
          donor_stock_kg_m2 = response[selected],
          local_similarity = local_similarity[group[selected]],
          regional_similarity = regional_similarity[group[selected]],
          broad_similarity = broad_similarity[group[selected]],
          local_standardized_distance = local_distance[group[selected]],
          regional_standardized_distance = regional_distance[group[selected]],
          broad_standardized_distance = broad_distance[group[selected]],
          local_overlap_variables = local_overlap_n[group[selected]],
          regional_overlap_variables = regional_overlap_n[group[selected]],
          broad_overlap_variables = broad_overlap_n[group[selected]],
          local_overlap_fraction = local_overlap[group[selected]],
          regional_overlap_fraction = regional_overlap[group[selected]],
          broad_overlap_fraction = broad_overlap[group[selected]],
          combined_similarity = raw_weight[selected],
          normalized_weight_within_source = normalized[selected],
          stringsAsFactors = FALSE
        )
        for (block in c("local", "regional", "broad")) {
          sim <- switch(
            block, local = local_similarity[group],
            regional = regional_similarity[group], broad = broad_similarity[group]
          )
          dist <- switch(
            block, local = local_distance[group],
            regional = regional_distance[group], broad = broad_distance[group]
          )
          overlap <- switch(
            block, local = local_overlap[group],
            regional = regional_overlap[group], broad = broad_overlap[group]
          )
          this_support[[paste(source, block, sep = "|")]] <- data.frame(
            target_id = target_id, interval_code = depth$interval_code,
            donor_source = source, scale_block = block,
            n_donors = length(group), similarity_mean = mean(sim),
            similarity_max = max(sim),
            similarity_combined_weighted_mean = sum(normalized * sim),
            standardized_distance_mean = mean(dist),
            standardized_distance_min = min(dist),
            standardized_distance_combined_weighted_mean = sum(normalized * dist),
            variables_in_block = length(model$blocks[[block]]),
            minimum_overlap_variables = model$minimum_overlap[[block]],
            overlap_fraction_mean = mean(overlap),
            overlap_fraction_combined_weighted_mean = sum(normalized * overlap),
            donors_similarity_ge_0_5 = sum(sim >= 0.5),
            fraction_similarity_ge_0_5 = mean(sim >= 0.5),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    n_sources <- length(source_means)
    if (!n_sources) {
      prediction_rows[[d]] <- empty_prediction_row(target_id, depth, "NUMERICAL_WEIGHT_FAILURE")
      next
    }
    source_pool_weight <- 1 / n_sources
    mixture_weights <- mixture_weights * source_pool_weight
    mixture_quantiles <- weighted_quantile(mixture_values, mixture_weights)
    estimate <- mean(source_means)
    within_variance <- mean(source_variances)
    between_variance <- mean((source_means - estimate)^2)
    total_variance <- within_variance + between_variance
    combined_support <- mean(source_support)
    maximum_combined <- max(source_maximum)
    status <- if (maximum_combined < support_warning_threshold) "LOW_SUPPORT" else "OK"

    prediction_rows[[d]] <- data.frame(
      target_id = target_id, interval_code = depth$interval_code,
      top_cm = depth$top_cm, bottom_cm = depth$bottom_cm,
      prediction_kg_m2 = estimate,
      mixture_p10_kg_m2 = mixture_quantiles[[1L]],
      mixture_median_kg_m2 = mixture_quantiles[[2L]],
      mixture_p90_kg_m2 = mixture_quantiles[[3L]],
      predictive_spread_sd_kg_m2 = sqrt(total_variance),
      within_source_sd_kg_m2 = sqrt(within_variance),
      between_source_sd_kg_m2 = sqrt(between_variance),
      n_sources = n_sources, n_donors = sum(lengths(source_groups[names(source_means)])),
      mean_effective_donors_per_source = mean(source_effective_n),
      local_support = mean(source_block_weighted$local),
      regional_support = mean(source_block_weighted$regional),
      broad_support = mean(source_block_weighted$broad),
      local_overlap_fraction = mean(source_overlap_weighted$local),
      regional_overlap_fraction = mean(source_overlap_weighted$regional),
      broad_overlap_fraction = mean(source_overlap_weighted$broad),
      combined_support = combined_support,
      maximum_combined_similarity = maximum_combined,
      local_exponent = depth$local, regional_exponent = depth$regional,
      broad_exponent = depth$broad, status = status,
      missing_covariates = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )

    if (keep_details) {
      source_table <- do.call(rbind, this_sources)
      source_table$equal_source_pool_weight <- source_pool_weight
      source_rows[[d]] <- source_table
      if (length(this_donors)) {
        donor_table <- do.call(rbind, this_donors)
        donor_table$equal_pool_contribution <-
          donor_table$normalized_weight_within_source * source_pool_weight
        donor_rows[[d]] <- donor_table
      }
      if (length(this_support)) {
        support_table <- do.call(rbind, this_support)
        pooled <- do.call(rbind, lapply(split(support_table, support_table$scale_block), function(z) {
          data.frame(
            target_id = target_id, interval_code = depth$interval_code,
            donor_source = "POOLED_EQUAL_SOURCE", scale_block = z$scale_block[[1L]],
            n_donors = sum(z$n_donors), similarity_mean = mean(z$similarity_mean),
            similarity_max = mean(z$similarity_max),
            similarity_combined_weighted_mean = mean(z$similarity_combined_weighted_mean),
            standardized_distance_mean = mean(z$standardized_distance_mean),
            standardized_distance_min = mean(z$standardized_distance_min),
            standardized_distance_combined_weighted_mean =
              mean(z$standardized_distance_combined_weighted_mean),
            variables_in_block = z$variables_in_block[[1L]],
            minimum_overlap_variables = z$minimum_overlap_variables[[1L]],
            overlap_fraction_mean = mean(z$overlap_fraction_mean),
            overlap_fraction_combined_weighted_mean =
              mean(z$overlap_fraction_combined_weighted_mean),
            donors_similarity_ge_0_5 = sum(z$donors_similarity_ge_0_5),
            fraction_similarity_ge_0_5 = mean(z$fraction_similarity_ge_0_5),
            stringsAsFactors = FALSE
          )
        }))
        support_rows[[d]] <- rbind(support_table, pooled)
      }
    }
  }

  list(
    predictions = do.call(rbind, prediction_rows),
    sources = if (length(source_rows)) do.call(rbind, source_rows) else NULL,
    donors = if (length(donor_rows)) do.call(rbind, donor_rows) else NULL,
    blocks = if (length(support_rows)) do.call(rbind, support_rows) else NULL
  )
}

bind_non_null <- function(items, name) {
  values <- lapply(items, `[[`, name)
  values <- values[!vapply(values, is.null, logical(1))]
  if (!length(values)) return(data.frame())
  rownames_out <- do.call(rbind, values)
  rownames(rownames_out) <- NULL
  rownames_out
}

predict_targets <- function(
    targets, id_column, model, depth_weights, min_source_donors,
    support_warning_threshold, top_n, keep_details = TRUE, progress = FALSE) {
  if (!nrow(targets)) {
    prediction_template <- empty_prediction_row(
      "__template__", depth_weights[1, , drop = FALSE], "TEMPLATE", ""
    )[0, , drop = FALSE]
    return(list(
      predictions = prediction_template,
      sources = data.frame(), donors = data.frame(), blocks = data.frame()
    ))
  }
  results <- vector("list", nrow(targets))
  for (i in seq_len(nrow(targets))) {
    if (progress && (i == 1L || i %% 250L == 0L || i == nrow(targets))) {
      cat("Predicting target ", i, " of ", nrow(targets), "\n", sep = "")
    }
    results[[i]] <- predict_one_target(
      targets[i, , drop = FALSE], as.character(targets[[id_column]][[i]]),
      model, depth_weights, min_source_donors, support_warning_threshold,
      top_n, keep_details
    )
  }
  list(
    predictions = bind_non_null(results, "predictions"),
    sources = bind_non_null(results, "sources"),
    donors = bind_non_null(results, "donors"),
    blocks = bind_non_null(results, "blocks")
  )
}

safe_write_csv <- function(data, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop("Output already exists (use --overwrite=true to replace Part 16 only): ", path, call. = FALSE)
  }
  utils::write.csv(data, path, row.names = FALSE, na = "")
}

safe_write_lines <- function(lines, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop("Output already exists (use --overwrite=true to replace Part 16 only): ", path, call. = FALSE)
  }
  writeLines(lines, path)
}

map_small_raster <- function(
    raster_path, output_path, model, depth_weights, min_source_donors,
    support_warning_threshold, overwrite = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("terra is required only when --raster is supplied.", call. = FALSE)
  }
  if (!file.exists(raster_path)) stop("Raster does not exist: ", raster_path, call. = FALSE)
  if (file.exists(output_path) && !overwrite) {
    stop("Raster output already exists; use --overwrite=true to replace Part 16 only.", call. = FALSE)
  }
  raster <- terra::rast(raster_path)
  absent <- setdiff(model$covariates, names(raster))
  if (length(absent)) {
    stop("Raster lacks model covariates: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  values <- terra::values(raster[[model$covariates]], mat = TRUE)
  eligible <- rep(TRUE, nrow(values))
  for (block in names(model$blocks)) {
    eligible <- eligible & rowSums(is.finite(values[, model$blocks[[block]], drop = FALSE])) >=
      model$minimum_overlap[[block]]
  }
  interval_output_names <- unlist(lapply(depth_weights$interval_code, function(code) {
    c(paste0("carbon_", code), paste0("support_", code), paste0("spread_", code))
  }))
  cumulative_output_names <- c("carbon_0_30cm", "carbon_0_100cm", "carbon_0_200cm")
  output_names <- c(interval_output_names, cumulative_output_names)
  output_values <- matrix(NA_real_, nrow(values), length(output_names), dimnames = list(NULL, output_names))
  cells <- which(eligible)
  if (!length(cells)) stop("Raster has no cells meeting the minimum overlap in every block.", call. = FALSE)
  for (j in seq_along(cells)) {
    if (j == 1L || j %% 250L == 0L || j == length(cells)) {
      cat("Mapping complete raster cell ", j, " of ", length(cells), "\n", sep = "")
    }
    cell <- cells[[j]]
    row <- as.data.frame(as.list(values[cell, model$covariates]), check.names = FALSE)
    row$core_uid <- paste0("RASTER_CELL:", cell)
    result <- predict_one_target(
      row, row$core_uid, model, depth_weights, min_source_donors,
      support_warning_threshold, top_n = 0L, keep_details = FALSE
    )$predictions
    for (d in seq_len(nrow(depth_weights))) {
      columns <- (d - 1L) * 3L + seq_len(3L)
      output_values[cell, columns] <- c(
        result$prediction_kg_m2[[d]], result$maximum_combined_similarity[[d]],
        result$predictive_spread_sd_kg_m2[[d]]
      )
    }
  }
  carbon_columns <- seq(1L, nrow(depth_weights) * 3L, by = 3L)
  output_values[, "carbon_0_30cm"] <- rowSums(output_values[, carbon_columns[1:2], drop = FALSE])
  output_values[, "carbon_0_100cm"] <- rowSums(output_values[, carbon_columns[1:4], drop = FALSE])
  output_values[, "carbon_0_200cm"] <- rowSums(output_values[, carbon_columns[1:5], drop = FALSE])
  layers <- lapply(seq_len(ncol(output_values)), function(j) {
    layer <- raster[[1L]]
    terra::values(layer) <- output_values[, j]
    layer
  })
  output <- do.call(c, layers)
  names(output) <- output_names
  terra::writeRaster(output, output_path, overwrite = overwrite, datatype = "FLT4S")
  output_path
}

run_self_test <- function() {
  set.seed(470316)
  n_a <- 20L
  n_b <- 4L
  ids <- c(paste0("A", seq_len(n_a)), paste0("B", seq_len(n_b)), "TARGET_NEAR", "TARGET_FAR")
  role <- c(rep("donor", n_a + n_b), "target", "target")
  base <- c(seq(-1, 1, length.out = n_a), seq(-0.8, 0.8, length.out = n_b), 0, 8)
  points <- data.frame(
    core_uid = ids, record_role = role,
    local_wetness = base,
    local_surface_variability = base * 0.6 + rep(c(-0.2, 0.2), length.out = length(base)),
    regional_river_distance = base * 0.7 + rep(c(-0.1, 0.1), length.out = length(base)),
    broad_climate = base * 0.4 + rep(c(-0.2, 0, 0.2), length.out = length(base)),
    broad_rare_gradient = c(rep(0, n_a + n_b - 1L), 5, 0, 0),
    stringsAsFactors = FALSE
  )
  points$local_surface_variability[points$core_uid %in% c("A1", "TARGET_NEAR")] <- NA_real_
  shard_dir <- tempfile("hierarchical_similarity_shards_")
  dir.create(shard_dir)
  utils::write.csv(
    points[seq_len(13L), ],
    file.path(shard_dir, "HierarchicalSimilarity_DonorCovariates_A.csv"), row.names = FALSE
  )
  utils::write.csv(
    points[14:nrow(points), ],
    file.path(shard_dir, "HierarchicalSimilarity_DonorCovariates_B.csv"), row.names = FALSE
  )
  utils::write.csv(data.frame(note = "must be ignored"), file.path(shard_dir, "manifest.csv"), row.names = FALSE)
  sharded_points <- read_point_csvs(shard_dir)
  if (nrow(sharded_points) != nrow(points) || attr(sharded_points, "point_shard_count") != 2L) {
    stop("Self-test failed: point-shard directory ingestion was not selective and complete.", call. = FALSE)
  }
  categorical_test <- points
  categorical_test$local_landform_categorical_code <- 1L
  rejected_categorical <- tryCatch({
    prepare_points(
      categorical_test, "core_uid", "record_role",
      c("target", "community", "validation", "prediction")
    )
    FALSE
  }, error = function(e) grepl("one-hot", conditionMessage(e), ignore.case = TRUE))
  if (!rejected_categorical) {
    stop("Self-test failed: a raw categorical code entered continuous distance.", call. = FALSE)
  }
  stocks <- data.frame(
    core_uid = ids[seq_len(n_a + n_b)],
    database_code = c(rep("SOURCE_A", n_a), rep("SOURCE_B", n_b)),
    stringsAsFactors = FALSE
  )
  for (d in seq_len(nrow(DEPTH_SPEC))) {
    code <- DEPTH_SPEC$interval_code[[d]]
    stocks[[paste0("carbon_stock_", code, "_kg_m2")]] <-
      c(rep(10 * d, n_a), rep(30 * d, n_b))
    stocks[[paste0("carbon_complete_", code)]] <- TRUE
  }
  stocks$carbon_complete_100_200cm[stocks$core_uid == "B4"] <- FALSE

  points_info <- prepare_points(
    points, "core_uid", "record_role", c("target", "community", "validation", "prediction")
  )
  stock_long <- prepare_stock_long(stocks, "core_uid", FALSE)
  model <- prepare_model(points_info, stock_long, "core_uid", min_block_overlap_fraction = 0.50)
  result <- predict_targets(
    points[points_info$target, ], "core_uid", model, DEPTH_SPEC,
    min_source_donors = 3L, support_warning_threshold = 0.10,
    top_n = 3L, keep_details = TRUE
  )
  near <- result$predictions[result$predictions$target_id == "TARGET_NEAR", ]
  far <- result$predictions[result$predictions$target_id == "TARGET_FAR", ]
  if (nrow(near) != 5L || any(abs(near$prediction_kg_m2 - 20 * seq_len(5L)) > 1e-8)) {
    stop("Self-test failed: equal-source pooling is not invariant to source sample size.", call. = FALSE)
  }
  if (!all(near$n_sources == 2L)) stop("Self-test failed: expected two pooled sources.", call. = FALSE)
  if (any(!near$status %in% c("OK", "LOW_SUPPORT")) ||
      any(abs(near$local_overlap_fraction - 0.5) > 1e-8)) {
    stop("Self-test failed: a permissible missing target band was not handled pairwise.", call. = FALSE)
  }
  rare_scale <- model$scaling[model$scaling$covariate == "broad_rare_gradient", ]
  if (nrow(rare_scale) != 1L || rare_scale$scale_method != "SD_fallback") {
    stop("Self-test failed: zero-MAD covariate did not use the donor SD fallback.", call. = FALSE)
  }
  if (!all(far$maximum_combined_similarity < near$maximum_combined_similarity)) {
    stop("Self-test failed: cross-scale support did not fall for an extrapolative target.", call. = FALSE)
  }
  deep_sources <- result$sources[result$sources$target_id == "TARGET_NEAR" &
                                   result$sources$interval_code == "100_200cm", ]
  if (deep_sources$n_donors[deep_sources$donor_source == "SOURCE_B"] != 3L) {
    stop("Self-test failed: strict depth completeness was not applied.", call. = FALSE)
  }
  if (!nrow(result$donors) || !nrow(result$blocks)) {
    stop("Self-test failed: donor or block diagnostics are empty.", call. = FALSE)
  }
  cat("Self-test passed: robust scaling (including SD fallback), pairwise missing-band overlap,\n")
  cat("product similarity, strict depth completeness, equal-source pooling, low-support\n")
  cat("detection, and diagnostic tables all behaved as expected.\n")
  invisible(TRUE)
}

main <- function() {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (isTRUE(cli$help)) {
    usage()
    return(invisible(TRUE))
  }
  if (isTRUE(cli$self_test)) return(run_self_test())

  points_path <- first_nonempty(cli$points, Sys.getenv("HMS_POINTS_CSV", unset = ""))
  stocks_path <- first_nonempty(cli$stocks, Sys.getenv("HMS_STOCKS_CSV", unset = ""))
  if (is.null(points_path) || is.null(stocks_path)) {
    usage()
    stop("--points and --stocks are required outside self-test mode.", call. = FALSE)
  }
  output_dir <- first_nonempty(
    cli$output, Sys.getenv("HMS_OUTPUT_DIR", unset = ""), DEFAULT_OUTPUT_DIR
  )
  weights_path <- first_nonempty(cli$weights, Sys.getenv("HMS_WEIGHTS_CSV", unset = ""))
  raster_path <- first_nonempty(cli$raster, Sys.getenv("HMS_RASTER", unset = ""))
  id_column <- first_nonempty(cli$id_column, "core_uid")
  role_column <- first_nonempty(cli$role_column)
  target_roles <- strsplit(
    first_nonempty(cli$target_roles, "target,community,validation,prediction"), ",", fixed = TRUE
  )[[1L]]
  top_n <- as.integer(first_nonempty(cli$top_n, "5"))
  min_source_donors <- as.integer(first_nonempty(cli$min_source_donors, "3"))
  min_block_overlap_fraction <- as.numeric(first_nonempty(cli$min_block_overlap_fraction, "0.50"))
  support_warning_threshold <- as.numeric(first_nonempty(cli$support_warning_threshold, "0.10"))
  allow_inferred <- parse_boolean(cli$allow_inferred_completeness, FALSE)
  overwrite <- parse_boolean(cli$overwrite, FALSE)
  if (!is.finite(top_n) || top_n < 1L || !is.finite(min_source_donors) || min_source_donors < 1L ||
      !is.finite(min_block_overlap_fraction) || min_block_overlap_fraction <= 0 ||
      min_block_overlap_fraction > 1 ||
      !is.finite(support_warning_threshold) || support_warning_threshold < 0 ||
      support_warning_threshold > 1) {
    stop("Invalid top-n, minimum donor, block-overlap, or support threshold setting.", call. = FALSE)
  }

  points <- read_point_csvs(points_path)
  point_shard_count <- attr(points, "point_shard_count")
  point_shard_manifest <- attr(points, "point_shard_manifest")
  stocks <- read_csv(stocks_path)
  depth_weights <- load_depth_weights(weights_path)
  points_info <- prepare_points(
    points, id_column, role_column, target_roles,
    allow_no_targets = !is.null(raster_path)
  )
  stock_long <- prepare_stock_long(stocks, id_column, allow_inferred)
  model <- prepare_model(
    points_info, stock_long, id_column,
    min_block_overlap_fraction = min_block_overlap_fraction
  )
  targets <- points_info$data[points_info$target, , drop = FALSE]
  results <- predict_targets(
    targets, id_column, model, depth_weights, min_source_donors,
    support_warning_threshold, top_n, keep_details = TRUE, progress = nrow(targets) > 50L
  )

  if (grepl("(^|[/\\\\])route_l([/\\\\]|$)", output_dir, ignore.case = TRUE)) {
    stop("Refusing to write this research pilot inside outputs/route_l.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  diagnostics <- data.frame(
    diagnostic = c(
      "point_csv_shards", "points_input", "targets_predicted", "stock_rows_input",
      "complete_stock_rows", "matched_donor_cores", "unmatched_stock_core_ids",
      "local_covariates_used", "regional_covariates_used", "broad_covariates_used",
      "minimum_block_overlap_fraction", "covariates_using_sd_fallback",
      "prediction_rows_ok", "prediction_rows_low_support", "prediction_rows_failed",
      "strict_completeness_flags", "equal_source_pooling", "route_l_modified"
    ),
    value = c(
      point_shard_count, nrow(points), nrow(targets), nrow(stock_long),
      sum(stock_long$complete_interval),
      length(model$donor_ids), length(model$unmatched_stock_ids),
      length(model$blocks$local), length(model$blocks$regional), length(model$blocks$broad),
      min_block_overlap_fraction, sum(model$scaling$scale_method == "SD_fallback"),
      sum(results$predictions$status == "OK"), sum(results$predictions$status == "LOW_SUPPORT"),
      sum(!results$predictions$status %in% c("OK", "LOW_SUPPORT")),
      !allow_inferred, TRUE, FALSE
    ),
    stringsAsFactors = FALSE
  )

  safe_write_csv(results$predictions, file.path(output_dir, "point_predictions.csv"), overwrite)
  safe_write_csv(results$sources, file.path(output_dir, "source_predictions.csv"), overwrite)
  safe_write_csv(results$donors, file.path(output_dir, "top_donors.csv"), overwrite)
  safe_write_csv(results$blocks, file.path(output_dir, "block_support.csv"), overwrite)
  safe_write_csv(model$scaling, file.path(output_dir, "robust_scaling.csv"), overwrite)
  safe_write_csv(depth_weights, file.path(output_dir, "depth_exponents.csv"), overwrite)
  safe_write_csv(diagnostics, file.path(output_dir, "run_diagnostics.csv"), overwrite)
  safe_write_csv(point_shard_manifest, file.path(output_dir, "input_point_shards.csv"), overwrite)

  raster_output <- NULL
  if (!is.null(raster_path)) {
    raster_output <- map_small_raster(
      raster_path, file.path(output_dir, "pilot_prediction_stack.tif"),
      model, depth_weights, min_source_donors, support_warning_threshold, overwrite
    )
  }
  safe_write_lines(c(
    "Hierarchical multi-scale similarity micro-pilot",
    "",
    "Research status: exploratory; not adopted for Route L.",
    "Distance within each block is the mean squared donor-target difference after donor-median scaling by MAD, with donor SD fallback only when MAD is unusable.",
    paste0("Distances use pairwise-complete variables and require at least ceiling(",
           min_block_overlap_fraction, " * retained variables), with an absolute minimum of one variable, in every block."),
    "Block similarity is exp(-0.5 * distance). Depth similarity is the product of block similarities raised to the prespecified depth exponents.",
    "Predictions are first calculated within each donor source and then pooled with equal source weights.",
    "Only explicit complete, finite, non-negative donor stocks are used at each depth unless inferred completeness was explicitly enabled.",
    "Reported spread is the standard deviation of the equal-source donor mixture; it is diagnostic, not a calibrated confidence or credible interval.",
    paste0("Support warning threshold (maximum combined donor similarity): ", support_warning_threshold, "."),
    paste0("Raster output: ", if (is.null(raster_output)) "not requested" else raster_output),
    "Raster carbon_0_30cm, carbon_0_100cm and carbon_0_200cm are sums of interval estimates; no cumulative uncertainty layer is claimed.",
    "No authoritative Route L file was read, replaced, or modified."
  ), file.path(output_dir, "RUN_NOTES.txt"), overwrite)

  cat("Completed isolated hierarchical similarity pilot.\n")
  cat("Targets: ", nrow(targets), "; donors: ", length(model$donor_ids), "\n", sep = "")
  cat("Output: ", normalizePath(output_dir, mustWork = FALSE), "\n", sep = "")
  invisible(results)
}

if (sys.nframe() == 0L) main()
