#!/usr/bin/env Rscript

# Evaluate hierarchical-similarity point predictions against measured community
# intervals. This is deliberately separate from model fitting: measured carbon
# never enters Earth Engine or the similarity calculation.

options(stringsAsFactors = FALSE)

DEFAULT_MODEL_DIR <- file.path(
  "outputs", "covariate_analysis", "part16_hierarchical_similarity_micro_pilot"
)
DEFAULT_TRUTH <- file.path(
  "outputs", "covariate_analysis", "part7_vertical_profile",
  "community_depth_segments.csv"
)

parse_cli <- function(args) {
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      out$help <- TRUE
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

safe_read <- function(path) {
  if (!file.exists(path)) stop("Input does not exist: ", path, call. = FALSE)
  utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "-9999", "-9999.0")
  )
}

metric_row <- function(data, group, interval) {
  error <- data$prediction_kg_m2 - data$observed_kg_m2
  n <- sum(is.finite(error))
  data.frame(
    validation_group = group,
    interval_code = interval,
    n = n,
    mean_observed_kg_m2 = if (n) mean(data$observed_kg_m2) else NA_real_,
    mean_predicted_kg_m2 = if (n) mean(data$prediction_kg_m2) else NA_real_,
    bias_kg_m2 = if (n) mean(error) else NA_real_,
    mae_kg_m2 = if (n) mean(abs(error)) else NA_real_,
    rmse_kg_m2 = if (n) sqrt(mean(error^2)) else NA_real_,
    correlation = if (n >= 3L && stats::sd(data$observed_kg_m2) > 0 &&
                      stats::sd(data$prediction_kg_m2) > 0) {
      stats::cor(data$observed_kg_m2, data$prediction_kg_m2)
    } else {
      NA_real_
    },
    mean_max_similarity = if (n && "maximum_combined_similarity" %in% names(data)) {
      mean(data$maximum_combined_similarity, na.rm = TRUE)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

main <- function() {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (isTRUE(cli$help)) {
    cat(paste0(
      "Community validation for hierarchical similarity\n\n",
      "Optional:\n",
      "  --predictions=PATH  Default: ",
      file.path(DEFAULT_MODEL_DIR, "point_predictions.csv"), "\n",
      "  --truth=PATH        Default: ", DEFAULT_TRUTH, "\n",
      "  --output=PATH       Default: ", DEFAULT_MODEL_DIR, "\n",
      "  --overwrite=true    Replace prior validation outputs.\n"
    ))
    return(invisible(TRUE))
  }

  predictions_path <- if (!is.null(cli$predictions)) cli$predictions else
    file.path(DEFAULT_MODEL_DIR, "point_predictions.csv")
  truth_path <- if (!is.null(cli$truth)) cli$truth else DEFAULT_TRUTH
  output_dir <- if (!is.null(cli$output)) cli$output else DEFAULT_MODEL_DIR
  overwrite <- !is.null(cli$overwrite) && truthy(cli$overwrite)

  predictions <- safe_read(predictions_path)
  truth <- safe_read(truth_path)
  required_prediction <- c(
    "target_id", "interval_code", "prediction_kg_m2", "status"
  )
  required_truth <- c(
    "core_uid", "plot_id", "interval_code", "target_kg_m2",
    "complete_interval"
  )
  if (!all(required_prediction %in% names(predictions))) {
    stop(
      "Prediction table is missing: ",
      paste(setdiff(required_prediction, names(predictions)), collapse = ", "),
      call. = FALSE
    )
  }
  if (!all(required_truth %in% names(truth))) {
    stop(
      "Truth table is missing: ",
      paste(setdiff(required_truth, names(truth)), collapse = ", "),
      call. = FALSE
    )
  }

  truth <- truth[
    truthy(truth$complete_interval) & is.finite(truth$target_kg_m2),
    required_truth,
    drop = FALSE
  ]
  names(truth)[names(truth) == "core_uid"] <- "target_id"
  names(truth)[names(truth) == "target_kg_m2"] <- "observed_kg_m2"
  joined <- merge(
    predictions, truth,
    by = c("target_id", "interval_code"), all = FALSE, sort = TRUE
  )
  joined <- joined[
    joined$status %in% c("OK", "LOW_SUPPORT") &
      is.finite(joined$prediction_kg_m2) &
      is.finite(joined$observed_kg_m2),
    , drop = FALSE
  ]
  if (!nrow(joined)) {
    stop("No complete community observations matched usable predictions.", call. = FALSE)
  }
  joined$residual_kg_m2 <- joined$prediction_kg_m2 - joined$observed_kg_m2
  joined$validation_group <- ifelse(
    grepl("^FS", joined$plot_id, ignore.case = TRUE), "FS_ONLY", "OTHER_COMMUNITY"
  )

  # Build 0-30 cm only where both measured intervals and both predictions exist.
  shallow <- joined[joined$interval_code %in% c("0_15cm", "15_30cm"), ]
  split_shallow <- split(shallow, shallow$target_id)
  cumulative <- lapply(split_shallow, function(x) {
    if (!all(c("0_15cm", "15_30cm") %in% x$interval_code)) return(NULL)
    x <- x[match(c("0_15cm", "15_30cm"), x$interval_code), ]
    data.frame(
      target_id = x$target_id[[1L]],
      plot_id = x$plot_id[[1L]],
      interval_code = "0_30cm",
      prediction_kg_m2 = sum(x$prediction_kg_m2),
      observed_kg_m2 = sum(x$observed_kg_m2),
      status = if (any(x$status == "LOW_SUPPORT")) "LOW_SUPPORT" else "OK",
      maximum_combined_similarity = min(
        x$maximum_combined_similarity, na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })
  cumulative <- cumulative[!vapply(cumulative, is.null, logical(1))]
  cumulative <- if (length(cumulative)) do.call(rbind, cumulative) else data.frame()
  if (nrow(cumulative)) {
    cumulative$residual_kg_m2 <-
      cumulative$prediction_kg_m2 - cumulative$observed_kg_m2
    cumulative$validation_group <- ifelse(
      grepl("^FS", cumulative$plot_id, ignore.case = TRUE),
      "FS_ONLY", "OTHER_COMMUNITY"
    )
  }

  metric_data <- joined
  if (nrow(cumulative)) {
    shared <- union(names(metric_data), names(cumulative))
    for (name in setdiff(shared, names(metric_data))) metric_data[[name]] <- NA
    for (name in setdiff(shared, names(cumulative))) cumulative[[name]] <- NA
    metric_data <- rbind(metric_data[, shared], cumulative[, shared])
  }
  groups <- list(
    ALL_COMMUNITY = metric_data,
    FS_ONLY = metric_data[metric_data$validation_group == "FS_ONLY", , drop = FALSE]
  )
  metrics <- do.call(rbind, lapply(names(groups), function(group_name) {
    group_data <- groups[[group_name]]
    intervals <- unique(group_data$interval_code)
    do.call(rbind, lapply(intervals, function(interval) {
      metric_row(
        group_data[group_data$interval_code == interval, , drop = FALSE],
        group_name, interval
      )
    }))
  }))
  rownames(metrics) <- NULL

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- list(
    community_validation_rows = joined,
    community_validation_0_30cm = cumulative,
    community_validation_metrics = metrics
  )
  for (name in names(outputs)) {
    path <- file.path(output_dir, paste0(name, ".csv"))
    if (file.exists(path) && !overwrite) {
      stop("Output exists; use --overwrite=true: ", path, call. = FALSE)
    }
    utils::write.csv(outputs[[name]], path, row.names = FALSE, na = "")
  }

  comparison_plot <- file.path(
    output_dir, "community_validation_observed_vs_predicted.png"
  )
  support_plot <- file.path(
    output_dir, "community_prediction_support_heatmap.png"
  )
  if ((!file.exists(comparison_plot) && !file.exists(support_plot)) || overwrite) {
    grDevices::png(comparison_plot, width = 1800, height = 900, res = 170)
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    interval_titles <- c("0_15cm" = "0–15 cm", "15_30cm" = "15–30 cm")
    for (interval in names(interval_titles)) {
      plot_data <- joined[joined$interval_code == interval, , drop = FALSE]
      limits <- range(
        c(plot_data$observed_kg_m2, plot_data$prediction_kg_m2),
        finite = TRUE
      )
      padding <- max(diff(limits) * 0.08, 0.25)
      limits <- limits + c(-padding, padding)
      graphics::plot(
        plot_data$observed_kg_m2, plot_data$prediction_kg_m2,
        xlim = limits, ylim = limits, asp = 1,
        xlab = "Observed carbon stock (kg C/m²)",
        ylab = "Similarity prediction (kg C/m²)",
        main = interval_titles[[interval]],
        pch = ifelse(plot_data$validation_group == "FS_ONLY", 19, 17),
        col = ifelse(plot_data$status == "LOW_SUPPORT", "#D73027", "#2166AC"),
        cex = 1.25
      )
      graphics::abline(0, 1, lty = 2, col = "#555555", lwd = 2)
      graphics::text(
        plot_data$observed_kg_m2, plot_data$prediction_kg_m2,
        labels = plot_data$plot_id, pos = 3, cex = 0.68
      )
      graphics::legend(
        "topleft", c("Usable support", "Low support"),
        col = c("#2166AC", "#D73027"), pch = 19, bty = "n", cex = 0.8
      )
    }
    grDevices::dev.off()

    support_data <- predictions
    target_order <- sort(unique(support_data$target_id))
    depth_order <- c("0_15cm", "15_30cm", "30_50cm", "50_100cm", "100_200cm")
    support_matrix <- matrix(
      NA_real_, nrow = length(target_order), ncol = length(depth_order),
      dimnames = list(sub("^COMMUNITY:", "", target_order), depth_order)
    )
    for (i in seq_len(nrow(support_data))) {
      row <- match(support_data$target_id[[i]], target_order)
      column <- match(support_data$interval_code[[i]], depth_order)
      if (!is.na(row) && !is.na(column)) {
        support_matrix[row, column] <-
          support_data$maximum_combined_similarity[[i]]
      }
    }
    palette <- grDevices::colorRampPalette(
      c("#B2182B", "#FDDDBC", "#D1E5F0", "#2166AC")
    )(100)
    grDevices::png(support_plot, width = 1500, height = 900, res = 170)
    graphics::par(mar = c(6, 5, 3, 5))
    graphics::image(
      x = seq_len(ncol(support_matrix)),
      y = seq_len(nrow(support_matrix)),
      z = t(support_matrix),
      zlim = c(0, 1), col = palette, axes = FALSE,
      xlab = "Depth interval", ylab = "",
      main = "Maximum donor similarity by core and depth"
    )
    graphics::axis(
      1, at = seq_len(ncol(support_matrix)),
      labels = c("0–15", "15–30", "30–50", "50–100", "100–200"),
      las = 2
    )
    graphics::axis(
      2, at = seq_len(nrow(support_matrix)),
      labels = rownames(support_matrix), las = 2
    )
    for (row in seq_len(nrow(support_matrix))) {
      for (column in seq_len(ncol(support_matrix))) {
        value <- support_matrix[row, column]
        if (is.finite(value)) {
          graphics::text(
            column, row, sprintf("%.2f", value),
            col = if (value < 0.25 || value > 0.75) "white" else "black",
            cex = 0.72
          )
        }
      }
    }
    graphics::mtext(
      "Red indicates weak environmental analogues; blue indicates stronger support",
      side = 1, line = 4.7, cex = 0.8
    )
    grDevices::dev.off()
  }

  cat("Community validation complete.\n")
  print(metrics, row.names = FALSE)
  cat("Output: ", normalizePath(output_dir, mustWork = FALSE), "\n", sep = "")
  invisible(metrics)
}

if (sys.nframe() == 0L) main()
