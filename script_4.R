
# ==========================

# Part 4: prepare the prior and compare mapped stocks

# ===============
library(terra)
library(readr)
library(dplyr)
source("config.R")
source("helpers/grid_helper.R")

# 1. Inputs and the shared analysis grid
output_folder <- "outputs/script4"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
shared <- load_analysis_grid()
template <- shared$template
analysis_mask <- shared$mask

# Figure settings: NULL uses one common range across all methods AND depths.
# Set numeric limits here if you want to fix the scale between workshop runs.
map_colours <- hcl.colors(100, "Viridis")
prediction_limits <- NULL
sd_limits <- NULL
panel_columns <- 3

# Transfer external values to the shared cell centres once.
align_to_grid <- function(x) {
  if (!compareGeom(x, template, stopOnError = FALSE)) {
    x <- project(x, template, method = "near")
  }
  mask(x, analysis_mask)
}

# 2. Prepare paired prior stock and uncertainty
prepared_priors <- list()
for (depth in intersect(names(prior_mean_files), names(prior_width_files))) {
  source_mean <- rast(prior_mean_files[[depth]])
  source_width <- rast(prior_width_files[[depth]])
  compareGeom(source_mean, source_width)
  pair <- align_to_grid(c(source_mean, source_width))
  names(pair) <- c("prior_mean", "interval_width_90")
  prior_sd <- pair$interval_width_90 / (2 * qnorm(0.95))
  usable <- is.finite(pair$prior_mean) & pair$prior_mean > 0 &
    is.finite(prior_sd) & prior_sd > 0
  pair <- mask(pair, usable, maskvalues = 0)
  prior_sd <- mask(prior_sd, usable, maskvalues = 0)
  log_variance <- log1p((prior_sd / pair$prior_mean)^2)
  log_mean <- log(pair$prior_mean) - log_variance / 2
  prior <- c(pair$prior_mean, prior_sd, log_mean, log_variance)
  names(prior) <- c("prior_mean", "prior_sd", "log_mean", "log_variance")
  writeRaster(prior, file.path(output_folder, paste0("prior_moments_0_", depth, "cm.tif")),
              datatype = "FLT8S", overwrite = TRUE)
  # Keep the familiar filenames; the width layer retains its original units.
  writeRaster(pair$prior_mean, file.path(output_folder, basename(prior_mean_files[[depth]])),
              overwrite = TRUE)
  writeRaster(pair$interval_width_90,
              file.path(output_folder, basename(prior_width_files[[depth]])), overwrite = TRUE)
  prepared_priors[[depth]] <- prior
  cat("\nPrior 0–", depth, "cm: Q95-Q05 -> approximate SD -> lognormal moments.\n", sep = "")
  print(global(prior, c("min", "max", "notNA"), na.rm = TRUE))
  png(file.path(output_folder, paste0("prior_0_", depth, "cm.png")),
      width = 1400, height = 700, res = 150)
  plot(prior[[1:2]], col = map_colours,
       main = c("Prior mean (kg C/m²)", "Approximate prior SD (kg C/m²)"))
  dev.off()
  plot(prior[[1:2]], col = map_colours,
       main = c("Prior mean", "Approximate prior SD"))
}

# 3. Read the saved prediction and uncertainty products
# Everything below uses the shared grid. No field stocks or fitted models change.
make_map_record <- function(prediction, variance, method, uncertainty, scale,
                            n_cores = NA_integer_, selected = FALSE) {
  if (is.null(variance)) {
    variance <- rast(template)
    values(variance) <- NA_real_
    uncertainty <- paste("Unavailable. Expected definition:", uncertainty)
  }
  list(prediction = prediction, variance = variance, info = data.frame(
    method = method, modelling_scale = scale, uncertainty_definition = uncertainty,
    n_cores = n_cores, selected_regression = selected,
    weighting_power = NA_real_, neighbours = NA_real_,
    CV_error_log_variance = NA_real_, prediction_log_variance = NA_real_,
    CV_RMSE_log = NA_real_, CV_RMSE_kg_C_m2 = NA_real_, AICc = NA_real_
  ))
}

# Regression uncertainty is derived from the EXISTING lm objects in Script 3.
# For a new point: prediction variance = SE(fitted mean)^2 + residual variance.
# This assumes the fitted raw-stock model with independent, constant-variance
# errors. It does not include model-selection or spatial error uncertainty.
covariate_file <- "outputs/script3/covariates/terrain_water.tif"
grid_covariates <- NULL
if (file.exists(covariate_file)) {
  grid_covariates <- as.data.frame(
    values(align_to_grid(rast(covariate_file)), mat = TRUE)
  )
}

comparisons <- list()
summary_rows <- list()
guide_rows <- list()
depth_groups <- c(paste0("0_", mapping_depths_cm), "Total")

for (depth in depth_groups) {
  suffix <- if (depth == "Total") "total" else paste0(depth, "cm")
  methods <- list()

  # Area-based, IDW and kriging: read carbon-scale variances already exported.
  for (method in c("area_based", "idw", "ok")) {
    prediction_file <- file.path("outputs/script2", paste0(method, "_", suffix, ".tif"))
    variance_file <- file.path("outputs/script2", paste0(method, "_variance_", suffix, ".tif"))
    if (!file.exists(prediction_file)) {
      cat("No saved", method, "prediction for", depth, "- skipped.\n")
      next
    }
    prediction <- rast(prediction_file)
    compareGeom(prediction, template)
    variance <- NULL
    if (file.exists(variance_file)) {
      variance <- rast(variance_file)
      compareGeom(variance, template)
      variance <- ifel(variance >= 0, variance, NA)
    } else {
      cat("No saved uncertainty for", method, depth, "\n")
    }
    uncertainty <- switch(method,
      area_based = paste(
        "Lognormal approximation from constant-mean normal-model prediction variance:",
        "var(log stock)*(1+1/n); independent cores and fixed estimated parameters"
      ),
      idw = "Lognormal approximation from pooled held-out plot errors; working error variance and fixed calibration",
      ok = "Lognormal approximation from conditional kriging variance; inspect variogram warnings"
    )
    record <- make_map_record(
      mask(prediction, analysis_mask), variance, method, uncertainty, "Natural log carbon"
    )
    calibration_file <- file.path("outputs/script2", paste0(method, "_calibration_", suffix, ".csv"))
    if (file.exists(calibration_file)) {
      calibration <- read_csv(calibration_file, show_col_types = FALSE)
      record$info$n_cores <- calibration$n_cores[1]
      record$info$CV_RMSE_log <- calibration$RMSE_log[1]
      record$info$weighting_power <- calibration$power[1]
      record$info$neighbours <- calibration$neighbours[1]
      record$info$CV_error_log_variance <- calibration$log_error_variance[1]
      if ("prediction_log_variance" %in% names(calibration)) {
        record$info$prediction_log_variance <- calibration$prediction_log_variance[1]
      }
    }
    methods[[method]] <- record
  }

  # Include every saved Script 3 regression, marking the selected model.
  regression_folder <- file.path("outputs/script3", paste0("stock_", depth))
  models_file <- file.path(regression_folder, "fitted_models.rds")
  selected_file <- file.path(regression_folder, "selected_model.csv")
  evaluation_file <- file.path(regression_folder, "model_evaluation.csv")
  selected_model <- if (file.exists(selected_file)) {
    read_csv(selected_file, show_col_types = FALSE)$model[1]
  } else NA_character_
  evaluation <- if (file.exists(evaluation_file)) {
    read_csv(evaluation_file, show_col_types = FALSE)
  } else NULL

  if (file.exists(models_file)) {
    fitted_models <- readRDS(models_file)
    for (model_name in names(fitted_models)) {
      prediction_file <- file.path(regression_folder, model_name, "carbon_stock_kg_m2.tif")
      if (!file.exists(prediction_file)) next
      prediction <- align_to_grid(rast(prediction_file))
      fit <- fitted_models[[model_name]]
      variance <- NULL
      predictors <- all.vars(formula(fit))[-1]
      if (!is.null(grid_covariates) || !length(predictors)) {
        newdata <- if (length(predictors)) {
          grid_covariates[, predictors, drop = FALSE]
        } else data.frame(cell = seq_len(ncell(template)))
        fitted <- stats::predict(fit, newdata = newdata, se.fit = TRUE)
        variance <- rast(template)
        values(variance) <- as.numeric(fitted$se.fit)^2 + stats::sigma(fit)^2
        variance <- mask(variance, prediction)
      }
      method <- paste0("Regression_", model_name)
      record <- make_map_record(
        prediction, variance, method,
        "Normal-model prediction variance: SE(fitted mean)^2 + residual variance; no selection/spatial error correction",
        "Original carbon units (Script 3 unchanged)", nobs(fit),
        selected = isTRUE(model_name == selected_model)
      )
      if (!is.null(evaluation)) {
        score <- evaluation[evaluation$model == model_name, ]
        if (nrow(score)) {
          record$info$CV_RMSE_kg_C_m2 <- score$RMSE[1]
          record$info$AICc <- score$AICc[1]
        }
      }
      methods[[method]] <- record
      if (isTRUE(model_name == selected_model)) {
        # Preserve the existing aligned selected-regression filename.
        writeRaster(prediction,
                    file.path(output_folder, paste0("regression_selected_", suffix, ".tif")),
                    overwrite = TRUE)
      }
    }
  }

  # Prior uncertainty is approximate stock SD, not the original Q95-Q05 width.
  depth_key <- sub("^0_", "", depth)
  if (depth_key %in% names(prepared_priors)) {
    prior <- prepared_priors[[depth_key]]
    methods$Prior <- make_map_record(
      prior$prior_mean, prior$prior_sd^2, "Prior",
      "Source Q95-Q05 width divided by 2*qnorm(0.95) gives approximate SD",
      "Source carbon mean; moment-matched lognormal prior"
    )
  }
  if (!length(methods)) next

  # Use matching-depth predictions on the SAME valid cells for each comparison.
  predictions <- do.call(c, unname(lapply(methods, function(x) x$prediction)))
  variances <- do.call(c, unname(lapply(methods, function(x) x$variance)))
  names(predictions) <- names(variances) <- names(methods)
  common_cells <- app(is.finite(predictions), "sum") == nlyr(predictions)
  predictions <- mask(predictions, common_cells, maskvalues = 0)
  variances <- mask(variances, common_cells, maskvalues = 0)
  standard_deviations <- sqrt(variances)
  comparisons[[depth]] <- list(prediction = predictions, sd = standard_deviations)

  # Export prediction, variance and SD stacks for GIS use.
  writeRaster(predictions, file.path(output_folder, paste0("comparison_", depth, ".tif")),
              overwrite = TRUE)
  writeRaster(variances, file.path(output_folder, paste0("uncertainty_variance_", depth, ".tif")),
              overwrite = TRUE)
  writeRaster(standard_deviations, file.path(output_folder, paste0("uncertainty_sd_", depth, ".tif")),
              overwrite = TRUE)

  for (method in names(methods)) {
    prediction_values <- as.numeric(values(predictions[[method]]))
    prediction_values <- prediction_values[is.finite(prediction_values)]
    sd_values <- as.numeric(values(standard_deviations[[method]]))
    sd_values <- sd_values[is.finite(sd_values)]
    key <- paste(depth, method, sep = ":")
    summary_rows[[key]] <- data.frame(
      depth = depth, method = method,
      prediction_cells = length(prediction_values),
      mean_kg_C_m2 = if (length(prediction_values)) mean(prediction_values) else NA_real_,
      min_kg_C_m2 = if (length(prediction_values)) min(prediction_values) else NA_real_,
      max_kg_C_m2 = if (length(prediction_values)) max(prediction_values) else NA_real_,
      negative_prediction_cells = sum(prediction_values < 0),
      spatial_SD_kg_C_m2 = if (length(prediction_values) > 1) sd(prediction_values) else NA_real_,
      uncertainty_cells = length(sd_values),
      mean_pointwise_SD_kg_C_m2 = if (length(sd_values)) mean(sd_values) else NA_real_,
      min_pointwise_SD_kg_C_m2 = if (length(sd_values)) min(sd_values) else NA_real_,
      max_pointwise_SD_kg_C_m2 = if (length(sd_values)) max(sd_values) else NA_real_,
      mean_pointwise_variance = if (length(sd_values)) mean(sd_values^2) else NA_real_
    )
    guide_rows[[key]] <- cbind(depth = depth, methods[[method]]$info)
  }
}

# 4. Summary table: map variation and pointwise uncertainty are separate columns.
map_summary <- bind_rows(summary_rows) %>%
  left_join(bind_rows(guide_rows), by = c("depth", "method"))
if (!nrow(map_summary)) stop("No saved prediction maps found. Run Scripts 2 and 3 first.")
write_csv(map_summary, file.path(output_folder, "map_summary.csv"))
write_csv(bind_rows(guide_rows), file.path(output_folder, "map_methods.csv"))
cat("\nMap summary: means/SDs are kg C/m²; variances are (kg C/m²)².\n")
print(as.data.frame(map_summary[, c(
  "depth", "method", "prediction_cells", "mean_kg_C_m2", "min_kg_C_m2",
  "max_kg_C_m2", "mean_pointwise_SD_kg_C_m2"
)]), row.names = FALSE, digits = 4)
cat("Spatial SD describes variation across predicted cells; pointwise SD describes uncertainty.\n")
cat("Log RMSE and carbon-scale RMSE are separate columns and must not be ranked together.\n")
cat("Uncertainty definitions differ by method; these maps do not measure regional-total uncertainty.\n")
cat("Only saved results are shown; old outputs require rerunning the producing script to refresh.\n")

# 5. Plot prediction and uncertainty panels with consistent colours and limits.
if (is.null(prediction_limits)) {
  prediction_limits <- range(c(0, map_summary$min_kg_C_m2, map_summary$max_kg_C_m2), finite = TRUE)
}
if (is.null(sd_limits)) {
  sd_limits <- range(c(0, map_summary$max_pointwise_SD_kg_C_m2), finite = TRUE)
}
if (diff(prediction_limits) == 0) prediction_limits <- prediction_limits + c(0, 1)
if (diff(sd_limits) == 0) sd_limits <- sd_limits + c(0, 1)

draw_comparison <- function(maps, depth, quantity, limits) {
  columns <- min(panel_columns, nlyr(maps))
  previous_settings <- par(mfrow = c(ceiling(nlyr(maps) / columns), columns))
  on.exit(par(previous_settings))
  depth_label <- if (depth == "Total") "Total sampled profile" else {
    paste0("0-", sub("^0_", "", depth), " cm")
  }
  for (method in names(maps)) {
    label <- switch(method, area_based = "Area-based", idw = "IDW",
                    ok = "Ordinary kriging", Prior = "Prior", method)
    label <- gsub("_", " ", label)
    selected <- map_summary$selected_regression[
      map_summary$depth == depth & map_summary$method == method
    ]
    if (isTRUE(selected)) label <- paste0(label, " (selected)")
    title <- paste0(label, "\n", depth_label, " | ", quantity, " (kg C/m²)")
    if (any(is.finite(values(maps[[method]])))) {
      # Force a continuous scale even for uniform area-based/mean-only maps.
      plot(maps[[method]], type = "continuous", col = map_colours,
           range = limits, main = title, cex.main = 0.9,
           mar = c(3, 3, 4, 5))
    } else {
      plot.new()
      title(main = title, cex.main = 0.9)
      text(0.5, 0.5, "No saved / computable values")
    }
  }
}

for (depth in names(comparisons)) {
  for (quantity in c("prediction", "sd")) {
    maps <- comparisons[[depth]][[quantity]]
    limits <- if (quantity == "prediction") prediction_limits else sd_limits
    label <- if (quantity == "prediction") "Mean" else "Pointwise SD"
    filename <- if (quantity == "prediction") paste0("maps_", depth, ".png") else {
      paste0("uncertainty_", depth, ".png")
    }
    columns <- min(panel_columns, nlyr(maps))
    png(file.path(output_folder, filename), width = 750 * columns,
        height = 750 * ceiling(nlyr(maps) / columns), res = 150)
    draw_comparison(maps, depth, label, limits)
    dev.off()
    draw_comparison(maps, depth, label, limits)
  }
}
cat("\nFigures, raster stacks, map_summary.csv and map_methods.csv:", output_folder, "\n")

