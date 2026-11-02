# ==============================================

# ------- Spatial Mapping Continued ---------

# ===========================================

# Here, we use the same the core data and teh same AOI
#Import the data -

# Now we are going to add environmental variables to give more information
# to the spacial spatial models. These include:
#===============================================================
# Part 3 - Raster-based carbon models
#===============================================================
# Run each numbered section in order.
# Inputs: Script1 core stocks, the same AOI, and the downloaded covariate TIFFs.
# This script retains the stock definitions supplied by Script1.
# Area-average / IDW comparisons and carbon inventory totals belong in Script5.

# Install missing packages once:
# install.packages(c("terra", "sf", "readr", "dplyr", "ggplot2", "jsonlite"))

library(terra)
library(sf)
library(readr)
library(dplyr)
library(ggplot2)
library(jsonlite)


# ==============================================================
# Step 3A - Settings and core data
# ==============================================================

# EDIT THIS: the unzipped folder containing topography.tif, vegetation.tif, etc.
covariate_folder <- "Draft1/data/covariates"

# These are used if the Script2 spatial objects are not already loaded.
mapping_file <- "Draft1/data/mapping_data.csv"
aoi_file <- "Draft1/data/aoi.geojson"

# Start with 0-30 cm. Later change to stock_0_50, stock_0_100 or stock_Total.
stock_variable <- "stock_0_30"
output_root <- "outputs/script3"
output_folder <- file.path(output_root, stock_variable)

# Full raster checks read the rasters in chunks and may take time on large AOIs.
full_raster_checks <- TRUE
make_maps <- TRUE
show_plots <- TRUE

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

if (exists("core_proj") && exists("aoi_proj")) {
  message("Using core_proj and aoi_proj from Script2.")
  analysis_cores <- core_proj
  analysis_aoi <- aoi_proj
} else {
  message("Reading the Script1 mapping table and AOI.")
  cores <- read_csv(mapping_file, show_col_types = FALSE) %>%
    rename(
      core_id = `Core ID`,
      latitude = Latitude,
      longitude = Longitude,
      stock_0_30 = `C Stock 0-30cm`,
      stock_0_50 = `C Stock 0-50cm`,
      stock_0_100 = `C Stock 0-100cm`,
      stock_Total = `Total Carbon Stock`
    )
  analysis_cores <- st_as_sf(
    cores, coords = c("longitude", "latitude"), crs = 4326
  ) %>% st_transform(3978)
  analysis_aoi <- st_read(aoi_file, quiet = TRUE) %>% st_transform(3978)
}

if (!all(c("core_id", stock_variable) %in% names(analysis_cores))) {
  stop("Core data need core_id and the selected stock column: ", stock_variable)
}
if (anyNA(analysis_cores$core_id) || anyDuplicated(analysis_cores$core_id)) {
  stop("Each core needs one unique, non-missing core_id. Check the Core Log.")
}
if (any(!st_is_valid(analysis_aoi))) stop("Check invalid AOI geometry before continuing.")

print(st_drop_geometry(analysis_cores)[, c("core_id", stock_variable)])


# ==============================================================
# Step 3B - Read and check the raster groups
# ==============================================================

raster_files <- file.path(covariate_folder, c(
  "topography.tif"
))
group_names <- c(
  "topography"
)
names(raster_files) <- group_names

missing_files <- raster_files[!file.exists(raster_files)]
if (length(missing_files)) {
  stop("Check covariate_folder. Missing files:\n", paste(missing_files, collapse = "\n"))
}

band_manifest_file <- file.path(covariate_folder, "band_manifest.csv")
band_manifest <- if (file.exists(band_manifest_file)) {
  read_csv(band_manifest_file, show_col_types = FALSE)
} else NULL

settings_file <- file.path(covariate_folder, "settings.json")
if (file.exists(settings_file)) {
  extraction_settings <- fromJSON(settings_file)
  print(extraction_settings[intersect(
    c("analysis_scale", "export_cell_m", "context_radii_m", "start_year", "end_year", "radar_pass"),
    names(extraction_settings)
  )])
  file.copy(settings_file, file.path(output_folder, "extraction_settings.json"), overwrite = TRUE)
}

raster_groups <- lapply(raster_files, rast)

# Recover the exact band names from the export manifest when it is available.
for (group in group_names) {
  r <- raster_groups[[group]]
  if (!is.null(band_manifest)) {
    bands <- band_manifest %>%
      filter(.data$group == .env$group) %>% arrange(band_number)
    if (nrow(bands) != nlyr(r)) stop("Band count does not match the manifest: ", group)
    names(r) <- bands$band_name
  }
  raster_groups[[group]] <- r
  cat("\n", group, ":\n", paste(names(r), collapse = ", "), "\n", sep = "")
}

template <- rast(raster_groups$topography[[1]])
if (is.lonlat(template)) stop("Use the projected metre-based raster exports from Colab.")
for (group in group_names) {
  if (!compareGeom(template, raster_groups[[group]], stopOnError = FALSE)) {
    stop("Raster grids do not match: ", group, ". Check the exports before resampling.")
  }
}

all_band_names <- unlist(lapply(raster_groups, names), use.names = FALSE)
if (anyDuplicated(all_band_names)) stop("Raster band names must be unique across groups.")

analysis_cores <- st_transform(analysis_cores, crs(template))
analysis_aoi <- st_sf(geometry = st_union(st_transform(analysis_aoi, crs(template))))
core_vector <- vect(analysis_cores)
aoi_vector <- vect(analysis_aoi)

inside <- lengths(st_intersects(analysis_cores, analysis_aoi)) > 0
if (any(!inside)) {
  stop("These cores are outside the AOI: ", paste(analysis_cores$core_id[!inside], collapse = ", "))
}

# Use the same centre-based AOI mask for all rasters and predictions.
aoi_mask <- rasterize(aoi_vector, template, field = 1, background = NA, touches = FALSE)
aoi_cells <- global(aoi_mask, "notNA")[1, 1]
if (aoi_cells == 0) stop("The AOI and exported raster extent do not overlap.")

raster_inventory <- do.call(rbind, lapply(group_names, function(group) {
  r <- raster_groups[[group]]
  data.frame(
    group = group, file = raster_files[[group]], bands = nlyr(r),
    rows = nrow(r), columns = ncol(r), cell_width_m = res(r)[1], cell_height_m = res(r)[2]
  )
}))
print(raster_inventory)
write_csv(raster_inventory, file.path(output_folder, "raster_inventory.csv"))

if (full_raster_checks) {
  raster_checks <- do.call(rbind, lapply(group_names, function(group) {
    r <- mask(raster_groups[[group]], aoi_mask)
    check <- global(r, c("min", "max", "notNA"), na.rm = TRUE)
    check$min[!is.finite(check$min)] <- NA_real_
    check$max[!is.finite(check$max)] <- NA_real_
    data.frame(
      group = group, band = names(r), minimum = check$min, maximum = check$max,
      valid_AOI_cells = check$notNA, missing_AOI_pct = 100 * (1 - check$notNA / aoi_cells)
    )
  }))
  print(raster_checks)
  write_csv(raster_checks, file.path(output_folder, "raster_checks.csv"))
}

if (show_plots) {
  plot(raster_groups$topography[[1]], main = "Elevation and core locations")
  plot(aoi_vector, add = TRUE, border = "black")
  points(core_vector, pch = 21, bg = "red")
}


# ==============================================================
# Step 3C - Extract the raster values at each core
# ==============================================================

covariates <- rast(unname(raster_files))
names(covariates) <- all_band_names
core_values <- terra::extract(covariates, core_vector, method = "simple", ID = FALSE)

model_data <- data.frame(
  core_id = as.character(analysis_cores$core_id),
  stock = analysis_cores[[stock_variable]],
  core_values,
  check.names = FALSE
)

# Class IDs are categories; the default models below use probabilities instead.
if ("dw_class_mode" %in% names(model_data)) {
  model_data$dw_class_mode <- factor(model_data$dw_class_mode, levels = 0:8)
}

band_groups <- rep(group_names, times = vapply(raster_groups, nlyr, numeric(1)))
point_checks <- do.call(rbind, lapply(seq_along(all_band_names), function(i) {
  band <- all_band_names[i]
  x <- model_data[[band]]
  observed <- if (is.numeric(x)) x[is.finite(x)] else x[!is.na(x)]
  data.frame(
    group = band_groups[i], band = band,
    role = if (band_groups[i] == "quality") "QA" else "candidate_predictor",
    missing_cores = nrow(model_data) - length(observed),
    distinct_values_at_cores = length(unique(observed)),
    minimum = if (is.numeric(x) && length(observed)) min(observed) else NA_real_,
    maximum = if (is.numeric(x) && length(observed)) max(observed) else NA_real_
  )
}))

print(head(model_data))
print(point_checks)
write_csv(model_data, file.path(output_folder, "core_covariates.csv"))
write_csv(point_checks, file.path(output_folder, "core_covariate_checks.csv"))
saveRDS(model_data, file.path(output_folder, "core_covariates.rds"))

# Keep counts for checking data quality, not as default model predictors.
quality_columns <- names(raster_groups$quality)
print(model_data[, c("core_id", quality_columns), drop = FALSE])


# ==============================================================
# Step 3D - Define a small set of raster-based models
# ==============================================================

# EDIT THESE FORMULAS after inspecting the band names and core values above.
# Use one predictor per group first, followed by a small combined model.
# stock is the selected Script1 response; its depth bins are unchanged.
model_formulas <- list(
  Topography = stock ~ relative_elevation_120m,
  Vegetation = stock ~ ndmi_spring_median,
  Land_cover = stock ~ dw_tree_prob_mean,
  Radar = stock ~ vh_db_median,
  Water = stock ~ water_occurrence_pct,
  Climate = stock ~ mean_temperature_c,
  Combined = stock ~ relative_elevation_120m + ndmi_spring_median
)

selected_predictors <- unique(unlist(lapply(model_formulas, function(f) all.vars(f)[-1])))
available_predictors <- intersect(selected_predictors, names(model_data))
numeric_predictors <- available_predictors[vapply(
  model_data[available_predictors], is.numeric, logical(1)
)]
varying_predictors <- numeric_predictors[vapply(
  model_data[numeric_predictors], function(x) length(unique(x[is.finite(x)])) > 1,
  logical(1)
)]

# Exploratory correlations between the prespecified predictors.
if (length(varying_predictors) > 1) {
  correlations <- suppressWarnings(cor(
    model_data[varying_predictors], use = "pairwise.complete.obs"
  ))
  print(round(correlations, 2))
  write_csv(
    data.frame(predictor = rownames(correlations), correlations, check.names = FALSE),
    file.path(output_folder, "selected_predictor_correlations.csv")
  )
}

for (predictor in varying_predictors) {
  p <- ggplot(model_data, aes(x = .data[[predictor]], y = stock)) +
    geom_point(size = 3, na.rm = TRUE) +
    geom_text(aes(label = core_id), vjust = -0.7, size = 3, na.rm = TRUE) +
    labs(x = predictor, y = "Carbon stock (kg C/m²)", title = stock_variable) +
    theme_minimal()
  if (show_plots) print(p)
  ggsave(file.path(output_folder, paste0("relationship_", predictor, ".png")),
         p, width = 7, height = 5, dpi = 180)
}


# ==============================================================
# Step 3E - Two small helpers for fitting and validation
# ==============================================================

fit_linear_model <- function(formula, data) {
  design <- model.matrix(formula, data = data)
  if (nrow(design) <= ncol(design)) {
    stop("Too few cores for this number of model coefficients.")
  }
  if (qr(design)$rank < ncol(design)) {
    stop("Predictors are constant or linearly dependent in these training cores.")
  }
  lm(formula, data = data, na.action = na.fail)
}

leave_one_core_out <- function(formula, data) {
  result <- data.frame(
    core_id = data$core_id, observed = data$stock,
    predicted = NA_real_, error = NA_real_, message = ""
  )
  for (i in seq_len(nrow(data))) {
    training <- data[-i, , drop = FALSE]
    test <- data[i, , drop = FALSE]
    prediction <- tryCatch({
      model <- fit_linear_model(formula, training)
      as.numeric(predict(model, newdata = test))
    }, error = function(e) {
      result$message[i] <<- conditionMessage(e)
      NA_real_
    })
    result$predicted[i] <- prediction
    result$error[i] <- prediction - result$observed[i]
  }
  result
}


# ==============================================================
# Step 3F - Fit and assess each model separately
# ==============================================================

fitted_models <- list()
training_tables <- list()
model_notes <- list()

for (model_name in names(model_formulas)) {
  cat("\n--- ", model_name, " ---\n", sep = "")
  formula <- model_formulas[[model_name]]
  predictors <- all.vars(formula)[-1]
  required_columns <- c("stock", predictors)
  model_folder <- file.path(output_folder, model_name)
  dir.create(model_folder, showWarnings = FALSE)
  
  missing_bands <- setdiff(predictors, names(model_data))
  reason <- ""
  if (length(missing_bands)) {
    reason <- paste("Missing bands:", paste(missing_bands, collapse = ", "))
  } else if (any(predictors %in% quality_columns)) {
    reason <- "QA counts are excluded from these environmental models."
  } else if (!all(vapply(model_data[predictors], is.numeric, logical(1)))) {
    reason <- "This first version fits numeric predictors/probabilities; class IDs stay categorical."
  }
  if (nzchar(reason)) {
    message("Skipping ", model_name, ": ", reason)
    model_notes[[model_name]] <- data.frame(model = model_name, status = "Skipped", reason = reason)
    writeLines(reason, file.path(model_folder, "not_fitted.txt"))
    next
  }
  
  # Different models may have different available cores. Each records its own IDs.
  keep <- Reduce(`&`, lapply(model_data[required_columns], is.finite))
  training <- model_data[keep, c("core_id", required_columns), drop = FALSE]
  excluded <- model_data[!keep, c("core_id", required_columns), drop = FALSE]
  write_csv(excluded, file.path(model_folder, "excluded_cores.csv"))
  
  # Leave one residual degree of freedom even in a withheld-core training fit.
  if (nrow(training) < length(predictors) + 3) {
    reason <- "Too few complete cores for this model and leave-one-core-out validation."
    fit <- NULL
  } else {
    fit <- tryCatch(fit_linear_model(formula, training), error = function(e) {
      reason <<- conditionMessage(e)
      NULL
    })
  }
  if (is.null(fit)) {
    message("Skipping ", model_name, ": ", reason)
    model_notes[[model_name]] <- data.frame(model = model_name, status = "Skipped", reason = reason)
    writeLines(reason, file.path(model_folder, "not_fitted.txt"))
    next
  }
  
  print(summary(fit))
  validation <- leave_one_core_out(formula, training)
  complete_validation <- all(is.finite(validation$predicted))
  fit_statistics <- data.frame(
    model = model_name, target = stock_variable, n_cores = nrow(training),
    formula = paste(deparse(formula), collapse = " "),
    fitted_R2 = summary(fit)$r.squared,
    LOOCV_predictions = sum(is.finite(validation$predicted)),
    LOOCV_RMSE = if (complete_validation) sqrt(mean(validation$error^2)) else NA_real_,
    LOOCV_MAE = if (complete_validation) mean(abs(validation$error)) else NA_real_,
    validation_status = if (complete_validation) "Complete" else "Incomplete: inspect fold messages"
  )
  print(fit_statistics)
  
  coefficients <- data.frame(term = rownames(coef(summary(fit))), coef(summary(fit)), check.names = FALSE)
  fitted_values <- data.frame(
    core_id = training$core_id, observed = training$stock,
    fitted = as.numeric(fitted(fit)), residual = as.numeric(residuals(fit))
  )
  write_csv(training, file.path(model_folder, "training_data.csv"))
  write_csv(coefficients, file.path(model_folder, "coefficients.csv"))
  write_csv(fitted_values, file.path(model_folder, "fitted_values.csv"))
  write_csv(validation, file.path(model_folder, "leave_one_core_out.csv"))
  write_csv(fit_statistics, file.path(model_folder, "model_statistics.csv"))
  saveRDS(fit, file.path(model_folder, "model.rds"))
  
  diagnostic_data <- rbind(
    data.frame(core_id = training$core_id, observed = training$stock,
               prediction = as.numeric(fitted(fit)), type = "Fitted"),
    data.frame(core_id = validation$core_id, observed = validation$observed,
               prediction = validation$predicted, type = "Leave one core out")
  )
  diagnostic_plot <- ggplot(diagnostic_data, aes(observed, prediction)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) +
    geom_point(size = 3, na.rm = TRUE) +
    geom_text(aes(label = core_id), vjust = -0.7, size = 3, na.rm = TRUE) +
    facet_wrap(~ type) + theme_minimal() +
    labs(title = paste(model_name, stock_variable),
         x = "Observed stock (kg C/m²)", y = "Predicted stock (kg C/m²)")
  if (show_plots) print(diagnostic_plot)
  ggsave(file.path(model_folder, "model_diagnostics.png"), diagnostic_plot,
         width = 10, height = 5, dpi = 180)
  
  residual_plot <- ggplot(fitted_values, aes(fitted, residual)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
    geom_point(size = 3) + theme_minimal() +
    labs(title = model_name, x = "Fitted stock (kg C/m²)", y = "Residual (observed - fitted)")
  if (show_plots) print(residual_plot)
  ggsave(file.path(model_folder, "residuals.png"), residual_plot,
         width = 7, height = 5, dpi = 180)
  
  fitted_models[[model_name]] <- fit
  training_tables[[model_name]] <- training
  model_notes[[model_name]] <- data.frame(
    model = model_name, status = "Fitted",
    reason = if (complete_validation) "" else "Some withheld-core fits failed; inspect validation."
  )
}

# This is a processing log, not a ranked model comparison.
model_log <- do.call(rbind, model_notes)
write_csv(model_log, file.path(output_folder, "model_run_log.csv"))
if (!length(fitted_models)) stop("No models were fitted. Inspect the core checks and model run log.")


# ==============================================================
# Step 3G - Predict each fitted model across the AOI
# ==============================================================

if (make_maps) {
  for (model_name in names(fitted_models)) {
    cat("\nMapping ", model_name, "\n", sep = "")
    fit <- fitted_models[[model_name]]
    training <- training_tables[[model_name]]
    predictors <- all.vars(formula(fit))[-1]
    model_folder <- file.path(output_folder, model_name)
    
    # Supply only the bands used by this model; unrelated NA bands do not mask it.
    model_rasters <- mask(covariates[[predictors]], aoi_mask)
    prediction <- terra::predict(
      model_rasters, fit, fun = stats::predict, na.rm = TRUE,
      filename = file.path(model_folder, "carbon_stock_kg_m2.tif"),
      overwrite = TRUE, wopt = list(gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
    )
    
    # Flag values beyond the observed training range of any model predictor.
    outside_range <- model_rasters[[1]] < min(training[[predictors[1]]]) |
      model_rasters[[1]] > max(training[[predictors[1]]])
    if (length(predictors) > 1) {
      for (predictor in predictors[-1]) {
        outside_range <- outside_range |
          model_rasters[[predictor]] < min(training[[predictor]]) |
          model_rasters[[predictor]] > max(training[[predictor]])
      }
    }
    flags <- c(outside_range, prediction < 0)
    names(flags) <- c("outside_training_range", "negative_prediction")
    writeRaster(flags, file.path(model_folder, "prediction_flags.tif"),
                overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
    
    map_stats <- global(prediction, c("min", "max", "notNA"), na.rm = TRUE)
    flag_counts <- global(flags, "sum", na.rm = TRUE)
    map_summary <- data.frame(
      model = model_name, target = stock_variable,
      validation_status = if (nzchar(model_notes[[model_name]]$reason)) "Incomplete" else "Complete",
      minimum_kg_m2 = map_stats$min, maximum_kg_m2 = map_stats$max,
      valid_AOI_cells = map_stats$notNA,
      missing_AOI_pct = 100 * (1 - map_stats$notNA / aoi_cells),
      outside_training_range_cells = flag_counts$sum[1],
      negative_prediction_cells = flag_counts$sum[2]
    )
    print(map_summary)
    write_csv(map_summary, file.path(model_folder, "map_checks.csv"))
    
    # Raw predictions are retained; flags are exported separately for review.
    png(file.path(model_folder, "carbon_stock_map.png"), width = 1500, height = 1200, res = 150)
    plot(prediction, main = paste(model_name, stock_variable, "(kg C/m²)"), maxcell = 100000)
    plot(aoi_vector, add = TRUE, border = "black")
    points(core_vector, pch = 21, bg = "white")
    dev.off()
    if (show_plots) {
      plot(prediction, main = paste(model_name, stock_variable, "(kg C/m²)"), maxcell = 100000)
      plot(aoi_vector, add = TRUE, border = "black")
      points(core_vector, pch = 21, bg = "white")
    }
  }
}


# ==============================================================
# Step 3H - Save the analysis setup for later use
# ==============================================================

model_specifications <- data.frame(
  model = names(model_formulas), target = stock_variable,
  formula = vapply(model_formulas, function(f) paste(deparse(f), collapse = " "), character(1))
)
write_csv(model_specifications, file.path(output_folder, "model_specifications.csv"))
saveRDS(fitted_models, file.path(output_folder, "fitted_models.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_folder, "R_session_info.txt"))
cat("\nScript3 outputs saved to:", normalizePath(output_folder), "\n")
# Each model folder contains its own diagnostics, validation, fitted model and maps.
# Models can use different available cores; Script5 must align cases for comparisons.
# Synthetic core responses test the workflow, not real environmental relationships.
# A different stock_variable is a separate run; independently fitted depth maps may cross.
