# Four workshop local-method branches for Li-matched organic depth and carbon.
# Raw observations are points, not an invented fourth interpolated local surface.
library(terra)
library(sf)
library(gstat)
library(readr)
library(ggplot2)
source("config.R")
source("helpers/li_local_models.R")
out <- "outputs/script3d_li_suite"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
mode <- if (!is.null(fs_config$li_contact_mode)) fs_config$li_contact_mode else "all_eight_layer_evidence"
layers <- read.csv("data/field/core_layers.csv")
locations <- read.csv("data/field/core_locations.csv")
contacts <- read.csv("data/field/organic_mineral_contacts.csv")
threshold <- if (!is.null(fs_config$li_organic_soc_threshold_g_kg)) fs_config$li_organic_soc_threshold_g_kg else 80
all_field <- li_prepare_field_observations(layers, locations, contacts, mode,
  organic_soc_threshold_g_kg = threshold)
write_csv(all_field, file.path(out, "field_observations.csv"))
field <- all_field[all_field$canfit, , drop = FALSE]
if (nrow(field) < 2) stop("At least two eligible organic contacts are required")
covariates <- rast("data/covariates/all_predictors_1km.tif")
template <- rast("data/boundary/analysis_grid.tif")
stopifnot(compareGeom(covariates[[1]], template, stopOnError = FALSE))
native_cov <- read.csv("data/covariates/native_core_covariates.csv", check.names = FALSE)
stopifnot(!anyDuplicated(native_cov$core_id), setequal(native_cov$core_id, all_field$core_id))
native_cov <- native_cov[match(all_field$core_id, native_cov$core_id), , drop = FALSE]
stopifnot(identical(as.character(native_cov$core_id), as.character(all_field$core_id)))
write_csv(native_cov, file.path(out, "native_core_covariates_all_119.csv"))
coverage <- data.frame(predictor = setdiff(names(native_cov), "core_id"), stringsAsFactors = FALSE)
coverage$missing_all_cores <- vapply(native_cov[coverage$predictor], function(x) sum(!is.finite(x)), integer(1))
eligible_cov <- native_cov[match(field$core_id, native_cov$core_id), , drop = FALSE]
coverage$missing_eligible_cores <- vapply(eligible_cov[coverage$predictor], function(x) sum(!is.finite(x)), integer(1))
coverage$distinct_eligible_values <- vapply(eligible_cov[coverage$predictor],
  function(x) length(unique(x[is.finite(x)])), integer(1))
write_csv(coverage, file.path(out, "covariate_coverage.csv"))

# Scientific data checks: known measured contacts and zero labels are retained exactly.
stopifnot(!anyDuplicated(all_field$core_id), all(all_field$depth >= 0), all(all_field$carbon >= 0),
  all(all_field$carbon[all_field$depth == 0] == 0),
  all(all_field$depth[all_field$confirmed] == contacts$organic_horizon_depth_cm[
    match(all_field$core_id[all_field$confirmed], contacts$plot_id)]))
if (mode == "all_eight_layer_evidence") {
  stopifnot(nrow(field) == 8, sum(field$confirmed) == 2, sum(field$depth == 0) == 4,
    abs(field$depth[field$core_id == "FS-05"]-10) < 1e-10,
    abs(field$depth[field$core_id == "PM-02-A"]-15.2) < 1e-10)
}
cells <- which(is.finite(values(template, mat = FALSE)))
grid_points <- st_as_sf(data.frame(xyFromCell(template, cells)), coords = c("x", "y"), crs = crs(template))
suite <- list()
for (variable in c("depth", "carbon")) {
  folder <- file.path(out, variable)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  response <- field[[variable]]
  data <- data.frame(core_id = field$core_id, response = response,
    eligible_cov[setdiff(names(eligible_cov), "core_id")], check.names = FALSE)
  points <- st_transform(st_as_sf(data.frame(field, response = response),
    coords = c("longitude", "latitude"), crs = 4326), st_crs(crs(template)))
  write_csv(data.frame(field, response = response), file.path(folder, "raw_cores.csv"))
  st_write(points, file.path(folder, "raw_core_points.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  area_mean <- mean(response)
  area_cv <- vapply(seq_along(response), function(i) mean(response[-i]), numeric(1))
  area_map <- template * area_mean
  idw <- li_idw_models(data, points, grid_points)
  write_csv(idw$parameters, file.path(folder, "idw_parameters.csv"))
  write_csv(idw$tuning, file.path(folder, "idw_tuning.csv"))
  idw_map <- rast(template); values(idw_map) <- NA_real_
  values(idw_map)[cells] <- idw$mapped
  regression <- li_regression_models(data, fs_config$covariate_models)
  write_csv(regression$availability, file.path(folder, "candidate_availability.csv"))
  write_csv(regression$evaluation, file.path(folder, "model_evaluation.csv"))
  write_csv(regression$selection, file.path(folder, "selected_model.csv"))
  saveRDS(regression$fit, file.path(folder, "selected_model.rds"))
  saveRDS(regression$fits, file.path(folder, "fitted_models.rds"))
  regression_predictors <- all.vars(formula(regression$fit))[-1]
  regression_map <- if (!length(regression_predictors)) template * unname(coef(regression$fit)[1]) else {
    terra::predict(covariates[[regression_predictors]], regression$fit,
      fun = stats::predict, na.rm = TRUE)
  }
  maps <- mask(c(area_map, idw_map, regression_map), template)
  names(maps) <- c("area_based", "idw", "regression")
  writeRaster(maps, file.path(folder, "local_maps.tif"), overwrite = TRUE)
  flags <- maps < 0
  names(flags) <- paste0(names(maps), "_negative_prediction")
  outside <- template * 0
  for (p in regression_predictors) {
    outside <- outside | (covariates[[p]] < min(data[[p]]) | covariates[[p]] > max(data[[p]]))
  }
  names(outside) <- "regression_outside_training_range"
  writeRaster(c(flags, outside), file.path(folder, "prediction_flags.tif"), overwrite = TRUE)
  pred <- list(raw_cores = rep(NA_real_, nrow(field)), area_based = area_cv,
    idw = idw$heldout, regression = regression$heldout)
  predictions <- do.call(rbind, lapply(names(pred), function(method) data.frame(
    method = method, core_id = field$core_id, observed = response,
    predicted = pred[[method]], error = pred[[method]]-response,
    validation_type = if (method == "raw_cores") "Observed point data; no local prediction/validation" else
      if (method == "idw") idw$status else if (method == "regression")
        "Whole-core holdout for selected fixed model; AICc selection is not nested" else
        "Whole-core holdout area mean", stringsAsFactors = FALSE)))
  write_csv(predictions, file.path(folder, "local_model_predictions.csv"))
  all_cv <- do.call(rbind, lapply(names(regression$heldouts), function(model) data.frame(
    model = model, core_id = field$core_id, observed = response,
    predicted = regression$heldouts[[model]], error = regression$heldouts[[model]]-response)))
  write_csv(all_cv, file.path(folder, "regression_leave_one_core_out.csv"))
  summary <- do.call(rbind, lapply(names(pred), function(method) {
    score <- li_cv_metrics(response, pred[[method]])
    data.frame(variable = variable, method = method, n_cores = nrow(field),
      n_confirmed_contacts = sum(field$confirmed), n_layer_estimated_contacts = sum(!field$confirmed),
      observed_mean = mean(response), between_core_sd = sd(response),
      LOOCV_RMSE = score["RMSE"], LOOCV_MAE = score["MAE"], LOOCV_bias = score["bias"],
      status = if (nrow(field) < 4) "Descriptive small-n fallback; not a validated covariate model" else
        "Exploratory field/classification evidence; model selection/tuning not nested", row.names = NULL)
  }))
  write_csv(summary, file.path(folder, "model_summary.csv"))
  map_summary <- do.call(rbind, lapply(names(maps), function(method) {
    v <- values(maps[[method]], mat = FALSE); valid <- v[is.finite(v)]
    data.frame(method = method, valid_cells = length(valid), map_mean = mean(valid),
      minimum = min(valid), maximum = max(valid), negative_cells = sum(valid < 0))
  }))
  write_csv(map_summary, file.path(folder, "map_summary.csv"))
  stopifnot(all(is.finite(values(maps$area_based)[cells])),
    max(abs(values(maps$area_based)[cells]-mean(response))) < 1e-10,
    min(values(maps$idw)[cells], na.rm = TRUE) >= min(response)-1e-8,
    max(values(maps$idw)[cells], na.rm = TRUE) <= max(response)+1e-8)
  suite[[variable]] <- summary
}

# Cumulative carbon curve: zero-depth cores remain visible but do not add slope
# information or residual degrees of freedom to the positive-profile model.
positive <- field[field$depth > 0, , drop = FALSE]
if (nrow(positive) < 2) stop("Two positive complete organic profiles are needed for a descriptive curve")
curve <- lm(carbon ~ 0 + depth, data = positive)
curve_cv <- vapply(seq_len(nrow(positive)), function(i) as.numeric(predict(
  lm(carbon ~ 0 + depth, data = positive[-i, , drop = FALSE]), positive[i, , drop = FALSE])), numeric(1))
curve_dir <- file.path(out, "carbon_curve")
dir.create(curve_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(curve, file.path(curve_dir, "carbon_depth_curve.rds"))
write_csv(positive, file.path(curve_dir, "positive_profile_endpoints.csv"))
write_csv(data.frame(core_id = positive$core_id, observed = positive$carbon, predicted = curve_cv,
  error = curve_cv-positive$carbon), file.path(curve_dir, "leave_one_positive_core_out.csv"))
depth_seq <- seq(0, max(positive$depth), length.out = 101)
write_csv(data.frame(depth_cm = depth_seq, carbon_kg_m2 = as.numeric(coef(curve)[1])*depth_seq),
  file.path(curve_dir, "supported_cumulative_curve.csv"))
knots <- do.call(rbind, lapply(seq_len(nrow(field)), function(i) {
  id <- field$core_id[i]; limit <- field$depth[i]
  d <- layers[layers$plot_id == id & layers$depth_from < limit, , drop = FALSE]
  d <- d[order(d$depth_from), , drop = FALSE]
  segment_stock <- d$soc*d$bulk_density*(1-d$coarse_frag)*
    pmax(0, pmin(d$depth_to, limit)-d$depth_from)/100
  data.frame(core_id = id, depth_cm = c(0, pmin(d$depth_to, limit)),
    cumulative_stock_kg_m2 = c(0, cumsum(segment_stock)), confirmed = field$confirmed[i])
}))
write_csv(knots, file.path(curve_dir, "measured_organic_cumulative_curves.csv"))
curve_summary <- data.frame(n_cores_total = nrow(field), n_positive_profiles_fitted = nrow(positive),
  n_zero_profiles_displayed_only = sum(field$depth == 0),
  slope_kg_C_m2_per_cm = unname(coef(curve)[1]), residual_df = df.residual(curve),
  measured_depth_support_max_cm = max(positive$depth),
  holdout_RMSE_kg_C_m2 = sqrt(mean((curve_cv-positive$carbon)^2)),
  model = "Through-zero cumulative-carbon curve; positive-profile endpoints only",
  status = "Exploratory: two confirmed contacts plus any authorized layer-derived contacts")
write_csv(curve_summary, file.path(curve_dir, "model_summary.csv"))
p <- ggplot(knots, aes(depth_cm, cumulative_stock_kg_m2, colour = core_id)) +
  geom_line(data = knots[knots$depth_cm > 0 | knots$core_id %in% positive$core_id, ]) +
  geom_point(aes(shape = confirmed), size = 2.5) +
  geom_abline(intercept = 0, slope = coef(curve)[1], linetype = 2, colour = "black") +
  theme_minimal(base_size = 12) + labs(title = "Fort Severn cumulative organic carbon with depth",
    subtitle = paste(nrow(positive), "positive profiles fit; zero-depth cores shown at the origin"),
    x = "Organic depth (cm)", y = "Cumulative organic carbon (kg C/m²)", shape = "Confirmed contact",
    caption = "Dashed: pooled through-zero curve. Classification-derived contacts remain labelled; no extrapolation beyond measured support.")
ggsave(file.path(curve_dir, "cumulative_carbon_curves.png"), p, width = 9, height = 5.5, dpi = 180)
write_csv(do.call(rbind, suite), file.path(out, "all_local_model_summary.csv"))
writeLines(c(paste("Contact mode:", mode),
  "All field contact labels are reconstructed from raw measured layers and explicit contacts; no Route L model outputs are sourced.",
  "Four local method branches: raw measured points; area mean; workshop-tuned IDW; workshop AICc-selected regression.",
  "Zeros remain exact zero labels; negative regression extrapolations are retained and flagged, not silently clamped.",
  "Area/IDW/regression local maps do not themselves constitute Li prior updates.",
  "Confirmed versus layer-classified evidence remains distinct throughout.",
  "Current classification labels are exploratory; tuning and model selection do not provide independent-site validation."),
  file.path(out, "interpretation.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
message("Li local four-method suite completed in ", out)
