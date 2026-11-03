# Four Li update branches, matched to the four workshop local methods.
# Explicit adaptation: original-scale Gaussian conditioning/fusion retains
# true zero organic depths/stocks. Negative posterior means are flagged, not clipped.
library(terra)
library(sf)
library(gstat)
library(readr)
library(dplyr)
source("config.R")
source("helpers/li_local_models.R")
source("helpers/workshop_gaussian.R")
source("helpers/extract_prior_at_cores.R")
out <- "outputs/script5/li_suite"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
local_dir <- "outputs/script3d_li_suite"
all_field <- read_csv(file.path(local_dir, "field_observations.csv"), show_col_types = FALSE)
field <- all_field[all_field$canfit, ]
stopifnot(nrow(field) >= 2, !anyDuplicated(field$core_id), all(field$depth >= 0), all(field$carbon >= 0))
template <- rast("data/boundary/analysis_grid.tif")
pts <- project(vect(field, geom = c("longitude", "latitude"), crs = "EPSG:4326"), crs(template))
xy <- crds(pts)
covs <- read_csv(file.path(local_dir, "native_core_covariates_all_119.csv"), show_col_types = FALSE)
covs <- covs[match(field$core_id, covs$core_id), ]
methods <- c("raw_cores", "area_based", "idw", "regression")
summaries <- list(); manifest <- list()

export_result <- function(result, cells, folder, template) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  stopifnot(nrow(result) == length(cells), all(is.finite(result$posterior_mean)),
    all(result$posterior_sd >= 0), all(result$variance_reduction >= -1e-8 & result$variance_reduction <= 1+1e-8))
  table <- cbind(data.frame(cell = cells, xyFromCell(template, cells)), result)
  write_csv(table, file.path(folder, "predictions.csv"))
  r <- rast(template, nlyrs = ncol(result)); names(r) <- names(result)
  v <- matrix(NA_real_, ncell(r), ncol(result)); v[cells,] <- as.matrix(result); values(r) <- v
  writeRaster(r, file.path(folder, "posterior_maps.tif"), overwrite = TRUE)
  invisible(r)
}

for (variable in c("depth", "carbon")) {
  folder <- file.path(out, variable)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  response <- field[[variable]]
  sf_points <- st_as_sf(pts); sf_points$response <- response
  model_data <- data.frame(core_id = field$core_id, response = response,
    covs[setdiff(names(covs), "core_id")], check.names = FALSE)
  native <- c(rast(paste0("data/priors/li_",variable,"_mean.tif")),
    rast(paste0("data/priors/li_",variable,"_sd.tif")))
  names(native) <- c("prior_mean", "prior_sd")
  native_at_core <- extract(native, project(pts, crs(native)), ID = FALSE)
  at_core <- fill_missing_core_priors(native, pts, native_at_core, field$core_id,
    file.path(folder, "core_prior_extraction_audit.csv"), allow_zero_mean = TRUE)
  write_csv(bind_cols(field, at_core), file.path(folder, "core_prior_audit.csv"))
  prior <- mask(project(native, template, method = "near"), template)
  v <- values(prior, mat = TRUE)
  cells <- which(is.finite(v[,1]) & v[,1] >= 0 & is.finite(v[,2]) & v[,2] > 0)
  stopifnot(length(cells) > 0)
  grid_xy <- xyFromCell(template, cells)
  local_maps <- rast(file.path(local_dir, variable, "local_maps.tif"))
  local_values <- values(local_maps, mat = TRUE)[cells,,drop=FALSE]
  idw_settings <- read_csv(file.path(local_dir, variable, "idw_parameters.csv"), show_col_types = FALSE)
  reg_settings <- read_csv(file.path(local_dir, variable, "selected_model.csv"), show_col_types = FALSE)
  scales <- c(area_based = sd(response), idw = idw_settings$RMSE[1], regression = reg_settings$RMSE[1])
  stopifnot(all(is.finite(scales) & scales > 0))
  obs_sd <- if (variable == "depth") fs_config$li_depth_observation_sd_cm else 0
  results <- list()
  results$raw_cores <- workshop_gaussian_gp(xy, response, at_core$prior_mean, at_core$prior_sd,
    grid_xy, v[cells,1], v[cells,2], fs_config$range_m, observation_sd = obs_sd)
  results$raw_cores$update_applied <- 1
  results$raw_cores$local_negative_or_missing <- 0
  for (method in methods[-1]) {
    # Workshop Script 5 uses physical local estimates; an unsupported negative
    # regression prediction is not silently converted to zero or used as evidence.
    local_mean <- local_values[,method]
    supported <- is.finite(local_mean) & local_mean >= 0
    z <- gaussian_results(v[cells,1], v[cells,2], v[cells,1], v[cells,2]^2)
    z$local_estimate <- local_mean
    z$data_weight <- 0
    if (any(supported)) z[supported,] <- workshop_gaussian_fusion(
      v[cells[supported],1], v[cells[supported],2], local_mean[supported], scales[[method]])
    z$update_applied <- as.integer(supported)
    z$local_negative_or_missing <- as.integer(!supported)
    results[[method]] <- z
  }

  # Outer holdouts refit local models, re-tune IDW and reselect regression.
  # Each working-error estimate uses only the outer training cores.
  cv <- list(); selection <- list()
  for (i in seq_len(nrow(field))) {
    tr <- setdiff(seq_len(nrow(field)), i)
    train <- model_data[tr,,drop=FALSE]
    idw <- li_idw_models(train, sf_points[tr,], sf_points[i,])
    reg <- li_regression_models(train, fs_config$covariate_models)
    local_prediction <- c(area_based = mean(response[tr]), idw = idw$mapped[1],
      regression = as.numeric(predict(reg$fit, model_data[i,,drop=FALSE])))
    error <- c(area_based = sd(response[tr]), idw = idw$parameters$RMSE[1], regression = reg$selection$RMSE[1])
    prediction <- list(raw_cores = workshop_gaussian_gp(xy[tr,,drop=FALSE], response[tr],
      at_core$prior_mean[tr], at_core$prior_sd[tr], xy[i,,drop=FALSE],
      at_core$prior_mean[i], at_core$prior_sd[i], fs_config$range_m, observation_sd = obs_sd))
    for (method in methods[-1]) {
      supported <- is.finite(local_prediction[[method]]) && local_prediction[[method]] >= 0
      prediction[[method]] <- if (supported) workshop_gaussian_fusion(at_core$prior_mean[i],
        at_core$prior_sd[i], local_prediction[[method]], error[[method]]) else
        gaussian_results(at_core$prior_mean[i], at_core$prior_sd[i], at_core$prior_mean[i], at_core$prior_sd[i]^2)
    }
    for (method in methods) {
      p <- prediction[[method]]
      cv[[paste(i,method)]] <- data.frame(method = method, core_id = field$core_id[i],
        observed = response[i], prior = p$prior_mean, predicted = p$posterior_mean,
        lower95 = p$lower95, upper95 = p$upper95, error = p$posterior_mean-response[i],
        local_prediction = if (method == "raw_cores") NA_real_ else local_prediction[[method]],
        local_supported = if (method == "raw_cores") TRUE else local_prediction[[method]] >= 0,
        validation_type = "Whole-core outer holdout; IDW tuning, regression selection and error scales refit within training data")
    }
    selection[[i]] <- data.frame(core_id = field$core_id[i], regression_model = reg$selected,
      idw_power = idw$parameters$idp[1], idw_neighbours = idw$parameters$nmax[1],
      area_error_sd = error[["area_based"]], idw_error_sd = error[["idw"]], regression_error_sd = error[["regression"]])
  }
  cv <- bind_rows(cv)
  write_csv(cv, file.path(folder, "nested_leave_one_core_out.csv"))
  write_csv(bind_rows(selection), file.path(folder, "nested_fold_choices.csv"))

  for (method in methods) {
    result <- results[[method]]
    export_result(result, cells, file.path(folder,method), template)
    cvi <- cv[cv$method == method,]
    write_csv(cvi, file.path(folder,method,"leave_one_core_out.csv"))
    summary <- data.frame(variable = variable, method = method,
      units = if (variable == "depth") "cm" else "kg C/m2",
      n_cores = nrow(field), confirmed_contacts = sum(field$confirmed),
      layer_estimated_contacts = sum(!field$confirmed), zero_observations = sum(response == 0),
      prior_cells = length(cells), update_applied_cells = sum(result$update_applied),
      prior_retained_cells = sum(!result$update_applied), negative_posterior_mean_cells = sum(result$posterior_mean < 0),
      mean_prior = mean(result$prior_mean), mean_updated = mean(result$posterior_mean),
      prior_RMSE = sqrt(mean((cvi$prior-cvi$observed)^2)),
      nested_update_RMSE = sqrt(mean(cvi$error^2)), nested_update_MAE = mean(abs(cvi$error)),
      interval_coverage = mean(cvi$observed >= cvi$lower95 & cvi$observed <= cvi$upper95),
      working_error_sd = if (method == "raw_cores") obs_sd else scales[[method]],
      response_model = "Original-scale Gaussian sensitivity; ensemble spreads are not calibrated errors",
      selected_regression = if (method == "regression") reg_settings$model[1] else NA_character_)
    write_csv(summary, file.path(folder,method,"summary.csv"))
    summaries[[paste(variable,method)]] <- summary
    manifest[[paste(variable,method)]] <- data.frame(variable = variable, method = method,
      local_evidence = if (method == "raw_cores") file.path(local_dir,variable,"raw_cores.csv") else
        file.path(local_dir,variable,"local_maps.tif"),
      update = file.path(folder,method,"posterior_maps.tif"), completed = TRUE,
      scale = "Original-scale Gaussian; true zeros retained")
  }
  # Range/observation-spread sensitivity for the raw spatial branch.
  scenarios <- expand.grid(range_m = c(1000,2000,4000),
    observation_sd = if (variable == "depth") c(2.5,5,10) else 0)
  sensitivity <- bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
    z <- workshop_gaussian_gp(xy,response,at_core$prior_mean,at_core$prior_sd,
      grid_xy,v[cells,1],v[cells,2],scenarios$range_m[i],scenarios$observation_sd[i])
    data.frame(scenarios[i,], mean_updated = mean(z$posterior_mean),
      negative_mean_cells = sum(z$posterior_mean < 0), mean_sd = mean(z$posterior_sd))
  }))
  write_csv(sensitivity, file.path(folder,"raw_cores/sensitivity.csv"))
}

# Carry the local carbon-depth curve through each depth update only on measured
# support. This diagnostic is separate from the direct organic-stock updates.
curve <- readRDS(file.path(local_dir,"carbon_curve/carbon_depth_curve.rds"))
max_depth <- max(field$depth)
for (method in methods) {
  depth <- rast(file.path(out,"depth",method,"posterior_maps.tif"))
  supported <- depth$posterior_mean >= 0 & depth$posterior_mean <= max_depth
  stock <- ifel(supported, depth$posterior_mean*unname(coef(curve)[1]), NA)
  names(stock) <- "local_curve_carbon_at_updated_depth"
  names(supported) <- "within_measured_curve_support"
  writeRaster(c(stock,supported),file.path(out,"depth",method,"curve_derived_carbon.tif"),overwrite=TRUE)
}
write_csv(bind_rows(summaries), file.path(out,"suite_summary.csv"))
write_csv(bind_rows(manifest), file.path(out,"completion_manifest.csv"))
writeLines(c("Four branches were fitted and updated for each Li variable: raw cores, area average, IDW and regression.",
  "All eight FS inputs retain confirmed-versus-layer-classified contact labels; four zero depths/stocks are retained.",
  "Gaussian conditioning/fusion is applied on the original scale to handle true zeros; this differs from workshop lognormal stocks.",
  "Regression predictions below zero or missing are flagged and the prior is retained in those cells.",
  "Negative Gaussian posterior means are also flagged and left visible; no hidden clipping is performed.",
  "FS-02/FS-04 need audited paired-source-pixel prior extraction within 500 m. Valid zero-valued depth priors are preserved.",
  "Depth and carbon remain dependent marginal sensitivity analyses; they are not multiplied into a joint posterior.",
  "The fitted local carbon-depth curve is mapped from every updated depth only within measured depth support.",
  "Outer holdouts refit local models, IDW tuning, regression selection and working error scales using training cores only.",
  "Classification error, prior calibration, spatial range and whole-inventory uncertainty remain outside these conditional intervals."),
  file.path(out,"interpretation.txt"))
message("Completed all eight Li prior updates (four methods x two variables).")
