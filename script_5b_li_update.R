# Apply the workshop Script 5 spatial update separately to Li depth and organic carbon.
# These are marginal sensitivity analyses. Depth and carbon priors are dependent:
# do not multiply the two posteriors or add their uncertainty as independent pools.
library(terra)
library(readr)
library(dplyr)
source("config.R")
source("helpers/workshop_gp.R")
out <- "outputs/script5/li"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
local <- read_csv("outputs/script3b_li_local/li_field_observations.csv", show_col_types = FALSE)
template <- rast("data/boundary/analysis_grid.tif")
points <- project(vect(local, geom = c("longitude", "latitude"), crs = "EPSG:4326"), crs(template))
for (variable in c("depth", "carbon")) {
  folder <- file.path(out, variable)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  native <- c(rast(paste0("data/priors/li_", variable, "_mean.tif")),
    rast(paste0("data/priors/li_", variable, "_sd.tif")))
  names(native) <- c("prior_mean", "prior_sd")
  at_core <- extract(native, project(points, crs(native)), ID = FALSE)
  prior <- mask(project(native, template, method = "near"), template)
  pv <- values(prior, mat = TRUE)
  cells <- which(is.finite(pv[,1]) & pv[,1] > 0 & is.finite(pv[,2]) & pv[,2] > 0)
  observed <- if (variable == "depth") local$peat_depth_cm else local$organic_carbon_kg_m2
  eligible <- is.finite(observed) & observed > 0 & is.finite(at_core$prior_mean) &
    at_core$prior_mean > 0 & is.finite(at_core$prior_sd) & at_core$prior_sd > 0
  audit <- bind_cols(local, at_core) %>% mutate(usable_positive_prior = eligible)
  write_csv(audit, file.path(folder, "core_prior_audit.csv"))
  if (sum(eligible) < 2) stop("Less than two contacts have usable Li priors for ", variable)
  # Depth observation spread is an assumption, not a measured laboratory error.
  # Carbon uses the workshop's zero measurement-error baseline.
  obs_sd_log <- if (variable == "depth") sqrt(log1p((fs_config$li_depth_observation_sd_cm / observed[eligible])^2)) else 0
  args <- list(core_xy = crds(points)[eligible,,drop=FALSE], observed = observed[eligible],
    core_mean = at_core$prior_mean[eligible], core_sd = at_core$prior_sd[eligible],
    new_xy = xyFromCell(template, cells), new_mean = pv[cells,1], new_sd = pv[cells,2],
    range_m = fs_config$range_m, observation_log_sd = obs_sd_log)
  result <- do.call(workshop_gp, args)
  result$nearest_contact_m <- apply(distance_matrix(args$new_xy, args$core_xy), 1, min)
  result$cell <- cells
  write_csv(result, file.path(folder, "predictions.csv"))
  map_cols <- setdiff(names(result), "cell")
  r <- rast(template, nlyrs = length(map_cols)); names(r) <- map_cols
  v <- matrix(NA_real_, ncell(r), length(map_cols)); v[cells,] <- as.matrix(result[map_cols]); values(r) <- v
  writeRaster(r, file.path(folder, "posterior_maps.tif"), overwrite = TRUE)
  png(file.path(folder, "posterior_mean.png"), width = 1200, height = 1000, res = 150)
  plot(r$posterior_mean, main = paste("Li", variable, "— workshop update (2 contacts)"))
  points(points, pch = 21, bg = "white"); dev.off()
  # Leave a whole core out; there is no regional calibration or independent-site claim.
  idx <- which(eligible)
  cv <- bind_rows(lapply(seq_along(idx), function(i) {
    train <- setdiff(seq_along(idx), i)
    pred <- workshop_gp(args$core_xy[train,,drop=FALSE], args$observed[train],
      args$core_mean[train], args$core_sd[train], args$core_xy[i,,drop=FALSE],
      args$core_mean[i], args$core_sd[i], range_m = fs_config$range_m,
      observation_log_sd = rep(obs_sd_log, length.out = length(idx))[train])
    data.frame(core_id = local$core_id[idx[i]], observed = observed[idx[i]],
      prior = pred$prior_mean, predicted = pred$posterior_mean,
      lower95 = pred$lower95, upper95 = pred$upper95,
      error = pred$posterior_mean-observed[idx[i]])
  }))
  write_csv(cv, file.path(folder, "leave_one_core_out.csv"))
  # Show dependence on the assumed range and (depth only) observation spread.
  scenarios <- expand.grid(range_m = c(1000, 2000, 4000),
    observation_sd = if (variable == "depth") c(2.5, 5, 10) else 0)
  sensitivity <- bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
    a <- args; a$range_m <- scenarios$range_m[i]
    a$observation_log_sd <- if (variable == "depth") sqrt(log1p((scenarios$observation_sd[i]/a$observed)^2)) else 0
    z <- do.call(workshop_gp, a)
    data.frame(scenarios[i,], mean_updated = mean(z$posterior_mean), mean_sd = mean(z$posterior_sd))
  }))
  write_csv(sensitivity, file.path(folder, "sensitivity.csv"))
  write_csv(data.frame(variable = variable, units = if (variable == "depth") "cm" else "kg C/m2",
    n_confirmed_contacts = sum(eligible), prediction_cells = length(cells),
    mean_prior = mean(result$prior_mean), mean_updated = mean(result$posterior_mean),
    LOOCV_RMSE = sqrt(mean(cv$error^2)), status = "Exploratory marginal update; uncalibrated ensemble spread"),
    file.path(folder, "summary.csv"))
}

# Evaluate the fitted local carbon-depth curve only on its measured support.
# This is an alternative diagnostic, not a second independent update.
depth <- rast(file.path(out, "depth/posterior_maps.tif"))$posterior_mean
fit <- readRDS("outputs/script3b_li_local/carbon_depth_curve.rds")
max_depth <- max(local$peat_depth_cm)
curve_carbon <- ifel(depth <= max_depth, depth*unname(coef(fit)[1]), NA)
names(curve_carbon) <- "local_curve_carbon_within_measured_depth_support"
support <- ifel(!is.na(depth), depth <= max_depth, NA)
names(support) <- "within_measured_curve_support"
writeRaster(c(curve_carbon, support), file.path(out, "local_curve_at_updated_depth.tif"), overwrite = TRUE)
writeLines(c("Li depth and carbon updates use the same workshop exponential log-Gaussian method as Sothe.",
  "They are separate marginal comparisons; Li carbon is derived from depth and is not independent evidence.",
  "Do not combine these posteriors into a joint inventory uncertainty.",
  "Li SD rasters describe disagreement across base models, not calibrated predictive error.",
  "Only two confirmed peat contacts support this exploratory run. All eight field stocks remain in Sothe.",
  "Local cumulative-curve mapping is NA beyond 30.73 cm rather than silently inventing a deep peat tail.",
  "Peat depths and organic carbon do not include a separate mineral-soil carbon pool."),
  file.path(out, "interpretation.txt"))
