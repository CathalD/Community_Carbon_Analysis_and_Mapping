# Part 2.3: Bayesian carbon updates
# Each method updates the same prior separately; they reuse the same cores and
# must not be applied sequentially as independent evidence.
library(terra)
library(readr)
library(dplyr)
source("config.R")
source("helpers/grid_helper.R")

# 1. Inputs and settings
target_depth_cm <- bayesian_depth_cm
stock_column <- paste0("C Stock 0-", target_depth_cm, "cm")
range_m <- 2000
observation_log_sd <- 0.00  # Demonstration assumption; set from field/lab error.
chunk_size <- 5000L
map_methods <- c("area_based", "idw")
output_folder <- file.path("outputs/script5", paste0("stock_0_", target_depth_cm))
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# 2. Read the shared grid and prior already prepared by Script 4
shared <- load_analysis_grid()
template <- shared$template
analysis_mask <- shared$mask
prior <- rast(file.path("outputs/script4", paste0("prior_moments_0_", target_depth_cm, "cm.tif")))
compareGeom(prior, template)
prior <- mask(prior, analysis_mask)
prior_values <- values(prior, mat = TRUE)
eligible_cells <- which(apply(is.finite(prior_values), 1, all) &
                          prior_values[, "prior_mean"] > 0 &
                          prior_values[, "log_variance"] > 0)
grid_xy <- xyFromCell(template, eligible_cells)
grid_mean <- prior_values[eligible_cells, "prior_mean"]
grid_sd <- prior_values[eligible_cells, "prior_sd"]
grid_log_mean <- prior_values[eligible_cells, "log_mean"]
grid_log_var <- prior_values[eligible_cells, "log_variance"]

raw_cores <- read_csv(mapping_file, show_col_types = FALSE)
cores <- data.frame(
  core_id = raw_cores[["Core ID"]], plot_id = raw_cores[["Plot ID"]],
  longitude = raw_cores$Longitude, latitude = raw_cores$Latitude,
  stock = raw_cores[[stock_column]]
)
cores <- cores[is.finite(cores$stock) & cores$stock > 0, ]
core_points <- project(vect(cores, geom = c("longitude", "latitude"), crs = coordinate_crs),
                       crs(template))
core_xy <- crds(core_points)

# Lookup the source prior at the actual core coordinates, not at moved centres.
# The final predictions still use the shared analysis grid.
depth_key <- as.character(target_depth_cm)
native_prior <- c(rast(prior_mean_files[[depth_key]]),
                  rast(prior_width_files[[depth_key]]) / (2 * qnorm(0.95)))
names(native_prior) <- c("prior_mean", "prior_sd")
at_cores <- extract(native_prior, project(core_points, crs(native_prior)), ID = FALSE)
core_log_var <- log1p((at_cores$prior_sd / at_cores$prior_mean)^2)
core_log_mean <- log(at_cores$prior_mean) - core_log_var / 2
cores$prior_mean <- at_cores$prior_mean
cores$prior_sd <- at_cores$prior_sd
cat("\nTarget:", stock_column, "Cores:", nrow(cores),
    "Plots:", length(unique(cores$plot_id)), "\n")
cat("Depth stocks retain Script 1's completed-profile-at-sampled-maximum assumption.\n")
cat("Raw-core range:", range_m, "m; observation log SD:", observation_log_sd, "\n")
print(cores)

# Common export of log-Gaussian posterior moments in carbon units.
# Intervals are pointwise and conditional on fixed model/error parameters.
save_posterior <- function(log_mean, log_variance, method) {
  folder <- file.path(output_folder, method)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  tiny_negative <- log_variance < 0 & log_variance > -1e-8
  log_variance[tiny_negative] <- 0
  if (any(log_variance < 0, na.rm = TRUE)) {
    warning(method, ": negative posterior log variance; affected cells exported as NA.")
    log_variance[log_variance < 0] <- NA_real_
  }
  results <- data.frame(
    cell = eligible_cells, x = grid_xy[, 1], y = grid_xy[, 2],
    prior_mean = grid_mean, prior_sd = grid_sd,
    posterior_log_mean = log_mean, posterior_log_var = log_variance,
    posterior_mean = exp(log_mean + log_variance / 2),
    posterior_median = exp(log_mean),
    posterior_variance = expm1(log_variance) * exp(2 * log_mean + log_variance),
    posterior_sd = sqrt(expm1(log_variance) * exp(2 * log_mean + log_variance)),
    lower95 = exp(log_mean - qnorm(0.975) * sqrt(log_variance)),
    upper95 = exp(log_mean + qnorm(0.975) * sqrt(log_variance))
  )
  results$change_from_prior <- results$posterior_mean - results$prior_mean
  results$log_variance_reduction <- 1 - log_variance / grid_log_var
  fields <- c("posterior_mean", "posterior_variance", "posterior_log_mean", "posterior_log_var")
  maps <- rast(template, nlyrs = length(fields))
  names(maps) <- fields
  cells <- matrix(NA_real_, ncell(template), length(fields))
  cells[eligible_cells, ] <- as.matrix(results[, fields])
  values(maps) <- cells
  maps <- mask(maps, analysis_mask)
  writeRaster(maps[[1]], file.path(folder, "posterior_mean.tif"), overwrite = TRUE)
  writeRaster(maps[[2]], file.path(folder, "posterior_variance.tif"), overwrite = TRUE)
  writeRaster(maps[[3:4]], file.path(folder, "posterior_log_moments.tif"),
              datatype = "FLT8S", overwrite = TRUE)
  write_csv(results, file.path(folder, "predictions.csv"))
  titles <- c(paste(method, "posterior mean (kg C/m²)"), "Posterior variance ((kg C/m²)²)")
  png(file.path(folder, "posterior_maps.png"), width = 1600, height = 800, res = 150)
  plot(maps[[1:2]], main = titles)
  dev.off()
  plot(maps[[1:2]], main = titles)
  cat("\n", method, "posterior summary:\n")
  print(data.frame(
    method = method, cells = sum(is.finite(results$posterior_mean)),
    mean_prior = mean(results$prior_mean),
    mean_posterior = mean(results$posterior_mean, na.rm = TRUE),
    mean_pointwise_sd = mean(results$posterior_sd, na.rm = TRUE)
  ))
  results
}

# 3. Raw-core spatial update: retain the existing exponential covariance model
distance_matrix <- function(a, b) {
  sqrt(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2)
}
core_distances <- distance_matrix(core_xy, core_xy)
core_log_sd <- sqrt(core_log_var)
K <- outer(core_log_sd, core_log_sd) * exp(-core_distances / range_m)
observation_var <- rep(observation_log_sd^2, nrow(cores))
V <- K + diag(observation_var, nrow = nrow(cores))
residual_log <- log(cores$stock) - core_log_mean
R <- chol(V)
alpha <- backsolve(R, forwardsolve(t(R), residual_log))
core_posterior_log_mean <- numeric(length(eligible_cells))
core_posterior_log_var <- numeric(length(eligible_cells))
for (start in seq(1, length(eligible_cells), by = chunk_size)) {
  rows <- start:min(start + chunk_size - 1, length(eligible_cells))
  K_new <- outer(sqrt(grid_log_var[rows]), sqrt(core_log_var)) *
    exp(-distance_matrix(grid_xy[rows, , drop = FALSE], core_xy) / range_m)
  core_posterior_log_mean[rows] <- grid_log_mean[rows] + as.vector(K_new %*% alpha)
  projected <- forwardsolve(t(R), t(K_new))
  core_posterior_log_var[rows] <- grid_log_var[rows] - colSums(projected^2)
}
core_results <- save_posterior(core_posterior_log_mean, core_posterior_log_var, "raw_cores")
saveRDS(list(
  range_m = range_m, observation_log_sd = observation_log_sd,
  target_depth_cm = target_depth_cm, cores = cores, core_xy = core_xy,
  posterior_log_mean = core_posterior_log_mean,
  posterior_log_var = core_posterior_log_var
), file.path(output_folder, "raw_cores", "model.rds"))

# 4. Map updates: read calibrated log locations and error variances from Script 2
# Local errors and prior uncertainty are assumed independent for this working
# pointwise update. Cross-pixel error covariance is not estimated here.
for (method in map_methods) {
  local_file <- file.path("outputs/script2", paste0(method, "_log_moments_0_", target_depth_cm, "cm.tif"))
  local <- rast(local_file)
  compareGeom(local, template)
  local_values <- values(local, mat = TRUE)[eligible_cells, , drop = FALSE]
  local_log_mean <- local_values[, "log_mean"]
  local_error_var <- local_values[, "log_variance"]
  weight <- grid_log_var / (grid_log_var + local_error_var)
  post_mean <- grid_log_mean + weight * (local_log_mean - grid_log_mean)
  post_var <- (1 - weight) * grid_log_var
  cat("\n", method, "input error calibration:\n")
  calibration_file <- file.path("outputs/script2", paste0(method, "_calibration_0_", target_depth_cm, "cm.csv"))
  if (file.exists(calibration_file)) print(read_csv(calibration_file, show_col_types = FALSE))
  results <- save_posterior(post_mean, post_var, method)
}
cat("\nBayesian rasters:", output_folder, "\n")
cat("Reported intervals are pointwise; no AOI-total uncertainty is calculated.\n")
