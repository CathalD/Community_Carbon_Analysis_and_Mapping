# ===================================
# Part 5: Updating existing maps with local field data
# ===================================

library(terra)
library(readr)
library(dplyr)
source("config.R")
source("helpers/grid_helper.R")

# ----------------
# Step 1 - Set up the analysis and define the inputs
# ----------------

bayesian_depth_cm <- 100 # A matching prior and Script 2 maps must exist.
target_depth_cm <- bayesian_depth_cm
stock_column <- paste0("C Stock 0-", target_depth_cm, "cm")
map_methods <- c("area_based", "idw")
output_folder <- file.path("outputs/script5", paste0("stock_0_", target_depth_cm))
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Retained for the later raw-core spatial-conditioning example.
# Pixel-wise conjugate updating below does not use a geographic range.
range_m <- 2000
observation_log_sd <- 0 # Assumed raw-core measurement error, not map error.
chunk_size <- 5000L

# A. Prediction locations: the centres of the shared grid cells inside the AOI.
shared <- load_analysis_grid()
template <- shared$template
analysis_mask <- shared$mask

# B. Prior: the existing McMaster estimate BEFORE these local data are used.
# Script 4 copied paired stock/uncertainty values to cell centres using nearest
# source values. These are point predictions, not computed cell/block averages.
# It approximated SD = (Q95-Q05)/(2*qnorm(0.95)), then moment-matched a lognormal:
# v0 = log(1 + SD^2/mean^2); m0 = log(mean) - v0/2.
# The interval-to-SD conversion is an assumption, not a measured CV error.
prior_file <- file.path("outputs/script4", paste0("prior_moments_0_", target_depth_cm, "cm.tif"))
if (!file.exists(prior_file)) stop("Run Script 4 for the selected depth first: ", prior_file)
prior <- rast(prior_file)
compareGeom(prior, template)
stopifnot(all(c("prior_mean", "prior_sd", "log_mean", "log_variance") %in% names(prior)))
prior <- mask(prior, analysis_mask)
prior_values <- values(prior, mat = TRUE)
eligible_cells <- which(apply(is.finite(prior_values), 1, all) &
                          prior_values[, "prior_mean"] > 0 &
                          prior_values[, "prior_sd"] > 0 &
                          prior_values[, "log_variance"] > 0)
if (!length(eligible_cells)) stop("No cells have valid positive prior moments.")
grid_xy <- xyFromCell(template, eligible_cells)
grid_mean <- prior_values[eligible_cells, "prior_mean"]
grid_sd <- prior_values[eligible_cells, "prior_sd"]
grid_log_mean <- prior_values[eligible_cells, "log_mean"]
grid_log_var <- prior_values[eligible_cells, "log_variance"]

prediction_table <- data.frame(
  cell_id = eligible_cells, x = grid_xy[, 1], y = grid_xy[, 2],
  prior_mean = grid_mean, prior_sd = grid_sd,
  prior_log_mean = grid_log_mean, prior_log_variance = grid_log_var
)
predicted_points <- vect(prediction_table, geom = c("x", "y"),
                         crs = crs(template), keepgeom = TRUE)
# predicted_points is the map's prediction grid; core_points below are observations.

# C. Observations: Script 1 carbon stocks, unchanged.
# All supplied cores retain Script 1's capped-at-recorded-maximum convention.
raw_cores <- read_csv(mapping_file, show_col_types = FALSE)
if (!stock_column %in% names(raw_cores)) stop("Missing Script 1 stock column: ", stock_column)
cores <- data.frame(
  core_id = raw_cores[["Core ID"]], plot_id = raw_cores[["Plot ID"]],
  longitude = raw_cores$Longitude, latitude = raw_cores$Latitude,
  stock = raw_cores[[stock_column]]
)
cores <- cores[is.finite(cores$stock) & cores$stock > 0 &
                 is.finite(cores$longitude) & is.finite(cores$latitude), ]
if (nrow(cores) < 2L) stop("At least two positive, geolocated core stocks are needed.")
core_points <- project(vect(cores, geom = c("longitude", "latitude"), crs = coordinate_crs),
                       crs(template))
core_xy <- crds(core_points)

# Retain the native prior at actual observed coordinates for later raw-core work.
# This lookup does not move cores to grid centres or add them to the map update.
depth_key <- as.character(target_depth_cm)
if (!depth_key %in% intersect(names(prior_mean_files), names(prior_width_files))) {
  stop("No configured prior stock/uncertainty pair for depth ", target_depth_cm)
}
native_mean <- rast(prior_mean_files[[depth_key]])
native_width <- rast(prior_width_files[[depth_key]])
compareGeom(native_mean, native_width)
native_prior <- c(native_mean, native_width / (2 * qnorm(0.95)))
names(native_prior) <- c("prior_mean", "prior_sd")
at_cores <- extract(native_prior, project(core_points, crs(native_prior)), ID = FALSE)
core_log_var <- log1p((at_cores$prior_sd / at_cores$prior_mean)^2)
core_log_mean <- log(at_cores$prior_mean) - core_log_var / 2
cores$prior_mean <- at_cores$prior_mean
cores$prior_sd <- at_cores$prior_sd

# D. Local maps: estimates made from these same field cores in Script 2.
# Area-based = equal weighting of log stock, constant over the area.
# IDW = distance weighting of log stock, varying by prediction location.
# Read their saved log moments directly; do not log a carbon-scale SD.
local_maps <- list()
for (method in map_methods) {
  local_file <- file.path("outputs/script2", paste0(method, "_log_moments_0_", target_depth_cm, "cm.tif"))
  if (!file.exists(local_file)) stop("Run Script 2 first: ", local_file)
  local <- rast(local_file)
  compareGeom(local, template)
  if (!all(c("log_mean", "log_variance") %in% names(local))) {
    stop("Local map must contain log_mean and log_variance: ", local_file)
  }
  local_maps[[method]] <- mask(local, analysis_mask)
}
area_based <- local_maps[["area_based"]]
area_based_values <- values(area_based, mat = TRUE)[eligible_cells, , drop = FALSE]
area_based_points <- prediction_table[, c("cell_id", "x", "y")]
area_based_points$local_log_mean <- area_based_values[, "log_mean"]
area_based_points$local_log_variance <- area_based_values[, "log_variance"]
area_based_points$local_mean <- exp(area_based_points$local_log_mean +
                                    area_based_points$local_log_variance / 2)
area_based_points$local_sd <- sqrt(expm1(area_based_points$local_log_variance) *
  exp(2 * area_based_points$local_log_mean + area_based_points$local_log_variance))

# The raw-data mean/SD below describe the observed sample in carbon units.
# The mapped area estimate is a lognormal approximation and can differ from them.
field_summary <- data.frame(
  depth_cm = target_depth_cm, cores = nrow(cores), plots = length(unique(cores$plot_id)),
  arithmetic_stock_mean = mean(cores$stock), observed_stock_sd = sd(cores$stock),
  mean_log_stock = mean(log(cores$stock)), observed_log_sd = sd(log(cores$stock)),
  area_mean_log_variance = var(log(cores$stock)) / nrow(cores),
  area_new_location_log_variance = var(log(cores$stock)) * (1 + 1 / nrow(cores))
)
# s^2/n describes uncertainty in an IID area mean; s^2*(1+1/n) describes
# prediction at an unobserved location. Script 2 uses the latter working model.
# IDW uses its held-out whole-plot log-error variance. Zero measurement error
# does not set either of these map prediction-error variances to zero.
cat("\nTarget:", stock_column, "| Prediction cells:", length(eligible_cells), "\n")
print(field_summary)
print(cores)

# --------------------
# Step 2 - Normal-Normal conjugate updating
# --------------------

# Question: how should each cell's prior change given a local map estimate?
# Bayes' theorem: p(Y | local data) = p(local data | Y) * p(Y) / p(local data).
# Here Y = log(C), with C in kg C/m2. We assume:
#   Prior:      Y ~ Normal(m0, v0)
#   Likelihood: yhat | Y ~ Normal(Y, vL)
# m0/v0 are the prior log mean/variance; yhat/vL are the local estimate/error variance.
# v0 and vL are VARIANCES, not SDs. Smaller variance means larger precision (1/v).
#
# Posterior variance: v_post = 1 / (1/v0 + 1/vL)
# Posterior log mean: m_post = v_post * (m0/v0 + yhat/vL)
# Equivalently: w_local = v0/(v0+vL)
#              m_post = m0 + w_local*(yhat-m0)
#              v_post = v0*vL/(v0+vL)
#
# These are exact for an independent, unbiased likelihood for the same target.
# Using interpolation error maps as vL is an APPROXIMATE fusion assumption.
# Matching grids alone does not establish compatible depth/support, unbiased
# estimates, independent training data, or calibrated prediction-error variance.
# Area/IDW reuse the same cores: update the original prior separately with each;
# do not feed one posterior into the other or add the raw cores again.

normal_normal_update <- function(m0, v0, yhat, vL) {
  valid <- is.finite(m0) & is.finite(v0) & v0 > 0 &
    is.finite(yhat) & is.finite(vL) & vL >= 0
  w_local <- m_post <- v_post <- rep(NA_real_, length(m0))
  w_local[valid] <- v0[valid] / (v0[valid] + vL[valid])
  m_post[valid] <- m0[valid] + w_local[valid] * (yhat[valid] - m0[valid])
  v_post[valid] <- v0[valid] * (vL[valid] / (v0[valid] + vL[valid]))
  data.frame(posterior_log_mean = m_post, posterior_log_variance = v_post,
             local_weight = w_local)
}

# --------------------
# Step 3 - Apply the formula to every eligible prediction cell
# --------------------

posterior_maps <- list()
posterior_points <- list()
for (method in map_methods) {
  local_values <- values(local_maps[[method]], mat = TRUE)[eligible_cells, , drop = FALSE]
  yhat <- local_values[, "log_mean"]
  vL <- local_values[, "log_variance"]
  updated <- normal_normal_update(grid_log_mean, grid_log_var, yhat, vL)
  if (!any(is.finite(updated$posterior_log_mean))) stop("No valid local moments for ", method)
  m_post <- updated$posterior_log_mean
  v_post <- updated$posterior_log_variance

  # For a lognormal posterior: median = exp(m); mean = exp(m+v/2).
  # The most probable carbon-scale value (mode) is exp(m-v), a third quantity.
  # Variance = (exp(v)-1)*exp(2*m+v). Intervals below are pointwise and conditional
  # on the working model and fixed error estimates; they are not calibrated CV coverage.
  results <- cbind(prediction_table, data.frame(
    method = method, local_log_mean = yhat, local_log_variance = vL
  ), updated)
  results$posterior_mean <- exp(m_post + v_post / 2)
  results$posterior_median <- exp(m_post)
  results$posterior_mode <- exp(m_post - v_post)
  results$posterior_variance <- expm1(v_post) * exp(2 * m_post + v_post)
  results$posterior_sd <- sqrt(results$posterior_variance)
  results$lower95 <- exp(m_post - qnorm(0.975) * sqrt(v_post))
  results$upper95 <- exp(m_post + qnorm(0.975) * sqrt(v_post))

  fields <- c("posterior_mean", "posterior_sd", "posterior_variance",
              "posterior_log_mean", "posterior_log_variance",
              "posterior_median", "posterior_mode", "lower95", "upper95", "local_weight")
  maps <- rast(template, nlyrs = length(fields))
  names(maps) <- fields
  cells <- matrix(NA_real_, ncell(template), length(fields))
  cells[eligible_cells, ] <- as.matrix(results[, fields])
  values(maps) <- cells
  maps <- mask(maps, analysis_mask)
  posterior_maps[[method]] <- maps
  posterior_points[[method]] <- results

  folder <- file.path(output_folder, method)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  write_csv(results, file.path(folder, "predictions.csv"))
  writeRaster(maps, file.path(folder, "posterior_all_moments.tif"), datatype = "FLT8S", overwrite = TRUE)
  writeRaster(maps$posterior_mean, file.path(folder, "posterior_mean.tif"), overwrite = TRUE)
  writeRaster(maps$posterior_variance, file.path(folder, "posterior_variance.tif"), overwrite = TRUE)
  log_maps <- maps[[c("posterior_log_mean", "posterior_log_variance")]]
  names(log_maps) <- c("posterior_log_mean", "posterior_log_var")
  writeRaster(log_maps, file.path(folder, "posterior_log_moments.tif"), datatype = "FLT8S", overwrite = TRUE)
  png(file.path(folder, "posterior_maps.png"), width = 1500, height = 750, res = 150)
  plot(maps[[c("posterior_mean", "posterior_sd")]],
       main = c(paste(method, "posterior mean (kg C/m2)"), "Pointwise posterior SD (kg C/m2)"))
  dev.off()
}

# The requested posterior stock raster: current depth, area-based local estimate.
# At target_depth_cm=100 this is a 0-100 cm result; do not label it post_stock_30.
post_stock <- posterior_maps[["area_based"]]$posterior_mean
post_stock_sd <- posterior_maps[["area_based"]]$posterior_sd
post_stock_idw <- posterior_maps[["idw"]]$posterior_mean

# --------------------
# Step 4 - Visualize the distributions at one example cell
# --------------------

# Choose a cell whose prior mean is nearest the median valid prior mean.
example_candidates <- which(is.finite(area_based_values[, "log_mean"]) &
                              is.finite(area_based_values[, "log_variance"]) &
                              area_based_values[, "log_variance"] > 0)
if (length(example_candidates)) {
  example_row <- example_candidates[which.min(abs(grid_mean[example_candidates] -
                                                   median(grid_mean[example_candidates])))]
  example_cell <- eligible_cells[example_row]
  example <- posterior_points[["area_based"]][example_row, ]
  distribution_parameters <- data.frame(
    distribution = c("McMaster prior", "Area-based local estimate", "Updated posterior"),
    log_mean = c(example$prior_log_mean, example$local_log_mean, example$posterior_log_mean),
    log_variance = c(example$prior_log_variance, example$local_log_variance,
                     example$posterior_log_variance)
  )
  distribution_parameters$mean <- exp(distribution_parameters$log_mean +
                                        distribution_parameters$log_variance / 2)
  distribution_parameters$sd <- sqrt(expm1(distribution_parameters$log_variance) *
    exp(2 * distribution_parameters$log_mean + distribution_parameters$log_variance))
  colours <- c("#3268A8", "#D27B26", "#31865A")
  draw_distributions <- function() {
    previous <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 4, 1))
    on.exit(par(previous))
    hist(cores$stock, breaks = "Sturges", col = "grey85", border = "white",
         xlab = "Observed stock (kg C/m2)", main = paste0("Field data: 0-", target_depth_cm, " cm"))
    abline(v = mean(cores$stock), col = colours[2], lwd = 2)
    legend("topright", legend = c(sprintf("Sample mean = %.2f", mean(cores$stock)),
                                  sprintf("Sample SD = %.2f", sd(cores$stock))), bty = "n", cex = 0.8)
    mu <- distribution_parameters$log_mean
    sigma <- sqrt(distribution_parameters$log_variance)
    log_x <- seq(min(mu - 4 * sigma), max(mu + 4 * sigma), length.out = 1000)
    log_y <- vapply(1:3, function(j) dnorm(log_x, mu[j], sigma[j]), numeric(length(log_x)))
    matplot(log_x, log_y, type = "l", lty = 1, lwd = 2, col = colours,
            xlab = "Log carbon stock", ylab = "Density", main = paste("Log-scale update | cell", example_cell))
    legend("topleft", legend = c("Prior", "Local likelihood (normalized)", "Posterior"),
           col = colours, lty = 1, lwd = 2, bty = "n", cex = 0.75)
    carbon_x <- seq(0, max(qlnorm(0.995, mu, sigma)), length.out = 1500)
    carbon_y <- vapply(1:3, function(j) dlnorm(carbon_x, mu[j], sigma[j]), numeric(length(carbon_x)))
    matplot(carbon_x, carbon_y, type = "l", lty = 1, lwd = 2, col = colours,
            xlab = "Carbon stock (kg C/m2)", ylab = "Density", main = "Implied stock distributions")
    # Local curve is its lognormal working uncertainty model. It is not a
    # likelihood density in C; the normalized likelihood above is in log C.
    legend("topright", legend = c("Prior", "Local estimate", "Posterior"), col = colours,
           lty = 1, lwd = 2, bty = "n", cex = 0.8)
  }
  png(file.path(output_folder, "conjugate_distribution_example.png"),
      width = 1950, height = 700, res = 150)
  draw_distributions()
  dev.off()
  if (interactive()) draw_distributions()
  write_csv(distribution_parameters, file.path(output_folder, "example_cell_distributions.csv"))
  cat("\nDistribution example - cell", example_cell, "\n")
  print(distribution_parameters)
}

write_csv(field_summary, file.path(output_folder, "field_summary.csv"))
write_csv(prediction_table, file.path(output_folder, "prior_prediction_points.csv"))
write_csv(area_based_points, file.path(output_folder, "area_based_prediction_points.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_folder, "session_info.txt"))
cat("\nResults:", output_folder, "\n")
cat("post_stock = area-based posterior mean; post_stock_idw = IDW posterior mean.\n")
cat("These are alternative, approximate pointwise updates of the same original prior.\n")
cat("Regional-total uncertainty requires covariance or joint simulation.\n")
