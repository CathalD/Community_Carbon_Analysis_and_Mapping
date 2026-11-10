# ==============================================
# Part 2 - Spatial interpolation
# ==============================================

library(sf)
library(terra)
library(gstat)
library(readr)
library(ggplot2)
source("helpers/grid_helper.R")

# 1. Define settings and organize the data ----------------------------

# Coordinate systems: geographic field locations, then a projected analysis grid.
coordinate_crs <- "EPSG:4326"
analysis_crs <- "EPSG:3978"
cell_size_m <- 1000                          # Side length of a grid cell, in metres.
mapping_depths_cm <- c(30, 50, 100)

# Input files produced in Script 1, and the study boundary supplied by the user.
mapping_file <- "outputs/script1/mapping_data.csv"
aoi_file <- "data/boundary/aoi.geojson"

# Shared grid defining analysis and output size
# These values also occur in config.R for the later scripts; keep them aligned.
grid_folder <- "data/grid"
grid_file <- file.path(grid_folder, "template_1km.tif")
mask_file <- file.path(grid_folder, "aoi_mask_1km.tif")
cell_id_file <- file.path(grid_folder, "cell_ids_1km.tif")
projected_aoi_file <- file.path(grid_folder, "aoi_projected.gpkg")
# The helper builds the grid from cell_size_m; the filenames above retain the existing labels.

output_folder <- "outputs/script2"
testing_folder <- file.path(output_folder, "test")
ranking_folder <- file.path(output_folder, "parameter_rankings")
dir.create(testing_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(ranking_folder, recursive = TRUE, showWarnings = FALSE)
map_colours <- hcl.colors(100, "Viridis")
method_labels <- c(area_based = "Area-based", idw = "IDW", ok = "Matern kriging")

# Initial settings for the worked 0-30 cm example.
idw_power <- 2
idw_neighbours <- Inf                         # All available cores.
matern_scale_m <- 20000                       # gstat distance scale in metres.
matern_kappa <- 0.5                           # nu: smoothness.
matern_nugget_fraction <- 0                   # Keep fixed in the parameter comparison.

# Settings to compare later. Edit these grids as needed.
idw_powers <- c(0, 0.5, 1, 2, 3, 4)
idw_neighbour_counts <- c(2, 3, 4, Inf)
matern_scales_m <- c(5000, 10000, 20000, 40000, 80000)
matern_kappas <- c(0.5, 1.5, 2.5, 10)


# ============

depth_variables <- setNames(paste0("C Stock 0-", mapping_depths_cm, "cm"),
                            paste0("0_", mapping_depths_cm))
depth_variables <- c(depth_variables, Total = "Total Carbon Stock")
depth_labels <- setNames(paste0("0-", mapping_depths_cm, " cm"), paste0("0_", mapping_depths_cm))
depth_labels <- c(depth_labels, Total = "Total sampled profile")
cores <- read_csv(mapping_file, show_col_types = FALSE)
core_proj <- st_transform(st_as_sf(cores, coords = c("Longitude", "Latitude"),
                                 crs = coordinate_crs), analysis_crs)
shared <- create_analysis_grid()
template <- shared$template
analysis_mask <- shared$mask
pred_points <- st_as_sf(as.points(shared$cell_ids))
plot_fold <- match(core_proj[["Plot ID"]], unique(core_proj[["Plot ID"]]))


#----------------------------------
# Prediction points at cell centres
plot(st_geometry(st_as_sf(shared$aoi)), main = "Cores and prediction points")
plot(st_geometry(pred_points), add = TRUE, pch = ".", col = "grey80")
plot(st_geometry(core_proj), add = TRUE, pch = 19, col = "red")


# 2. Work through the three methods for 0-30 cm -----------------------

test_cores <- core_proj
test_cores$stock <- test_cores[[depth_variables[["0_30"]]]]
test_cores$log_stock <- log(test_cores$stock)
print(st_drop_geometry(test_cores)[, c("Core ID", "stock", "log_stock")])

# Each method predicts log carbon as a weighted combination of observed log stocks.
# Area: equal weights. IDW: inverse-distance weights.
# Kriging: weights from covariance among cores AND their relation to the target.
# Matern nu=0.5 gives exponential covariance
# Uncertainty here concerns prediction at an unmeasured location.


# Lognormal mean = exp(m+v/2); variance = (exp(v)-1)*exp(2*m+v).
# These are working lognormal approximations, conditional on the stated models.
export_maps <- function(m, v, method, depth, folder) {
  v <- rep_len(v, nrow(pred_points))
  points <- pred_points
  points$prediction <- exp(m + v / 2)
  points$variance <- expm1(v) * exp(2 * m + v)
  points$sd <- sqrt(points$variance)
  points$log_mean <- m
  points$log_variance <- v
  maps <- mask(rasterize(vect(points), template,
                         field = c("prediction", "variance", "sd", "log_mean", "log_variance")),
               analysis_mask)
  suffix <- if (depth == "Total") "total" else paste0(depth, "cm")
  writeRaster(maps$prediction, file.path(folder, paste0(method, "_", suffix, ".tif")), overwrite = TRUE)
  writeRaster(maps$variance, file.path(folder, paste0(method, "_variance_", suffix, ".tif")), overwrite = TRUE)
  writeRaster(maps$sd, file.path(folder, paste0(method, "_sd_", suffix, ".tif")), overwrite = TRUE)
  writeRaster(maps[[c("log_mean", "log_variance")]],
              file.path(folder, paste0(method, "_log_moments_", suffix, ".tif")),
              datatype = "FLT8S", overwrite = TRUE)
  titles <- c(paste(method_labels[[method]], depth_labels[[depth]], "mean (kg C/m2)"), "Pointwise SD (kg C/m2)")
  png(file.path(folder, paste0(method, "_", suffix, ".png")), width = 1400, height = 650, res = 150)
  plot(maps[[c("prediction", "sd")]], type = "continuous", col = map_colours, main = titles)
  dev.off()
  plot(maps[[c("prediction", "sd")]], type = "continuous", col = map_colours, main = titles)
  maps
}

# 2A. Area-based: equal weights --------------------------------------

n <- nrow(test_cores)
area_m <- rep(mean(test_cores$log_stock), nrow(pred_points))
area_v <- var(test_cores$log_stock) * (1 + 1 / n)
# s^2*(1+1/n): new-location prediction variance under an IID constant-mean model.
# s^2/n would describe uncertainty in the area mean instead.
area_30 <- export_maps(area_m, area_v, "area_based", "0_30", testing_folder)

# Predict withheld plots using the remaining cores' mean.
area_heldout <- numeric(n)
area_heldout_variance <- numeric(n)
for (fold in unique(plot_fold)) {
  training_y <- test_cores$log_stock[plot_fold != fold]
  area_heldout[plot_fold == fold] <- mean(training_y)
  area_heldout_variance[plot_fold == fold] <- var(training_y) * (1 + 1 / length(training_y))
}
cv_area_30 <- data.frame(observed = test_cores$log_stock, var1.pred = area_heldout,
                         var1.var = area_heldout_variance)

print(cv_area_30)


# 2B. IDW: inverse-distance weights ----------------------------------

# w_i = d_i^(-power) / sum(d_j^(-power)).
# With power=0 and all cores, this is equal weighting away from coincident cores.
idw_prediction_30 <- idw(log_stock ~ 1, test_cores, pred_points,
                         idp = idw_power, nmax = idw_neighbours, debug.level = 0)
cv_idw_30 <- krige.cv(log_stock ~ 1, test_cores, model = NULL,
                      nfold = plot_fold, nmax = idw_neighbours,
                      set = list(idp = idw_power), debug.level = 0)
idw_m <- idw_prediction_30$var1.pred
idw_v <- var(cv_idw_30$residual)
# IDW has no model variance: use pooled held-out log-error variance as an approximation.
# This describes error spread and excludes bias; the rankings show bias separately.
# It is constant on the log scale, not generally on the carbon scale.
idw_30 <- export_maps(idw_m, idw_v, "idw", "0_30", testing_folder)

# 2C. Matern kriging: covariance weights ------------------------------

sill_log <- var(test_cores$log_stock)
matern_model_30 <- vgm(psill = sill_log * (1 - matern_nugget_fraction),
                       model = "Mat", range = matern_scale_m,
                       nugget = sill_log * matern_nugget_fraction, kappa = matern_kappa)
# The scale is not a distance cutoff. At kappa=0.5: correlation=exp(-distance/scale).
# Semivariance gamma(h)=0.5*E[(Y(s)-Y(s+h))^2]; for a stationary model gamma(h)=C(0)-C(h).
maximum_distance <- max(as.numeric(st_distance(test_cores)))
cloud <- variogram(log_stock ~ 1, test_cores, cloud = TRUE,
                   cutoff = maximum_distance * 1.01)
semivariance_curve <- variogramLine(matern_model_30, maxdist = maximum_distance)
covariance_curve <- variogramLine(matern_model_30, maxdist = maximum_distance, covariance = TRUE)
curves <- rbind(data.frame(distance_km = covariance_curve$dist / 1000,
                           value = covariance_curve$gamma, quantity = "Covariance"),
                data.frame(distance_km = semivariance_curve$dist / 1000,
                           value = semivariance_curve$gamma, quantity = "Semivariance"))
pairs <- data.frame(distance_km = cloud$dist / 1000, value = cloud$gamma, quantity = "Semivariance")
covariance_plot <- ggplot(curves, aes(distance_km, value)) +
  geom_line(colour = "royalblue", linewidth = 1) +
  geom_point(data = pairs, colour = "grey50") + facet_wrap(~quantity, scales = "free_y") +
  labs(title = "Chosen Matern model and observed core pairs", x = "Distance (km)", y = "Log-stock units squared") +
  theme_minimal()
print(covariance_plot)
ggsave(file.path(testing_folder, "matern_covariance_0_30cm.png"), covariance_plot, width = 10, height = 4)

ok_prediction_30 <- krige(log_stock ~ 1, test_cores, pred_points,
                          model = matern_model_30, debug.level = 0)
cv_ok_30 <- krige.cv(log_stock ~ 1, test_cores, model = matern_model_30,
                     nfold = plot_fold, debug.level = 0)
ok_30 <- export_maps(ok_prediction_30$var1.pred, ok_prediction_30$var1.var,
                     "ok", "0_30", testing_folder)
# krige() supplies map variance; krige.cv() supplies withheld predictions and variances.
# These variances are conditional on the supplied covariance parameters.
# They do not include parameter-estimation uncertainty or uncertainty in regional totals.

# See predictions made while each plot was absent.
evaluation_30 <- rbind(data.frame(method = "Area-based", observed = cv_area_30$observed, predicted = cv_area_30$var1.pred),
                       data.frame(method = "IDW", observed = cv_idw_30$observed, predicted = cv_idw_30$var1.pred),
                       data.frame(method = "Matern kriging", observed = cv_ok_30$observed, predicted = cv_ok_30$var1.pred))
evaluation_plot <- ggplot(evaluation_30, aes(observed, predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(size = 3) +
  facet_wrap(~method) + coord_equal() + theme_minimal() +
  labs(title = "Initial 0-30 cm settings: withheld-plot predictions", x = "Observed log stock", y = "Predicted log stock")
print(evaluation_plot)
ggsave(file.path(testing_folder, "heldout_comparison_0_30cm.png"), evaluation_plot, width = 10, height = 4)

# See how all three methods weight the same cores at one unmeasured grid point.
xy <- st_coordinates(pred_points)
example <- which.min((xy[, 1] - mean(xy[, 1]))^2 + (xy[, 2] - mean(xy[, 2]))^2)
d <- as.numeric(st_distance(test_cores, pred_points[example, ]))
area_weights <- rep(1 / n, n)
idw_weights <- 1 / pmax(d, 1e-10)^idw_power
keep <- head(order(d), min(idw_neighbours, n))
idw_weights[-keep] <- 0
idw_weights <- idw_weights / sum(idw_weights)
D <- as.matrix(st_distance(test_cores))
K <- matrix(variogramLine(matern_model_30, dist_vector = as.vector(D), covariance = TRUE)$gamma, n, n)
k0 <- variogramLine(matern_model_30, dist_vector = d, covariance = TRUE)$gamma
kriging_weights <- solve(rbind(cbind(K, 1), c(rep(1, n), 0)), c(k0, 1))[1:n]
weights <- data.frame(core = rep(test_cores[["Core ID"]], 3),
                       weight = c(area_weights, idw_weights, kriging_weights),
                       method = rep(c("Area-based", "IDW", "Matern kriging"), each = n))
weight_plot <- ggplot(weights, aes(core, weight)) + geom_col(fill = "steelblue") +
  facet_wrap(~method) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Weights at the same prediction point", x = "Core", y = "Weight")
print(weight_plot)
ggsave(file.path(testing_folder, "weights_0_30cm.png"), weight_plot, width = 10, height = 5)

# 3. Compare parameters for 0-30 cm -----------------------------------

# This helper repeats the same held-out prediction calculation across each iteration
# Lower RMSE ranks first. Bias = predicted minus observed, on the log scale.
# Scores compare candidate settings; they are not independent final-model accuracy.
compare_parameters <- function(data, folds) {
  idw_table <- expand.grid(power = idw_powers, neighbours = idw_neighbour_counts)
  for (i in seq_len(nrow(idw_table))) {
    cv <- krige.cv(log_stock ~ 1, data, model = NULL, nfold = folds,
                    nmax = idw_table$neighbours[i], set = list(idp = idw_table$power[i]), debug.level = 0)
    error <- cv$var1.pred - cv$observed
    idw_table$RMSE_log[i] <- sqrt(mean(error^2))
    idw_table$MAE_log[i] <- mean(abs(error))
    idw_table$bias_log[i] <- mean(error)
  }
  idw_table <- idw_table[order(idw_table$RMSE_log), ]
  idw_table <- cbind(rank = seq_len(nrow(idw_table)), idw_table)
  matern_table <- expand.grid(scale_m = matern_scales_m, kappa = matern_kappas)
  for (i in seq_len(nrow(matern_table))) {
    # Fix total sill=1 for ranking: its overall scale cancels from ordinary-kriging weights.
    model <- vgm(1 - matern_nugget_fraction, "Mat", matern_table$scale_m[i],
                 nugget = matern_nugget_fraction, kappa = matern_table$kappa[i])
    cv <- krige.cv(log_stock ~ 1, data, model = model, nfold = folds, debug.level = 0)
    error <- cv$var1.pred - cv$observed
    matern_table$RMSE_log[i] <- sqrt(mean(error^2))
    matern_table$MAE_log[i] <- mean(abs(error))
    matern_table$bias_log[i] <- mean(error)
  }
  matern_table <- matern_table[order(matern_table$RMSE_log), ]
  matern_table <- cbind(rank = seq_len(nrow(matern_table)), matern_table)
  list(idw = idw_table, matern = matern_table)
}

# One plotting helper reuses the same layout for every depth.
plot_rankings <- function(rankings, depth) {
  idw_plot <- ggplot(rankings$idw, aes(factor(power), factor(neighbours), fill = RMSE_log)) +
    geom_tile() + geom_text(aes(label = sprintf("%.3f", RMSE_log))) +
    scale_fill_gradient(low = "#FFF7BC", high = "#D95F0E") +
    scale_y_discrete(labels = function(x) ifelse(x == "Inf", "All", x)) + theme_minimal() +
    labs(title = paste("IDW parameter comparison:", depth_labels[[depth]]), x = "Power", y = "Maximum neighbours", fill = "Log RMSE")
  matern_plot <- ggplot(rankings$matern, aes(factor(scale_m / 1000), factor(kappa), fill = RMSE_log)) +
    geom_tile() + geom_text(aes(label = sprintf("%.3f", RMSE_log))) +
    scale_fill_gradient(low = "#FFF7BC", high = "#D95F0E") + theme_minimal() +
    labs(title = paste("Matern parameter comparison:", depth_labels[[depth]]), x = "Distance scale (km)", y = "Kappa / nu", fill = "Log RMSE")
  print(idw_plot)
  print(matern_plot)
  suffix <- if (depth == "Total") "total" else paste0(depth, "cm")
  ggsave(file.path(ranking_folder, paste0("idw_ranking_", suffix, ".png")), idw_plot, width = 8, height = 5)
  ggsave(file.path(ranking_folder, paste0("matern_ranking_", suffix, ".png")), matern_plot, width = 8, height = 5)
}

rankings_30 <- compare_parameters(test_cores, plot_fold)
print(rankings_30$idw)
print(rankings_30$matern)
write_csv(rankings_30$idw, file.path(ranking_folder, "idw_ranking_0_30cm.csv"))
write_csv(rankings_30$matern, file.path(ranking_folder, "matern_ranking_0_30cm.csv"))
plot_rankings(rankings_30, "0_30")
# Power=0 gives equal weights among selected neighbours; only All gives the area mean.
# Compare Matern scale and kappa together: the same scale has different correlation extents at different kappa.

# 4. Repeat the parameter comparison for the other depths -------------

for (depth in setdiff(names(depth_variables), "0_30")) {
  depth_cores <- core_proj
  depth_cores$log_stock <- log(depth_cores[[depth_variables[[depth]]]])
  rankings <- compare_parameters(depth_cores, plot_fold)
  cat("\nDepth:", depth, "\n")
  print(rankings$idw)
  print(rankings$matern)
  suffix <- if (depth == "Total") "total" else paste0(depth, "cm")
  write_csv(rankings$idw, file.path(ranking_folder, paste0("idw_ranking_", suffix, ".csv")))
  write_csv(rankings$matern, file.path(ranking_folder, paste0("matern_ranking_", suffix, ".csv")))
  plot_rankings(rankings, depth)
}

# 5. REVIEW THE RANKINGS, THEN ENTER YOUR OWN SETTINGS -----------------

# These are editable starting values, NOT automatically selected winners.
# Inf means all available cores. Distances here are metres.
chosen_settings <- data.frame(
  depth = c("0_30", "0_50", "0_100", "Total"),
  idw_power = c(4, 3, 4, 4),
  idw_neighbours = c(Inf, Inf, Inf, Inf),
  matern_scale_m = c(40000, 80000, 40000, 40000),
  matern_kappa = c(2.5, 0.5, 2.5, 2.5)
)
print(chosen_settings)

# 6. Run YOUR settings and export the final maps ----------------------

# Run this section after editing chosen_settings. It replaces the final Script 2 products.
output_folder <- "outputs/script2"
write_csv(chosen_settings, file.path(output_folder, "chosen_settings.csv"))
for (i in seq_len(nrow(chosen_settings))) {
  settings <- chosen_settings[i, ]
  depth <- settings$depth
  depth_cores <- core_proj
  depth_cores$log_stock <- log(depth_cores[[depth_variables[[depth]]]])
  y <- depth_cores$log_stock
  suffix <- if (depth == "Total") "total" else paste0(depth, "cm")
  area_heldout <- numeric(length(y))
  area_heldout_variance <- numeric(length(y))
  for (fold in unique(plot_fold)) {
    training_y <- y[plot_fold != fold]
    area_heldout[plot_fold == fold] <- mean(training_y)
    area_heldout_variance[plot_fold == fold] <- var(training_y) * (1 + 1 / length(training_y))
  }
  cv_area <- data.frame(observed = y, var1.pred = area_heldout, var1.var = area_heldout_variance)
  cv_idw <- krige.cv(log_stock ~ 1, depth_cores, model = NULL, nfold = plot_fold,
                     nmax = settings$idw_neighbours, set = list(idp = settings$idw_power), debug.level = 0)
  idw_prediction <- idw(log_stock ~ 1, depth_cores, pred_points,
                        idp = settings$idw_power, nmax = settings$idw_neighbours, debug.level = 0)
  sill_log <- var(y)
  model <- vgm(sill_log * (1 - matern_nugget_fraction), "Mat", settings$matern_scale_m,
               nugget = sill_log * matern_nugget_fraction, kappa = settings$matern_kappa)
  cv_ok <- krige.cv(log_stock ~ 1, depth_cores, model = model, nfold = plot_fold, debug.level = 0)
  ok_prediction <- krige(log_stock ~ 1, depth_cores, pred_points, model = model, debug.level = 0)
  log_means <- list(area_based = rep(mean(y), nrow(pred_points)),
                    idw = idw_prediction$var1.pred, ok = ok_prediction$var1.pred)
  log_variances <- list(area_based = var(y) * (1 + 1 / length(y)),
                        idw = var(cv_idw$residual), ok = ok_prediction$var1.var)
  heldout <- list(area_based = cv_area, idw = cv_idw, ok = cv_ok)
  variance_sources <- c(area_based = "IID constant-mean new-location prediction variance: s_log^2*(1+1/n)",
                        idw = "Pooled withheld-plot log-error variance; working approximation; no bias correction",
                        ok = "Conditional Matern kriging log variance; mean across mapped cells recorded below")
  for (method in c("area_based", "idw", "ok")) {
    maps <- export_maps(log_means[[method]], log_variances[[method]], method, depth, output_folder)
    cv <- heldout[[method]]
    error <- cv$var1.pred - cv$observed
    # A compact results row keeps Script 4's existing reader up to date.
    calibration <- data.frame(method = method, depth = depth, n_cores = length(y),
      n_plots = length(unique(plot_fold)), power = if (method == "area_based") 0 else if (method == "idw") settings$idw_power else NA_real_,
      neighbours = if (method == "idw") settings$idw_neighbours else Inf,
      log_error_bias = mean(error), log_error_variance = var(error),
      prediction_log_variance = mean(log_variances[[method]]), variance_source = variance_sources[[method]],
      RMSE_log = sqrt(mean(error^2)), MAE_log = mean(abs(error)),
      matern_scale_m = if (method == "ok") settings$matern_scale_m else NA_real_,
      matern_kappa = if (method == "ok") settings$matern_kappa else NA_real_,
      nugget_fraction = if (method == "ok") matern_nugget_fraction else NA_real_)
    write_csv(calibration, file.path(output_folder, paste0(method, "_calibration_", suffix, ".csv")))
    write_csv(data.frame(core_id = depth_cores[["Core ID"]], plot_id = depth_cores[["Plot ID"]],
                         observed_log = cv$observed, predicted_log = cv$var1.pred, log_error = error,
                         prediction_log_variance = cv$var1.var),
              file.path(output_folder, paste0(method, "_heldout_", suffix, ".csv")))
    print(calibration[, c("method", "depth", "RMSE_log", "MAE_log", "log_error_bias")])
  }
}
# Final files use area_based / idw / ok labels expected by Scripts 4 and 5.
# Carbon maps: mean (kg C/m2), SD (kg C/m2), variance ((kg C/m2)^2).
# Log-moment bands stay named log_mean and log_variance for Script 5.
# Neither pixel SD nor pixel variance alone gives uncertainty in a regional total.


# 7. Summarize uncertainty and the distribution of map values ---------

# Re-read the final exports so these tables describe the saved maps.
# Stock variance uses (kg C/m2)^2. Stock SD uses kg C/m2.
# Log-error variance describes the spread of held-out errors; MSE also includes bias.
# R var() divides by n-1, while MSE divides by n: MSE=(n-1)/n*var(error)+bias^2.
# Spatial variation in map means is different from uncertainty at a single cell.
uncertainty_rows <- list()
distribution_rows <- list()
sd_distribution_30 <- list()
for (depth in chosen_settings$depth) {
  stock <- core_proj[[depth_variables[[depth]]]]
  suffix <- if (depth == "Total") "total" else paste0(depth, "cm")
  for (method in c("area_based", "idw", "ok")) {
    key <- paste(depth, method, sep = "_")
    map_mean <- as.numeric(values(rast(file.path(output_folder, paste0(method, "_", suffix, ".tif")))))
    map_sd <- as.numeric(values(rast(file.path(output_folder, paste0(method, "_sd_", suffix, ".tif")))))
    map_variance <- as.numeric(values(rast(file.path(output_folder, paste0(method, "_variance_", suffix, ".tif")))))
    cv_stats <- read_csv(file.path(output_folder, paste0(method, "_calibration_", suffix, ".csv")),
                         show_col_types = FALSE)
    uncertainty_rows[[key]] <- data.frame(
      depth = depth, method = method, n_cores = cv_stats$n_cores,
      observed_stock_variance = var(stock), observed_log_variance = var(log(stock)),
      spatial_variance_of_map_means = var(map_mean, na.rm = TRUE),
      mean_prediction_variance = mean(map_variance, na.rm = TRUE),
      mean_prediction_sd = mean(map_sd, na.rm = TRUE),
      log_error_variance = cv_stats$log_error_variance,
      MSE_log = cv_stats$RMSE_log^2, RMSE_log = cv_stats$RMSE_log,
      MAE_log = cv_stats$MAE_log, bias_log = cv_stats$log_error_bias,
      working_distribution = "Lognormal", variance_source = cv_stats$variance_source)

    # Summaries across cells: these quantiles are NOT a prediction interval at one cell.
    quantities <- list(map_mean = map_mean, prediction_sd = map_sd, prediction_variance = map_variance)
    units <- c(map_mean = "kg C/m2", prediction_sd = "kg C/m2", prediction_variance = "(kg C/m2)^2")
    for (quantity in names(quantities)) {
      z <- quantities[[quantity]]
      percentiles <- quantile(z, c(0.05, 0.5, 0.95), na.rm = TRUE)
      distribution_rows[[paste(key, quantity, sep = "_")]] <- data.frame(
        depth = depth, method = method, quantity = quantity, units = units[[quantity]],
        mean = mean(z, na.rm = TRUE), spatial_sd = sd(z, na.rm = TRUE),
        min = min(z, na.rm = TRUE), p05 = unname(percentiles[1]), median = unname(percentiles[2]),
        p95 = unname(percentiles[3]), max = max(z, na.rm = TRUE))
    }
    if (depth == "0_30") sd_distribution_30[[method]] <- data.frame(
      method = method_labels[[method]], prediction_sd = map_sd[is.finite(map_sd)])
  }
}
uncertainty_summary <- do.call(rbind, uncertainty_rows)
map_distribution_summary <- do.call(rbind, distribution_rows)
print(uncertainty_summary)
print(map_distribution_summary)
write_csv(uncertainty_summary, file.path(output_folder, "uncertainty_summary.csv"))
write_csv(map_distribution_summary, file.path(output_folder, "map_distribution_summary.csv"))

# Distribution of cell-level SDs across the 0-30 cm map.
# This histogram is not the carbon-stock distribution for any individual cell.
sd_distribution_plot <- ggplot(do.call(rbind, sd_distribution_30), aes(prediction_sd)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
  facet_wrap(~method, scales = "free_y") + theme_minimal() +
  labs(title = "Distribution of prediction SD across cells: 0-30 cm",
       x = "Pointwise prediction SD (kg C/m2)", y = "Number of cells")
print(sd_distribution_plot)
ggsave(file.path(output_folder, "prediction_sd_distribution_0_30cm.png"),
       sd_distribution_plot, width = 10, height = 4)
