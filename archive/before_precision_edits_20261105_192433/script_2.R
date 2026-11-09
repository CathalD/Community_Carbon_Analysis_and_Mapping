
# ==============================================

# Part 2 - Spatial Interpolation

# ==============================================

library(sf)
library(terra)
library(gstat)
library(dplyr)
library(readr)
library(ggplot2)



# ===============
# Step 1 - Set-up the mapping settings
#=================
# Set-up common settings
analysis_crs <- "EPSG:3978"
coordinate_crs <- "EPSG:4326"
grid_folder <- "data/grid"

# Input data
mapping_file <- "outputs/script1/mapping_data.csv" # This comes from script 1
aoi_file <- "data/boundary/aoi.geojson" # Boundary file added by user


# Set-up output folder
output_folder <- "outputs/script2"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

#Other settings:
# Adjust these limits together for all maps if results exceed them.
mapping_depths_cm <- c(30, 50, 100)
cell_size_m <- 1000
source("config.R")
# This function creates a blank grid over the boundary AOI using the cellsize defined above
source("helpers/grid_helper.R")

# Map colours and display limits
map_colours <- hcl.colors(100, palette = "Viridis")
mean_limits <- c(0, 100)       # kg C/m²
variance_limits <- c(0, 2500)  # (kg C/m²)²




# ===================
# Step 2 - Import organize data
#==============

cores <- read_csv(mapping_file, show_col_types = FALSE) %>%
  rename(core_id = `Core ID`, plot_id = `Plot ID`,
         latitude = Latitude, longitude = Longitude)
terra::projNetwork(FALSE)
core_proj <- st_as_sf(project(
  vect(cores, geom = c("longitude", "latitude"), crs = coordinate_crs), analysis_crs
))

# 2. One saved grid, mask and set of cell centres
shared <- create_analysis_grid()
r <- shared$template
analysis_mask <- shared$mask
aoi_proj <- st_as_sf(shared$aoi)
pred_points <- st_as_sf(as.points(shared$cell_ids, values = TRUE, na.rm = TRUE))
grid_sf <- st_as_sf(as.polygons(
  shared$cell_ids, aggregate = FALSE, values = TRUE, na.rm = TRUE
))
st_write(aoi_proj, file.path(output_folder, "aoi_polygon.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
depth_variables <- setNames(
  paste0("C Stock 0-", mapping_depths_cm, "cm"), paste0("0_", mapping_depths_cm)
)
depth_variables <- c(depth_variables, Total = "Total Carbon Stock")



# ================================
# Step 3 - Run Spatial interpolation

# Note: Carbon stocks are always log transformed when modelling - This accounts for 2 things
# 1. The data cannot be negative (you can't measure a negative amount of carbon)
# 2. SOil carbon data tends to be skewed, this often corrects it


# First lets run 1 and then we can loop it

depth <- "0_100"
data <- sf::st_drop_geometry(core_proj)
core_proj$stock <- data[[depth_variables[[depth]]]]
core_proj$log_stock <- log(data$stock)
print(core_proj)


idw_1 <- gstat::idw(formula = log_stock ~ 1, locations = core_proj, newdata = pred_points, nmax = Inf, idp = 0)

plot(idw_1)


# So, we first we create a single function that will handle our log-transformed data
# Convert log-scale to carbon prediction and variance.
export_maps <- function(log_mean, log_variance, method, depth) { # AI helped write this function
  
  log_variance <- rep_len(log_variance, nrow(pred_points))
                                                                  #First convert data to natural log
  moments <- data.frame(
    prediction = exp(log_mean + log_variance / 2),
    variance = expm1(log_variance) *
      exp(2 * log_mean + log_variance),
    log_mean = log_mean,
    log_variance = log_variance
  )
  
  maps <- rast(r, nlyrs = ncol(moments))
  names(maps) <- names(moments)
  
  cells <- matrix(NA_real_, ncell(r), ncol(moments))
  cells[pred_points$grid_id, ] <- as.matrix(moments)
  values(maps) <- cells
  
  maps <- mask(maps, analysis_mask)
  
  suffix <- if (depth == "Total") {
    "total"
  } else {
    paste0(depth, "cm")
  }
  
  writeRaster(
    maps[["prediction"]],
    file.path(output_folder, paste0(method, "_", suffix, ".tif")),
    overwrite = TRUE
  )
  
  writeRaster(
    maps[["variance"]],
    file.path(
      output_folder,
      paste0(method, "_variance_", suffix, ".tif")
    ),
    overwrite = TRUE
  )
  
  writeRaster(
    maps[[c("log_mean", "log_variance")]],
    file.path(
      output_folder,
      paste0(method, "_log_moments_", suffix, ".tif")
    ),
    datatype = "FLT8S",
    overwrite = TRUE
  )
  
  # Plot predicted carbon and variance
  draw_maps <- function() {
    
    previous_settings <- par(mfrow = c(1, 2))
    on.exit(par(previous_settings))
    
    terra::plot(
      maps[["prediction"]],
      col = map_colours,
      range = mean_limits,
      main = paste(method, depth, "mean (kg C/m²)")
    )
    
    terra::plot(
      maps[["variance"]],
      col = map_colours,
      range = variance_limits,
      main = "Predictive variance ((kg C/m²)²)"
    )
  }
  
  print(global(maps[[1:2]], c("min", "max"), na.rm = TRUE))
  
  # Save the figure.
  png(
    file.path(output_folder, paste0(method, "_", suffix, ".png")),
    width = 1600,
    height = 800,
    res = 150
  )
  
  draw_maps()
  dev.off()
  
  # Display the same figure in RStudio.
  draw_maps()
  
  invisible(maps)
}


# Export the variance as well
export_calibration <- function(
    data, predicted, method, depth, power, neighbours) {             #AI helped write this function
  
  error <- predicted - data$log_stock
  
  heldout <- data.frame(
    core_id = data$core_id,
    plot_id = data$plot_id,
    observed = data$stock,
    observed_log = data$log_stock,
    predicted_log = predicted,
    log_error = error
  )
  
  calibration <- data.frame(
    method = method,
    depth = depth,
    n_cores = nrow(data),
    n_plots = length(unique(data$plot_id)),
    power = power,
    neighbours = neighbours,
    log_error_bias = mean(error),
    log_error_variance = var(error),
    RMSE_log = sqrt(mean(error^2)),
    scale = "log C; calibration parameters treated as fixed",
    depth_convention =
      "Script 1 stocks; completed profile at sampled maximum assumed"
  )
  
  print(calibration)
  
  suffix <- if (depth == "Total") {
    "total"
  } else {
    paste0(depth, "cm")
  }
  
  write_csv(
    calibration,
    file.path(
      output_folder,
      paste0(method, "_calibration_", suffix, ".csv")
    )
  )
  
  write_csv(
    heldout,
    file.path(
      output_folder,
      paste0(method, "_heldout_", suffix, ".csv")
    )
  )
  
  calibration
}


# =======================================
# Analysis 1 - Area-based mapping
# ================================

# All of the spatial interpolation methods follow the same general formula

# C(Si) = E W * C (So)
# C(Si) = Carbon stock prediction at location i
# W(i) = Weight given to observation C(So)
# C(So) = Observed carbon stock at location o

# The package gstat applies using the function gstat::idw()

# 1. Settings
depth <- "0_100"
area_power <- 0
area_neighbours <- Inf

area_polygons <- file.path(output_folder, "area_based.gpkg")

# ================================
# 2. Prepare observations
# ================================

# Keep the projected coordinates and Script 1 stocks.
area_cores <- core_proj

area_cores$stock <- area_cores[[depth_variables[[depth]]]]

# Natural-log transformation.
area_cores$log_stock <- log(area_cores$stock)

print(
  sf::st_drop_geometry(area_cores)[
    , c("core_id", "plot_id", "stock", "log_stock")
  ]
)


# ================================
# 3. Predict at the shared grid centres
# ================================

# Power zero and all neighbours give equal weights
# at grid locations that do not coincide with observed cores.
idw_1 <- gstat::idw(
  formula = log_stock ~ 1,
  locations = area_cores,
  newdata = pred_points,
  nmax = area_neighbours,
  idp = area_power,
  debug.level = 0
)

# These predictions are on the natural-log scale.
predicted_log <- idw_1$var1.pred

print(summary(predicted_log))


# ================================
# 4. Estimate held-out prediction errors
# ================================

heldout_predictions <- rep(NA_real_, nrow(area_cores))

for (plot_id in unique(area_cores$plot_id)) {
  
  test <- area_cores$plot_id == plot_id
  
  heldout_fit <- gstat::idw(
    formula = log_stock ~ 1,
    locations = area_cores[!test, ],
    newdata = area_cores[test, ],
    nmax = area_neighbours,
    idp = area_power,
    debug.level = 0
  )
  
  heldout_predictions[test] <- heldout_fit$var1.pred
}

# Error = predicted minus observed log stock.
log_error <- heldout_predictions - area_cores$log_stock
log_error_bias <- mean(log_error)
log_error_variance <- var(log_error)
RMSE_log <- sqrt(mean(log_error^2))

print(data.frame(
  core_id = area_cores$core_id,
  observed_log = area_cores$log_stock,
  predicted_log = heldout_predictions,
  log_error = log_error
))

print(data.frame(
  depth = depth,
  cores = nrow(area_cores),
  log_error_bias = log_error_bias,
  log_error_variance = log_error_variance,
  RMSE_log = RMSE_log
))


# ================================
# 5. Define the working log predictions
# ================================

# Correct the estimated average log error.
log_mean <- predicted_log - log_error_bias
log_variance <- log_error_variance


# ================================
# 6. Convert to carbon units
# ================================

median_stock <- exp(log_mean)

mean_stock <- exp(log_mean + log_variance / 2)

variance_stock <- mean_stock^2 * expm1(log_variance)

sd_stock <- sqrt(variance_stock)

# Area-based results are constant across this grid.
print(data.frame(
  depth = depth,
  cores = nrow(area_cores),
  median_kg_C_m2 = median_stock[1],
  mean_kg_C_m2 = mean_stock[1],
  SD_kg_C_m2 = sd_stock[1],
  variance = variance_stock[1]
))


# ================================
# 7. Plot and export the rasters
# ================================

# Uses the existing shared-grid export function,
# colour palette and colour limits.
maps <- export_maps(
  log_mean = log_mean,
  log_variance = log_variance,
  method = "area_based",
  depth = depth
)

plot(maps)

# ================================
# 8. Export the error calculations
# ================================

calibration <- export_calibration(
  data = area_cores,
  predicted = heldout_predictions,
  method = "area_based",
  depth = depth,
  power = area_power,
  neighbours = area_neighbours
)


# ================================
# 9. Export the polygon result
# ================================

area <- aoi_proj
area$depth <- depth
area$method <- "area_based"
area$mean_stock_kg_m2 <- mean_stock[1]

sf::st_write(
  area,
  area_polygons,
  layer = paste0("area_based_", depth),
  delete_layer = TRUE,
  quiet = TRUE
)



# ==================================
# Now we just loop this for the other depths

# ================================
# Area-based mapping: all depths
# ================================

# Settings
area_power <- 0
area_neighbours <- Inf
area_polygons <- file.path(output_folder, "area_based.gpkg")

# Store each depth's rasters for use afterwards.
area_maps <- list()


for (depth in names(depth_variables)) {
  
  # ================================
  # 1. Prepare observations
  # ================================
  
  area_cores <- core_proj
  
  area_cores$stock <- area_cores[[depth_variables[[depth]]]]
  
  area_cores$log_stock <- log(area_cores$stock)
  
  cat(
    "\nArea-based:", depth,
    "| Cores:", nrow(area_cores),
    "| Distance power:", area_power, "\n"
  )
  
  
  # ================================
  # 2. Predict at the shared grid centres
  # ================================
  
  idw_1 <- gstat::idw(
    formula = log_stock ~ 1,
    locations = area_cores,
    newdata = pred_points,
    nmax = area_neighbours,
    idp = area_power,
    debug.level = 0
  )
  
  predicted_log <- idw_1$var1.pred
  
  
  # ================================
  # 3. Estimate held-out prediction errors
  # ================================
  
  heldout_predictions <- rep(NA_real_, nrow(area_cores))
  
  for (plot_id in unique(area_cores$plot_id)) {
    
    test <- area_cores$plot_id == plot_id
    
    heldout_fit <- gstat::idw(
      formula = log_stock ~ 1,
      locations = area_cores[!test, ],
      newdata = area_cores[test, ],
      nmax = area_neighbours,
      idp = area_power,
      debug.level = 0
    )
    
    heldout_predictions[test] <- heldout_fit$var1.pred
  }
  
  log_error <- heldout_predictions - area_cores$log_stock
  
  log_error_bias <- mean(log_error)
  log_error_variance <- var(log_error)
  RMSE_log <- sqrt(mean(log_error^2))
  
  
  # ================================
  # 4. Define the working log predictions
  # ================================
  
  log_mean <- predicted_log - log_error_bias
  
  # Retain the current held-out-error approximation.
  log_variance <- log_error_variance
  
  
  # ================================
  # 5. Convert to carbon units
  # ================================
  
  median_stock <- exp(log_mean)
  
  mean_stock <- exp(log_mean + log_variance / 2)
  
  variance_stock <- mean_stock^2 * expm1(log_variance)
  
  sd_stock <- sqrt(variance_stock)
  
  print(data.frame(
    depth = depth,
    cores = nrow(area_cores),
    median_kg_C_m2 = median_stock[1],
    mean_kg_C_m2 = mean_stock[1],
    SD_kg_C_m2 = sd_stock[1],
    variance = variance_stock[1],
    log_error_variance = log_error_variance,
    RMSE_log = RMSE_log
  ))
  
  
  # ================================
  # 6. Plot and export the rasters
  # ================================
  
  maps <- export_maps(
    log_mean = log_mean,
    log_variance = log_variance,
    method = "area_based",
    depth = depth
  )
  
  area_maps[[depth]] <- maps
  
  
  # ================================
  # 7. Export the error calculations
  # ================================
  
  calibration <- export_calibration(
    data = area_cores,
    predicted = heldout_predictions,
    method = "area_based",
    depth = depth,
    power = area_power,
    neighbours = area_neighbours
  )
  
  
  # ================================
  # 8. Export the polygon result
  # ================================
  
  area <- aoi_proj
  area$depth <- depth
  area$method <- "area_based"
  area$mean_stock_kg_m2 <- mean_stock[1]
  
  sf::st_write(
    area,
    area_polygons,
    layer = paste0("area_based_", depth),
    delete_layer = TRUE,
    quiet = TRUE
  )
  
} # End depth loop



# ================================
# Analysis 2 - IDW mapping
# ================================
# Here build of this, but we are adjusting the parameters
# Now the distance from the prediction to the field sampled is directly
# weighted to this distance, for each prediction point
# We test this for each depth and solve for the least amount of RMSE

# Settings
idw_powers <- c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5)
idw_neighbours <- c(2, 3, 4, 5)

tune_idw <- TRUE


# IDW prediction function
predict_idw_log_stock <- function(data, newdata, power, neighbours) {
  
  gstat::idw(
    log_stock ~ 1,
    locations = data,
    newdata = newdata,
    idp = power,
    nmax = min(neighbours, nrow(data)),
    debug.level = 0
  )$var1.pred
}


# Predictions with each whole plot withheld.
heldout_idw_predictions <- function(data, power, neighbours) {
  
  predicted <- rep(NA_real_, nrow(data))
  
  for (plot in unique(data$plot_id)) {
    
    test <- data$plot_id == plot
    
    predicted[test] <- predict_idw_log_stock(
      data = data[!test, ],
      newdata = data[test, ],
      power = power,
      neighbours = neighbours
    )
  }
  
  predicted
}


# Select IDW settings using held-out log-scale RMSE.
choose_idw <- function(data) {
  
  if (!tune_idw || length(unique(data$plot_id)) < 2) {
    
    return(data.frame(
      idp = fixed_idw_power,
      nmax = min(fixed_idw_neighbours, nrow(data)),
      RMSE_log = NA_real_
    ))
  }
  
  parameters <- expand.grid(
    idp = idw_powers,
    nmax = idw_neighbours
  )
  
  parameters <- parameters[
    parameters$nmax <= nrow(data) - 1,
    ,
    drop = FALSE
  ]
  
  if (!nrow(parameters)) {
    
    return(data.frame(
      idp = fixed_idw_power,
      nmax = 1,
      RMSE_log = NA_real_
    ))
  }
  
  parameters$RMSE_log <- vapply(
    seq_len(nrow(parameters)),
    function(i) {
      
      predicted <- heldout_idw_predictions(
        data,
        power = parameters$idp[i],
        neighbours = parameters$nmax[i]
      )
      
      sqrt(mean((predicted - data$log_stock)^2))
    },
    numeric(1)
  )
  
  parameters[which.min(parameters$RMSE_log), , drop = FALSE]
}


settings_by_depth <- list()

for (depth in names(depth_variables)) {
  
  # Prepare the same Script 1 stocks on the same projected coordinates.
  data <- core_proj
  data$stock <- data[[depth_variables[[depth]]]]
  data <- data[is.finite(data$stock) & data$stock > 0, ]
  data$log_stock <- log(data$stock)
  
  cat("\nIDW:", depth, "| Cores:", nrow(data), "\n")
  
  if (length(unique(data$plot_id)) < 2) {
    warning("At least two plots are needed for error calibration.")
    next
  }
  
  # Select settings for the final map using all retained cores.
  selected <- choose_idw(data)
  print(selected)
  
  settings_by_depth[[depth]] <- data.frame(
    Depth = depth,
    Variable = if (depth == "Total") {
      "stock_Total"
    } else {
      paste0("stock_", depth)
    },
    N_Cores = nrow(data),
    idp = selected$idp,
    nmax = selected$nmax,
    RMSE_log = selected$RMSE_log
  )
  
  # Estimate errors on outer held-out plots.
  # Retune within each training fold to keep test plots out of tuning.
  heldout_predictions <- rep(NA_real_, nrow(data))
  
  for (plot in unique(data$plot_id)) {
    
    test <- data$plot_id == plot
    training <- data[!test, ]
    
    fold_settings <- choose_idw(training)
    
    heldout_predictions[test] <- predict_idw_log_stock(
      data = training,
      newdata = data[test, ],
      power = fold_settings$idp,
      neighbours = fold_settings$nmax
    )
  }
  
  calibration <- export_calibration(
    data = data,
    predicted = heldout_predictions,
    method = "idw",
    depth = depth,
    power = selected$idp,
    neighbours = selected$nmax
  )
  
  # Predict at the same shared grid centres.
  predicted_log <- predict_idw_log_stock(
    data = data,
    newdata = pred_points,
    power = selected$idp,
    neighbours = selected$nmax
  )
  
  maps <- export_maps(
    log_mean = predicted_log - calibration$log_error_bias,
    log_variance = calibration$log_error_variance,
    method = "idw",
    depth = depth
  )
}

best_by_depth <- bind_rows(settings_by_depth)

write_csv(
  best_by_depth,
  file.path(output_folder, "idw_best_parameters_by_depth.csv")
)


# ================================
# Step 3.3 Kriging
# ================================

# Building off this - Now we replace the inverse distance weighting with
# a measure of semi-variance between each prediction point and the observed data
# We then model this
#Start with the 30 cm depth
# 1. Settings
depth <- "0_30"

variogram_models <- c("Sph", "Exp", "Gau")
variogram_distance_classes <- 5
variogram_cutoff_fraction <- 1

# Store each depth's rasters.
ok_maps <- list()


# ================================
# 2. Prepare observations
# ================================

for (depth in setdiff(names(depth_variables), "0_30")) {
  
  ok_cores <- core_proj
  
  ok_cores$stock <- ok_cores[[depth_variables[[depth]]]]
  
  ok_cores$log_stock <- log(ok_cores$stock)

cat(
  "\nOrdinary kriging:", depth,
  "| Cores:", nrow(ok_cores), "\n"
)

# Depth-specific filenames prevent overwriting other depths.
suffix <- if (depth == "Total") {
  "total"
} else {
  paste0(depth, "cm")
}


# ================================
# 3. Calculate the variogram
# ================================

maximum_distance <- max(
  as.numeric(sf::st_distance(ok_cores))
)

cutoff_distance <- maximum_distance * variogram_cutoff_fraction

# Semivariance for individual pairs of observed cores.
ok_cloud <- gstat::variogram(
  log_stock ~ 1,
  locations = ok_cores,
  cloud = TRUE,
  cutoff = cutoff_distance * 1.000001
)

# Average semivariance within distance classes.
ok_variogram <- gstat::variogram(
  log_stock ~ 1,
  locations = ok_cores,
  cutoff = cutoff_distance * 1.000001,
  width = cutoff_distance / variogram_distance_classes
)

# np = number of pairs contributing to each bin.
print(ok_variogram)


# ================================
# 4. Plot the cloud and bin averages
# ================================

cloud_plot <- ggplot() +
  
  geom_point(
    data = as.data.frame(ok_cloud),
    aes(x = dist / 1000, y = gamma),
    colour = "grey55",
    alpha = 0.6,
    size = 2
  ) +
  
  geom_point(
    data = as.data.frame(ok_variogram),
    aes(x = dist / 1000, y = gamma),
    colour = "#D95F02",
    size = 3.5
  ) +
  
  labs(
    title = paste("Log-carbon variogram:", depth),
    subtitle = "Grey: individual pairs. Orange: distance-bin averages.",
    x = "Distance between cores (km)",
    y = "Semivariance (log-stock units²)"
  ) +
  
  theme_minimal()

print(cloud_plot)

ggsave(
  file.path(
    output_folder,
    paste0("kriging_variogram_cloud_", suffix, ".png")
  ),
  cloud_plot,
  width = 9,
  height = 6,
  dpi = 200
)


# ================================
# 5. Fit and display the variogram model
# ================================

ok_model <- gstat::fit.variogram(
  ok_variogram,
  model = gstat::vgm(variogram_models)
)

print(ok_model)

# Read any fitting warnings printed in the console.
curve_data <- gstat::variogramLine(
  ok_model,
  maxdist = cutoff_distance,
  n = 300
)

fitted_plot <- cloud_plot +
  
  geom_line(
    data = curve_data,
    aes(x = dist / 1000, y = gamma),
    colour = "#2166AC",
    linewidth = 1
  )

print(fitted_plot)

ggsave(
  file.path(
    output_folder,
    paste0("kriging_variogram_fit_", suffix, ".png")
  ),
  fitted_plot,
  width = 9,
  height = 6,
  dpi = 200
)


# ================================
# 6. Predict at the shared grid centres
# ================================

ok_prediction <- gstat::krige(
  formula = log_stock ~ 1,
  locations = ok_cores,
  newdata = pred_points,
  model = ok_model
)

# Both quantities are on the natural-log scale.
log_mean <- ok_prediction$var1.pred
log_variance <- ok_prediction$var1.var

# Remove tiny negative values caused by numerical rounding.
log_variance[
  log_variance < 0 & log_variance > -1e-8
] <- 0

# Larger negative variances cannot be interpreted.
log_variance[log_variance < 0] <- NA_real_

cat(
  "Cells without valid log variance:",
  sum(!is.finite(log_variance)), "\n"
)

print(summary(data.frame(
  log_prediction = log_mean,
  log_prediction_variance = log_variance
)))


# ================================
# 7. Convert, plot and export
# ================================

# The existing function converts log moments to
# carbon-scale means and variances, then plots and exports.
ok_rasters <- export_maps(
  log_mean = log_mean,
  log_variance = log_variance,
  method = "ok",
  depth = depth
)

# Save the empirical variogram and fitted parameters.
readr::write_csv(
  as.data.frame(ok_variogram),
  file.path(
    output_folder,
    paste0("kriging_empirical_variogram_", suffix, ".csv")
  )
)

readr::write_csv(
  as.data.frame(ok_model),
  file.path(
    output_folder,
    paste0("kriging_variogram_model_", suffix, ".csv")
  )
)

saveRDS(
  ok_model,
  file.path(
    output_folder,
    paste0("kriging_variogram_model_", suffix, ".rds")
  )
)

ok_maps[[depth]] <- ok_rasters

} # End depth loop
