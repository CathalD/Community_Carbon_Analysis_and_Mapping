#==========================================
#----- Script 5: Bayesian carbon mapping ---
#==========================================

library(terra)
library(readr)
library(dplyr)
library(ggplot2)

#==========================================
#----- 5A: Settings and inputs -------------
#==========================================

target_depth_cm <- 100

stock_column <- paste0(
  "C Stock 0-", target_depth_cm, "cm"
)

output_folder <- file.path(
  "outputs/Testing_Bayes_Models",
  paste0("stock_0_", target_depth_cm)
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

#------------------------------------
input_files <- c(
  field = "outputs/script1/mapping_data.csv",
  boundary = "data/boundary/aoi.geojson",
  prior_mean = "outputs/script4/mcmaster_0_100cm.tif",
  prior_uncertainty = "outputs/script4/mcmaster_uncertainty_0_100cm.tif",
  grid = "outputs/script2/idw_0_100cm.tif",
  area_based_map = "outputs/script2/area_based.gpkg"
)

# Make a new blank grid:
# First build a grid using terra
#----- Create a new 1 km analysis grid -----
terra::projNetwork(FALSE)
aoi <- vect(input_files[["boundary"]])
aoi <- project(aoi, "EPSG:3978")

template <- rast(
  xmin = floor(xmin(aoi) / 1000) * 1000,
  xmax = ceiling(xmax(aoi) / 1000) * 1000,
  ymin = floor(ymin(aoi) / 1000) * 1000,
  ymax = ceiling(ymax(aoi) / 1000) * 1000,
  resolution = 1000,
  crs = crs(aoi)
)

#----- Convert grid cells to polygons -----

grid_sf <- sf::st_as_sf(
  as.polygons(
    template,
    aggregate = FALSE,
    values = FALSE,
    na.rm = FALSE
  )
)

aoi_sf <- sf::st_as_sf(aoi)

#----- Plot the grid and boundary -----

ggplot() +
  geom_sf(
    data = grid_sf,
    fill = "lightblue",
    colour = "grey80",
    linewidth = 0.2
  ) +
  geom_sf(
    data = aoi_sf,
    fill = NA,
    colour = "black",
    linewidth = 0.5
  ) +
  labs(title = "New 1 km analysis grid") +
  theme_minimal()

#--------------------------------------------


raw_cores <- read_csv(
  input_files[["field"]],
  show_col_types = FALSE
)

required_columns <- c(
  "Core ID", "Plot ID",
  "Longitude", "Latitude",
  stock_column
)

cores <- data.frame(
  core_id = as.character(raw_cores[["Core ID"]]),
  plot_id = as.character(raw_cores[["Plot ID"]]),
  longitude = raw_cores[["Longitude"]],
  latitude = raw_cores[["Latitude"]],
  stock = raw_cores[[stock_column]]
)

# -----------------------------
prior_mean <- rast(
  input_files[["prior_mean"]]
)

prior_uncertainty <- rast(
  input_files[["prior_uncertainty"]]
)
names(prior_mean) <- "prior_mean"
names(prior_uncertainty) <- "prior_interval_width_90"
# Prepare the conversion in 5B.
# These are interval widths, in kg C/m².
prior_uncertainty_type <- "pi90_width"
prior_mean
prior_uncertainty

#------------------------------

# Use the new grid and boundary created above.

core_points <- vect(
  cores,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

core_points <- project(
  core_points,
  crs(template)
)

print(cores[, c("core_id", "plot_id", "stock")])

cat("\nTarget depth: 0–", target_depth_cm, " cm\n", sep = "")
cat("Field cores:", nrow(cores), "\n")
cat("Sampling plots:", length(unique(cores$plot_id)), "\n")
cat("Grid resolution:", res(template), "metres\n")
cat("Stock units: kg C/m²\n")

# ===================================================================

# Prepare the prior data

# Convert Quantile Uncertainty to SD
prior_sd <- prior_uncertainty / (2 * qnorm(0.95))

# create a new dataframe with update mean and sd
native_prior <- c(prior_mean, prior_sd)

names(native_prior) <- c(
  "prior_mean",
  "prior_sd"
)

# Extract values at core points
at_cores <- terra::extract(
  native_prior,
  project(core_points, crs(native_prior)),
  method = "simple",
  ID = FALSE
)

at_cores


#-=============
# Resample the prior are the analysis grid scale

# Keep each prior mean and SD paired when transferring to the grid.
grid_prior <- project(
  native_prior,
  template,
  method = "near"
)
# MAke sure its in the correct crs
grid_prior <- mask(
  grid_prior,
  project(aoi, crs(grid_prior)),
  touches = FALSE
)
#MAsk any zero value
valid_grid <- (
  is.finite(grid_prior[["prior_mean"]]) &
    grid_prior[["prior_mean"]] > 0 &
    is.finite(grid_prior[["prior_sd"]]) &
    grid_prior[["prior_sd"]] > 0
)

grid_prior <- mask(
  grid_prior,
  valid_grid,
  maskvalues = 0
)

writeRaster(
  grid_prior,
  file.path(output_folder, "prepared_prior.tif"),
  overwrite = TRUE
)
grid_prior_df <- as.data.frame(grid_prior, xy = TRUE)
ggplot(grid_prior_df, aes(x = x, y = y, fill = prior_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +  # A clean, high-contrast color palette
  coord_equal()  


#=============================
# Back to the core data
print(cores)

#==========================================
# Prepare the area-based grid

area_based <- vect(
  input_files[["area_based_map"]],
  layer = "area_based_0_100"
)

# Match the prior's coordinate system.
area_based <- project(
  area_based,
  crs(grid_prior)
)

# Transfer polygon stocks onto the exact prior grid.
# Each cell takes the value of the polygon at its centre.
area_based_grid <- rasterize(
  area_based,
  grid_prior[["prior_mean"]],
  field = "mean_stock_kg_m2",
  background = NA,
  touches = FALSE
)

names(area_based_grid) <- "area_based_stock"

# Retain only cells with valid prior values.
area_based_grid <- mask(
  area_based_grid,
  grid_prior[["prior_mean"]]
)


writeRaster(
  area_based_grid,
  file.path(output_folder, "area_based_stock.tif"),
  overwrite = TRUE
)

#----- Plot -----

area_based_df <- as.data.frame(
  area_based_grid,
  xy = TRUE,
  na.rm = TRUE
)

ggplot(
  area_based_df,
  aes(x = x, y = y, fill = area_based_stock)
) +
  geom_raster() +
  scale_fill_viridis_c(name = "kg C/m²") +
  coord_equal() +
  labs(
    title = paste0(
      "Area-based carbon stock: 0–",
      target_depth_cm, " cm"
    ),
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()



# =============

# Now we shall prep the Inverse distanc weighted analysis

idw_based = rast("outputs/script2/idw_0_100cm.tif")  
  
# Align with the prior grid.
idw_based_grid <- project(
  idw_based,
  grid_prior[["prior_mean"]],
  method = "bilinear"
)

names(idw_based_grid) <- "idw_stock"

# Retain only cells with valid prior values.
idw_based_grid <- mask(
  idw_based_grid,
  grid_prior[["prior_mean"]]
)


writeRaster(
  idw_based_grid,
  file.path(output_folder, "idw_stock.tif"),
  overwrite = TRUE
)

#----- Plot -----

idw_based_df <- as.data.frame(
  idw_based_grid,
  xy = TRUE,
  na.rm = TRUE
)

ggplot(
  idw_based_df,
  aes(x = x, y = y, fill = idw_stock)
) +
  geom_raster() +
  scale_fill_viridis_c(name = "kg C/m²") +
  coord_equal() +
  labs(
    title = "IDW carbon stock: 0–100 cm",
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()

# ==================
# Bayesian Updated


# Now we have our first analysis:

# We'll treat the area-base as the update to the prior carbon map


range_m <- 2000 # A seperate covariate that that controls the weight of the updated core value relative to distant (only used in raw core update)
observation_log_sd <- 0.00 # in cases where we have observer error, this can be used
chunk_size <- 5000L # For larger analyses, this control how many rasters we caluclate in 1 go


# ------- 
cores$prior_mean <- at_cores$prior_mean
cores$prior_sd <- at_cores$prior_sd

# ------------------
# Extract values from the centre of each grid
# Coordinates of the prediction cell centres.
grid_values <- values(
  grid_prior,
  mat = TRUE
)

# Here we only select prior mean and sd with values above 0
eligible_cells <- which(
  is.finite(grid_values[, 1]) &
    grid_values[, 1] > 0 &
    is.finite(grid_values[, 2]) &
    grid_values[, 2] > 0
)
grid_xy <- xyFromCell(
  grid_prior,
  eligible_cells
)

# Coordinates of the observed cores.
core_xy <- crds(
  project(core_points, crs(grid_prior))
)


# All carbon stock models needs to be analysed on log-normal scale

core_log_var <- log1p(
  (cores$prior_sd / cores$prior_mean)^2
)
core_log_mean <- (
  log(cores$prior_mean) -
    core_log_var / 2
)



#----- Prior distributions at the grid cells -----

grid_mean <- grid_values[eligible_cells, 1]
grid_sd <- grid_values[eligible_cells, 2]

grid_log_var <- log1p((grid_sd / grid_mean)^2)

grid_log_mean <- (log(grid_mean) - grid_log_var / 2)


#----- Field measurements on the log scale -----

observed_log <- log(cores$stock)

observation_var <- rep(
  observation_log_sd^2,
  nrow(cores)
)

# --- Bayesian Analysis 1 -------
# Baysian spatial model with Raw core data

# Step 1 - Calculate distance matrix

# HEre we accoutn for

# This uses the ordinary distance formula:
# d_{ij}=\sqrt{(x_i-x_j)^2+(y_i-y_j)^2}

#----- Calculate distances between locations -----

distance_matrix <- function(a, b) {
  
  sqrt(
    outer(a[, 1], b[, 1], "-")^2 +
      outer(a[, 2], b[, 2], "-")^2
  )
  
}

core_distances <- distance_matrix(
  core_xy,
  core_xy
)
core_distances

#----- Spatial correlation between cores -----

core_correlation <- exp(
  -core_distances / range_m
)

round(core_correlation, 2)

#----- Prior covariance between cores -----

core_log_sd <- sqrt(core_log_var)

K <- (
  outer(core_log_sd, core_log_sd) *
    core_correlation
)

round(K, 4)

data.frame(
  core_id = cores$core_id,
  prior_log_variance = core_log_var,
  covariance_diagonal = diag(K)
)

#----- Covariance of observed stocks -----

V <- K + diag(
  observation_var,
  nrow = nrow(cores)
)

round(V, 4)


# ====================
# Make a bayesian inference:
#----- Differences between observations and prior expectations -----
# First calculate the residuals
residuals_1 <- observed_log - core_log_mean
data.frame(
  core_id = cores$core_id,
  observed_log = observed_log,
  prior_log_mean = core_log_mean,
  residual_1= residuals_1
)


#----- Solve for the uncertainty-adjusted residuals -----

R <- chol(V)
alpha <- backsolve(
  R,
  forwardsolve(
    t(R),
    residuals_1
  )
)


# - --------
#----- Select a prediction cell near core 1 -----

distances_to_core_1 <- distance_matrix(
  grid_xy,
  core_xy[1, , drop = FALSE]
)

prediction_row <- which.min(
  distances_to_core_1[, 1]
)

prediction_xy <- grid_xy[
  prediction_row,
  ,
  drop = FALSE
]

prediction_xy
#----- Covariance between this prediction cell and the cores -----

K_new <- (
  outer(
    sqrt(grid_log_var[prediction_row]),
    sqrt(core_log_var)
  ) *
    exp(
      -distance_matrix(prediction_xy, core_xy) /
        range_m
    )
)

signif(K_new, 3)

#----- Posterior log-scale mean -----

mu_post <- (
  grid_log_mean[prediction_row] +
    as.numeric(K_new %*% alpha)
)


#----- Posterior log-scale variance -----

projected <- forwardsolve(
  t(R),
  t(K_new)
)

v_post <- (
  grid_log_var[prediction_row] -
    sum(projected^2)
)

# Only allow negligible negative rounding differences.
stopifnot(
  v_post >= -1e-8 *
    max(1, grid_log_var[prediction_row])
)

v_post <- max(0, v_post)
#----- Inspect the result in carbon units -----

data.frame(
  prior_mean = grid_mean[prediction_row],
  posterior_mean = exp(mu_post + v_post / 2),
  prior_log_variance = grid_log_var[prediction_row],
  posterior_log_variance = v_post
)

#-------
# Now scale for every core and pixel in AOI

#----- Differences between observations and prior expectations -----

residual_log <- observed_log - core_log_mean

# Solve V × alpha = residual_log.
R <- chol(V)

z <- forwardsolve(
  t(R),
  residual_log
)

alpha <- backsolve(R, z)
z
alpha


#----- Bayesian update using field cores -----

core_posterior_log_mean <- numeric(length(eligible_cells))
core_posterior_log_var <- numeric(length(eligible_cells))

for (start in seq(
  1,
  length(eligible_cells),
  by = chunk_size
)) {
  
  rows <- start:min(
    start + chunk_size - 1,
    length(eligible_cells)
  )
  
  # Covariance between these grid cells and the cores.
  K_new <- (
    outer(
      sqrt(grid_log_var[rows]),
      sqrt(core_log_var)
    ) *
      exp(
        -distance_matrix(
          grid_xy[rows, , drop = FALSE],
          core_xy
        ) / range_m
      )
  )
  
  # Updated log-scale mean.
  core_posterior_log_mean[rows] <- (
    grid_log_mean[rows] +
      as.vector(K_new %*% alpha)
  )
  
  # Calculate the reduction in log-scale variance.
  projected <- forwardsolve(
    t(R),
    t(K_new)
  )
  
  updated_var <- (
    grid_log_var[rows] -
      colSums(projected^2)
  )
  
  # Reject substantial negative variances.
  if (any(
    updated_var < -1e-8 * pmax(1, grid_log_var[rows])
  )) {
    stop("Invalid posterior variance; check the inputs.")
  }
  
  core_posterior_log_var[rows] <- pmax(0, updated_var)
}

#----- Results in kg C/m² -----

core_results <- data.frame(
  cell = eligible_cells,
  x = grid_xy[, 1],
  y = grid_xy[, 2],
  prior_mean = grid_mean,
  
  posterior_mean = exp(
    core_posterior_log_mean +
      core_posterior_log_var / 2
  ),
  
  posterior_median = exp(
    core_posterior_log_mean
  ),
  
  lower95 = exp(
    core_posterior_log_mean -
      qnorm(0.975) * sqrt(core_posterior_log_var)
  ),
  
  upper95 = exp(
    core_posterior_log_mean +
      qnorm(0.975) * sqrt(core_posterior_log_var)
  )
)

core_results$change_from_prior <- (
  core_results$posterior_mean -
    core_results$prior_mean
)

head(core_results)


ggplot(
  core_results,
  aes(x = x, y = y, fill = change_from_prior)
) +
  geom_raster() +
  geom_point(
    data = data.frame(
      x = core_xy[, 1],
      y = core_xy[, 2]
    ),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    colour = "black",
    size = 2
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Change\nkg C/m²"
  ) +
  coord_equal() +
  labs(
    title = "Field-core Bayesian update",
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()

#==========================================
#----- Complete the raw-core update --------
#==========================================

#----- Solve for the core information -----

residual_log <- observed_log - core_log_mean

R <- chol(V)

z <- forwardsolve(
  t(R),
  residual_log
)

alpha <- backsolve(R, z)


#----- Update the grid in batches -----

core_posterior_log_mean <- numeric(length(eligible_cells))
core_posterior_log_var <- numeric(length(eligible_cells))

for (start in seq(
  1,
  length(eligible_cells),
  by = chunk_size
)) {
  
  rows <- start:min(
    start + chunk_size - 1,
    length(eligible_cells)
  )
  
  # Covariance between these grid cells and the cores.
  K_new <- (
    outer(
      sqrt(grid_log_var[rows]),
      sqrt(core_log_var)
    ) *
      exp(
        -distance_matrix(
          grid_xy[rows, , drop = FALSE],
          core_xy
        ) / range_m
      )
  )
  
  # Posterior log-scale mean.
  core_posterior_log_mean[rows] <- (
    grid_log_mean[rows] +
      as.vector(K_new %*% alpha)
  )
  
  # Posterior log-scale variance.
  projected <- forwardsolve(
    t(R),
    t(K_new)
  )
  
  updated_var <- (
    grid_log_var[rows] -
      colSums(projected^2)
  )
  
  if (any(
    updated_var < -1e-8 * pmax(1, grid_log_var[rows])
  )) {
    stop("Invalid posterior variance; check the inputs.")
  }
  
  # Remove negligible negative rounding differences.
  core_posterior_log_var[rows] <- pmax(0, updated_var)
}


#----- Results in carbon units -----

core_results <- data.frame(
  cell = eligible_cells,
  x = grid_xy[, 1],
  y = grid_xy[, 2],
  prior_mean = grid_mean,
  prior_sd = grid_sd,
  
  posterior_mean = exp(
    core_posterior_log_mean +
      core_posterior_log_var / 2
  ),
  
  posterior_median = exp(
    core_posterior_log_mean
  ),
  
  posterior_sd = sqrt(
    expm1(core_posterior_log_var) *
      exp(
        2 * core_posterior_log_mean +
          core_posterior_log_var
      )
  ),
  
  lower95 = exp(
    core_posterior_log_mean -
      qnorm(0.975) * sqrt(core_posterior_log_var)
  ),
  
  upper95 = exp(
    core_posterior_log_mean +
      qnorm(0.975) * sqrt(core_posterior_log_var)
  )
)

core_results$interval_width95 <- (
  core_results$upper95 -
    core_results$lower95
)

core_results$change_from_prior <- (
  core_results$posterior_mean -
    core_results$prior_mean
)

core_results$log_variance_reduction <- (
  1 - core_posterior_log_var / grid_log_var
)

# Check that predictions are finite and log-variance decreased.
stopifnot(
  all(is.finite(as.matrix(core_results))),
  all(core_posterior_log_var <= grid_log_var + 1e-8)
)

head(core_results)


#----- Output folder -----

core_folder <- file.path(
  output_folder,
  "raw_cores"
)

dir.create(
  core_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


#----- Build the output raster -----

map_columns <- setdiff(
  names(core_results),
  c("cell", "x", "y")
)

core_maps <- rast(
  grid_prior[[1]],
  nlyrs = length(map_columns)
)

names(core_maps) <- map_columns

map_values <- matrix(
  NA_real_,
  nrow = ncell(core_maps),
  ncol = length(map_columns)
)

map_values[core_results$cell, ] <- as.matrix(
  core_results[, map_columns]
)

values(core_maps) <- map_values


#----- Export -----

writeRaster(
  core_maps,
  file.path(core_folder, "raw_core_maps.tif"),
  overwrite = TRUE
)

write_csv(
  core_results,
  file.path(core_folder, "raw_core_predictions.csv")
)

# Retain log-scale results and settings for later analysis.
saveRDS(
  list(
    range_m = range_m,
    observation_log_sd = observation_log_sd,
    target_depth_cm = target_depth_cm,
    cores = cores,
    core_xy = core_xy,
    observation_var = observation_var,
    posterior_log_mean = core_posterior_log_mean,
    posterior_log_var = core_posterior_log_var,
    predictions = core_results
  ),
  file.path(core_folder, "raw_core_results.rds")
)

#----- Plot the posterior mean -----

core_plot <- ggplot(
  core_results,
  aes(x = x, y = y, fill = posterior_mean)
) +
  geom_raster() +
  geom_point(
    data = data.frame(
      x = core_xy[, 1],
      y = core_xy[, 2]
    ),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_fill_viridis_c(name = "kg C/m²") +
  coord_equal() +
  labs(
    title = paste0(
      "Raw-core update: 0-",
      target_depth_cm,
      " cm"
    ),
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()

print(core_plot)

ggsave(
  file.path(core_folder, "raw_core_mean.png"),
  plot = core_plot,
  width = 8,
  height = 7,
  dpi = 300
)

plot(core_maps[["change_from_prior"]])
plot(core_maps[["posterior_sd"]])
plot(core_maps[["interval_width95"]])

#----- Descriptive map summary -----

core_summary <- data.frame(
  model = c("Prior", "Raw-core update"),
  
  prediction_cells = nrow(core_results),
  
  mean_stock_kg_m2 = c(
    mean(core_results$prior_mean),
    mean(core_results$posterior_mean)
  ),
  
  minimum_stock_kg_m2 = c(
    min(core_results$prior_mean),
    min(core_results$posterior_mean)
  ),
  
  maximum_stock_kg_m2 = c(
    max(core_results$prior_mean),
    max(core_results$posterior_mean)
  ),
  
  mean_pointwise_sd_kg_m2 = c(
    mean(core_results$prior_sd),
    mean(core_results$posterior_sd)
  )
)

print(core_summary)

write_csv(
  core_summary,
  file.path(core_folder, "raw_core_summary.csv")
)



# ==================

# Now lets update the Area-Based Model

#==========================================

area_field_mean <- mean(cores$stock)
area_field_sd <- sd(cores$stock)

stopifnot(
  nrow(cores) >= 3,
  all(is.finite(cores$stock) & cores$stock > 0),
  is.finite(area_field_sd),
  area_field_sd > 0
)

# Confirm that the prepared map uses the current core mean.
area_values <- values(area_based_grid, mat = FALSE)
area_values <- area_values[is.finite(area_values)]

stopifnot(
  length(area_values) > 0,
  all(abs(area_values - area_field_mean) < 1e-5)
)

data.frame(
  field_mean = area_field_mean,
  working_error_sd = area_field_sd
)






#==========================================
#----- Bayesian analysis: Area-based -------
#==========================================

area_folder <- file.path(output_folder, "area_based")

dir.create(
  area_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

#----- Field summaries -----

area_field_mean <- mean(cores$stock)
area_field_sd <- sd(cores$stock)

stopifnot(
  all(is.finite(cores$stock) & cores$stock > 0),
  is.finite(area_field_sd),
  area_field_sd > 0,
  compareGeom(
    grid_prior[[1]],
    area_based_grid,
    stopOnError = FALSE
  )
)

#----- Prepare this analysis's pixel table -----

area_results <- as.data.frame(
  c(grid_prior, area_based_grid),
  cells = TRUE,
  xy = TRUE,
  na.rm = TRUE
) %>%
  filter(
    is.finite(prior_mean), prior_mean > 0,
    is.finite(prior_sd), prior_sd > 0,
    is.finite(area_based_stock), area_based_stock > 0
  )

# Confirm that the map uses the current field mean.
stopifnot(
  nrow(area_results) > 0,
  all(
    abs(area_results$area_based_stock - area_field_mean) <
      1e-5
  )
)

#----- Prior on the log scale -----

area_prior_var <- log1p(
  (area_results$prior_sd / area_results$prior_mean)^2
)

area_prior_mean_log <- (
  log(area_results$prior_mean) -
    area_prior_var / 2
)

#----- Working likelihood -----

# Approximate log-error variance using mapped stock
# as the reference for the kg C/m² error scale.
area_error_var_log <- log1p(
  (area_field_sd / area_results$area_based_stock)^2
)

# Assumption: unbiased observation error on the log scale.
area_weight <- (
  area_prior_var /
    (area_prior_var + area_error_var_log)
)

#----- Posterior on the log scale -----

area_post_mean_log <- (
  area_prior_mean_log +
    area_weight *
    (
      log(area_results$area_based_stock) -
        area_prior_mean_log
    )
)

area_post_var_log <- (
  (1 - area_weight) * area_prior_var
)

#----- Return to carbon units -----

area_results$posterior_log_mean <- area_post_mean_log
area_results$posterior_log_var <- area_post_var_log
area_results$data_weight <- area_weight

area_results$posterior_mean <- exp(
  area_post_mean_log + area_post_var_log / 2
)

area_results$posterior_sd <- sqrt(
  expm1(area_post_var_log) *
    exp(2 * area_post_mean_log + area_post_var_log)
)

area_results$lower95 <- exp(
  area_post_mean_log -
    qnorm(0.975) * sqrt(area_post_var_log)
)

area_results$upper95 <- exp(
  area_post_mean_log +
    qnorm(0.975) * sqrt(area_post_var_log)
)

area_results$change_from_prior <- (
  area_results$posterior_mean -
    area_results$prior_mean
)

area_results$log_variance_reduction <- (
  1 - area_post_var_log / area_prior_var
)

stopifnot(
  all(is.finite(as.matrix(area_results))),
  all(
    area_post_var_log >= 0 &
      area_post_var_log <= area_prior_var
  )
)

#----- Create this analysis's output raster -----

area_map_columns <- c(
  "prior_mean", "prior_sd",
  "posterior_mean", "posterior_sd",
  "lower95", "upper95",
  "change_from_prior", "log_variance_reduction"
)

area_maps <- rast(
  grid_prior[[1]],
  nlyrs = length(area_map_columns)
)

names(area_maps) <- area_map_columns

area_values <- matrix(
  NA_real_,
  nrow = ncell(area_maps),
  ncol = length(area_map_columns)
)

area_values[area_results$cell, ] <- as.matrix(
  area_results[, area_map_columns]
)

values(area_maps) <- area_values

#----- Save -----

writeRaster(
  area_maps,
  file.path(area_folder, "posterior_maps.tif"),
  overwrite = TRUE
)

write_csv(
  area_results,
  file.path(area_folder, "predictions.csv")
)

area_summary <- data.frame(
  method = "Area-based",
  prediction_cells = nrow(area_results),
  mean_prior_stock_kg_m2 = mean(area_results$prior_mean),
  mean_updated_stock_kg_m2 =
    mean(area_results$posterior_mean),
  working_error_sd_kg_m2 = area_field_sd,
  error_source = "Between-core SD",
  model =
    "Approximate pointwise log-Gaussian fusion; no joint covariance",
  likelihood = "Log-unbiased error; plug-in log variance"
)

print(area_summary)

write_csv(
  area_summary,
  file.path(area_folder, "summary.csv")
)

#----- Plot -----

area_plot <- ggplot(
  area_results,
  aes(x, y, fill = posterior_mean)
) +
  geom_raster() +
  scale_fill_viridis_c(name = "kg C/m²") +
  coord_equal() +
  labs(
    title = paste0(
      "Area-based update: 0-",
      target_depth_cm,
      " cm"
    ),
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()

print(area_plot)

ggsave(
  file.path(area_folder, "posterior_mean.png"),
  plot = area_plot,
  width = 8,
  height = 7,
  dpi = 300
)




#==========================================
#----- Bayesian analysis: IDW --------------
#==========================================

idw_folder <- file.path(output_folder, "idw")

dir.create(
  idw_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

#----- Load the settings selected in Script 2 -----

idw_parameters <- read_csv(
  "outputs/script2/idw_best_parameters_by_depth.csv",
  show_col_types = FALSE
)

idw_parameters <- idw_parameters[
  idw_parameters$Variable ==
    paste0("stock_0_", target_depth_cm),
  , drop = FALSE
]

stopifnot(
  nrow(idw_parameters) == 1,
  idw_parameters$N_Cores == nrow(cores),
  all(is.finite(cores$stock) & cores$stock > 0),
  !anyNA(cores$plot_id),
  all(nzchar(cores$plot_id)),
  length(unique(cores$plot_id)) >= 3,
  compareGeom(
    grid_prior[[1]],
    idw_based_grid,
    stopOnError = FALSE
  )
)

#----- Estimate working prediction error -----

idw_points <- sf::st_as_sf(
  project(core_points, crs(grid_prior))
)

idw_heldout <- numeric(nrow(cores))

for (plot_id in unique(cores$plot_id)) {
  
  test <- cores$plot_id == plot_id
  train <- !test
  
  idw_heldout[test] <- gstat::idw(
    stock ~ 1,
    locations = idw_points[train, ],
    newdata = idw_points[test, ],
    idp = idw_parameters$idp,
    nmax = min(idw_parameters$nmax, sum(train)),
    debug.level = 0
  )$var1.pred
}

idw_error_checks <- data.frame(
  core_id = cores$core_id,
  plot_id = cores$plot_id,
  observed = cores$stock,
  predicted = idw_heldout,
  error = cores$stock - idw_heldout
)

idw_error_sd <- sqrt(
  mean(idw_error_checks$error^2)
)

stopifnot(
  is.finite(idw_error_sd),
  idw_error_sd > 0
)

#----- Prepare this analysis's pixel table -----

idw_results <- as.data.frame(
  c(grid_prior, idw_based_grid),
  cells = TRUE,
  xy = TRUE,
  na.rm = TRUE
) %>%
  filter(
    is.finite(prior_mean), prior_mean > 0,
    is.finite(prior_sd), prior_sd > 0,
    is.finite(idw_stock), idw_stock > 0
  )

stopifnot(nrow(idw_results) > 0)

#----- Prior on the log scale -----

idw_prior_var <- log1p(
  (idw_results$prior_sd / idw_results$prior_mean)^2
)

idw_prior_mean_log <- (
  log(idw_results$prior_mean) -
    idw_prior_var / 2
)

#----- Working likelihood -----

idw_error_var_log <- log1p(
  (idw_error_sd / idw_results$idw_stock)^2
)

# Assumption: unbiased observation error on the log scale.
idw_weight <- (
  idw_prior_var /
    (idw_prior_var + idw_error_var_log)
)

#----- Posterior on the log scale -----

idw_post_mean_log <- (
  idw_prior_mean_log +
    idw_weight *
    (
      log(idw_results$idw_stock) -
        idw_prior_mean_log
    )
)

idw_post_var_log <- (
  (1 - idw_weight) * idw_prior_var
)

#----- Return to carbon units -----

idw_results$posterior_log_mean <- idw_post_mean_log
idw_results$posterior_log_var <- idw_post_var_log
idw_results$data_weight <- idw_weight

idw_results$posterior_mean <- exp(
  idw_post_mean_log + idw_post_var_log / 2
)

idw_results$posterior_sd <- sqrt(
  expm1(idw_post_var_log) *
    exp(2 * idw_post_mean_log + idw_post_var_log)
)

idw_results$lower95 <- exp(
  idw_post_mean_log -
    qnorm(0.975) * sqrt(idw_post_var_log)
)

idw_results$upper95 <- exp(
  idw_post_mean_log +
    qnorm(0.975) * sqrt(idw_post_var_log)
)

idw_results$change_from_prior <- (
  idw_results$posterior_mean -
    idw_results$prior_mean
)

idw_results$log_variance_reduction <- (
  1 - idw_post_var_log / idw_prior_var
)

stopifnot(
  all(is.finite(as.matrix(idw_results))),
  all(
    idw_post_var_log >= 0 &
      idw_post_var_log <= idw_prior_var
  )
)

#----- Create this analysis's output raster -----

idw_map_columns <- c(
  "prior_mean", "prior_sd",
  "posterior_mean", "posterior_sd",
  "lower95", "upper95",
  "change_from_prior", "log_variance_reduction"
)

idw_maps <- rast(
  grid_prior[[1]],
  nlyrs = length(idw_map_columns)
)

names(idw_maps) <- idw_map_columns

idw_values <- matrix(
  NA_real_,
  nrow = ncell(idw_maps),
  ncol = length(idw_map_columns)
)

idw_values[idw_results$cell, ] <- as.matrix(
  idw_results[, idw_map_columns]
)

values(idw_maps) <- idw_values

#----- Save -----

writeRaster(
  idw_maps,
  file.path(idw_folder, "posterior_maps.tif"),
  overwrite = TRUE
)

write_csv(
  idw_results,
  file.path(idw_folder, "predictions.csv")
)

write_csv(
  idw_error_checks,
  file.path(idw_folder, "error_checks.csv")
)

idw_summary <- data.frame(
  method = "IDW",
  prediction_cells = nrow(idw_results),
  mean_prior_stock_kg_m2 = mean(idw_results$prior_mean),
  mean_updated_stock_kg_m2 =
    mean(idw_results$posterior_mean),
  working_error_sd_kg_m2 = idw_error_sd,
  error_source =
    "Held-out plot RMSE; IDW settings fixed, not nested validation",
  model =
    "Approximate pointwise log-Gaussian fusion; no joint covariance",
  likelihood = "Log-unbiased error; plug-in log variance"
)

print(idw_summary)

write_csv(
  idw_summary,
  file.path(idw_folder, "summary.csv")
)

#----- Plot -----

idw_plot <- ggplot(
  idw_results,
  aes(x, y, fill = posterior_mean)
) +
  geom_raster() +
  scale_fill_viridis_c(name = "kg C/m²") +
  coord_equal() +
  labs(
    title = paste0(
      "IDW update: 0-",
      target_depth_cm,
      " cm"
    ),
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()

print(idw_plot)

ggsave(
  file.path(idw_folder, "posterior_mean.png"),
  plot = idw_plot,
  width = 8,
  height = 7,
  dpi = 300
)









